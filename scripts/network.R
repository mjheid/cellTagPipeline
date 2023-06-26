suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))

# Command line options
option_list = list(
    make_option(c("--visualize"), action="store_true", default=FALSE, 
              help="Visualize filtering"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--start_node"), type="character", default="BiggestNode", 
              help="Start node to build Graph from.", metavar="character"),
    make_option(c("--save_progress_name"), type="character", default="test", 
              help="Name appended to outputs to be saved", metavar="character")
)

# read arguments and put into data.frame op
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Catch errors occuring. When catched, print error and exit programm.
tryCatch(
{

# Output directory. Assumes that data from previous steps is in here
DIR = opt$out

# Read complete celltag obj of previous steps.
bam.test.obj = readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

# If the pipeline has been run with library v1, get clone composition of v1
if (file.exists(paste0(DIR, "v1", opt$save_progress_name, ".RDS"))) {
    bam.v1.obj = readRDS(paste0(DIR, "v1", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v1 = bam.v1.obj@clone.composition$v1
}
# If the pipeline has been run with library v2, get clone composition of v2
if (file.exists(paste0(DIR, "v2", opt$save_progress_name, ".RDS"))) {
    bam.v2.obj = readRDS(paste0(DIR, "v2", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v2 = bam.v2.obj@clone.composition$v2
}
# If the pipeline has been run with library v3, get clone composition of v3
if (file.exists(paste0(DIR, "v3", opt$save_progress_name, ".RDS"))) {
    bam.v3.obj = readRDS(paste0(DIR, "v3", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v3 = bam.v3.obj@clone.composition$v3
}

# CellTag matrix to Linked List. Uses clone composition.
bam.test.obj = convertCellTagMatrix2LinkList(bam.test.obj)

# Create Network nodes from linked list
bam.test.obj = getNodesfromLinkList(bam.test.obj)

# Save cellTag object
saveRDS(bam.test.obj, paste0(DIR, opt$save_progress_name, ".RDS"))


# Do Network Visualization, if start node is set to BiggestNode start with that node
if (opt$start_node == "BiggestNode") {
    # Grab all unique source noes and how often they appear, order them, and use the first one to drawSubnet
    tmp = table(bam.test.obj@network.link.list$source)
    start_nodes = names(tmp)[order(-tmp)]
    bam.test.obj = drawSubnet(tag = start_nodes[1], overlay = "tag", celltag.obj = bam.test.obj)
} 
# Do Network Visualization with User Start -node
else {
    bam.test.obj = drawSubnet(tag = opt$start_node, overlay = "tag", celltag.obj = bam.test.obj)
}
# Save created network
saveNetwork(bam.test.obj@network, paste0(DIR, opt$save_progress_name, opt$start_node, "network.construction.html"))

# Get the data for ploting celltags per version
bar.data = bam.test.obj@celltag.aggr.final

# Reshape the data to long format
bar.data = bar.data %>%
  rownames_to_column(var = "uniqueID") %>%
  gather(key = "CellTag", value = "Clone", -uniqueID, na.rm = FALSE)

# Create a mapping dataframe. TODO: Pretty sure not needed anymore.
bar.data$CellTag = factor(bar.data$CellTag)
bar.data$Clone = factor(bar.data$Clone)

# Plot the bar chart
p = ggplot(data = bar.data) + 
  geom_bar(mapping = aes(x = CellTag, fill = factor(Clone)), position = "fill", show.legend = FALSE) + 
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()

ggsave(paste0(DIR, opt$save_progress_name, "_clone_barchart.pdf"), p)

},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)
