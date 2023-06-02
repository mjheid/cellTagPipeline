suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))

option_list = list(
    make_option(c("--visualize"), action="store_true", default=FALSE, 
              help="Visualize filtering"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--start_node"), type="character", default="BiggestNode", 
              help="Start node to build Graph from.", metavar="character"),
    make_option(c("--save_progress_name"), type="character", default="insert", 
              help="Name appended to outputs to be saved", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


tryCatch(
{


DIR <- opt$out

# Read the RDS file and get the object
bam.test.obj <- readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

if (file.exists(paste0(DIR, "v1", opt$save_progress_name, ".RDS"))) {
    bam.v1.obj <- readRDS(paste0(DIR, "v1", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v1 = bam.v1.obj@clone.composition$v1
}
if (file.exists(paste0(DIR, "v2", opt$save_progress_name, ".RDS"))) {
    bam.v2.obj <- readRDS(paste0(DIR, "v2", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v2 = bam.v2.obj@clone.composition$v2
}
if (file.exists(paste0(DIR, "v3", opt$save_progress_name, ".RDS"))) {
    bam.v3.obj <- readRDS(paste0(DIR, "v3", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v3 = bam.v3.obj@clone.composition$v3
}

bam.test.obj <- convertCellTagMatrix2LinkList(bam.test.obj)

bam.test.obj <- getNodesfromLinkList(bam.test.obj)

saveRDS(bam.test.obj, paste0(DIR, opt$save_progress_name, ".RDS"))


# Network Visualization
if (opt$start_node == "BiggestNode") {
    # Grab all unique source noes and how often they appear, order them, and use the first one to drawSubnet
    tmp = table(bam.test.obj@network.link.list$source)
    start_nodes = names(tmp)[order(-tmp)]
    bam.test.obj <- drawSubnet(tag = start_nodes[1], overlay = "tag", celltag.obj = bam.test.obj)
} else {
    bam.test.obj <- drawSubnet(tag = opt$start_node, overlay = "tag", celltag.obj = bam.test.obj)
}
saveNetwork(bam.test.obj@network, paste0(DIR, opt$save_progress_name, opt$start_node, "network.construction.html"))

# Get the data for ploting
library(ggplot2)
library(dplyr)
library(tidyr)

bar.data <- bam.test.obj@celltag.aggr.final

# Reshape the data to long format
bar.data <- bar.data %>%
  rownames_to_column(var = "uniqueID") %>%
  gather(key = "CellTag", value = "Clone", -uniqueID, na.rm = FALSE)

# Create a mapping dataframe for the lines
bar.data$CellTag <- factor(bar.data$CellTag)
bar.data$Clone <- factor(bar.data$Clone)

# Using ggplot to plot
# Plot the bar chart
p <- ggplot(data = bar.data) + 
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
