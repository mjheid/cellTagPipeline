suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))

option_list = list(
    make_option(c("--visualize"), action="store_true", default=FALSE, 
              help="Visualize filtering"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--start_node"), type="character", default="CellTagV1_2", 
              help="Start node to build Graph from.", metavar="character"),
    make_option(c("--save_progress_name"), type="character", default="test", 
              help="Name appended to outputs to be saved", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

DIR <- opt$out

# Read the RDS file and get the object
bam.test.obj <- readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

bam.test.obj <- convertCellTagMatrix2LinkList(bam.test.obj)

bam.test.obj <- getNodesfromLinkList(bam.test.obj)

# Simulate some additional data
# TODO: Give cells a cluster number, so we can visualize development??? How? Based on what? Found cellTags?
#       Cell identity? Cell barcode? UMI? How do i traceit back??? Why do they not have a fucking example
#       that showcases the actual usage of their libary, do they not want others to use their work???
additional_data <- data.frame(sample(1:10, size = length(rownames(bam.test.obj@celltag.aggr.final)), replace = TRUE), row.names = rownames(bam.test.obj@celltag.aggr.final))
colnames(additional_data) <- "Cluster"
# Add the data to the object
bam.test.obj <- addData2Nodes(bam.test.obj, additional_data)

# Network Visualization
bam.test.obj <- drawSubnet(tag = opt$start_node, overlay = "Cluster", celltag.obj = bam.test.obj)
saveNetwork(bam.test.obj@network, paste0(DIR, opt$save_progress_name, opt$start_node, "network.construction.html"))

# Get the data for ploting
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)

bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)

# Using ggplot to plot
pdf(paste0(DIR, opt$save_progress_name, "_clone_barchart.pdf"))
ggplot(data = bar.data) + 
  geom_bar(mapping = aes(x = CellTag, fill = factor(Clone)), position = "fill", show.legend = FALSE) + 
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()
dev.off()