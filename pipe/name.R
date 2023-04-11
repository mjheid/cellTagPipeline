library(CellTagR)

# Recount and generate collapsed matrix
bam.test.obj <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj, collapsed.rslt.file = "~/Desktop/collapsing_rslt.txt")
# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)

collapsed.rslt.dir <- "~/Desktop/star_collapse"
# Recount and generate collapsed matrix
bam.test.obj <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)

# Calling binarization
bam.test.obj <- SingleCellDataBinarization(bam.test.obj, 2)

MetricPlots(bam.test.obj)

# Read the RDS file and get the object
dt.mtx.whitelist.path <- system.file("extdata", "v1_whitelist.csv", package = "CellTagR")
bam.test.obj <- SingleCellDataWhitelist(bam.test.obj, dt.mtx.whitelist.path)

MetricPlots(bam.test.obj)

bam.test.obj <- MetricBasedFiltering(bam.test.obj, 20, comparison = "less")
bam.test.obj <- MetricBasedFiltering(bam.test.obj, 2, comparison = "greater")

MetricPlots(bam.test.obj)

bam.test.obj <- JaccardAnalysis(bam.test.obj)

# Call clones
bam.test.obj <- CloneCalling(celltag.obj = bam.test.obj, correlation.cutoff=0.7)
# Check them out!!
bam.test.obj@clone.composition[["v1"]]
bam.test.obj@clone.size.info[["v1"]]

# Read the RDS file and get the object
dt.mtx.path <- system.file("extdata", "bam_v123_obj.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)

bam.test.obj <- convertCellTagMatrix2LinkList(bam.test.obj)

bam.test.obj <- getNodesfromLinkList(bam.test.obj)

# Simulate some additional data
additional_data <- data.frame(sample(1:10, size = length(rownames(bam.test.obj@celltag.aggr.final)), replace = TRUE), row.names = rownames(bam.test.obj@celltag.aggr.final))
colnames(additional_data) <- "Cluster"
# Add the data to the object
bam.test.obj <- addData2Nodes(bam.test.obj, additional_data)

# Network Visualization
bam.test.obj <- drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", celltag.obj = bam.test.obj)
bam.test.obj@network

saveNetwork(bam.test.obj@network, "~/Desktop/presentation/Demo/hf1.d15.network.construction.html")


# Get the data for ploting
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)

bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)

# Using ggplot to plot
ggplot(data = bar.data) + 
  geom_bar(mapping = aes(x = CellTag, fill = factor(Clone)), position = "fill", show.legend = FALSE) + 
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()
