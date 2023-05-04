library(tidyverse)
library(foreach)
library(networkD3)

source("./scripts/function_source_for_network_visualization.R")
source("./scripts/function_source_for_network_construction.R")

# The called CellTag V1, V2 and V3 will be used to construct Lineage network graph.
# Input Celltag data should be Nx3 matrix. Each row represents cell name (cell barcode) and each column represents Celltag name.
# See the example (./input_data_for_network_construction/)

# load celltag table

mef.clones <- read_csv(file = "hf1.d15.v1.clones.csv")

d3.clones <- read_csv(file = "hf1.d15.v2.clones.csv")
  
d13.clones <- read_csv(file = "hf1.d15.v3.clones.csv")


colnames(mef.clones)[1] <- "CellTagV1"

colnames(d3.clones)[1] <- "CellTagV2"

colnames(d13.clones)[1] <- "CellTagV3"

clone.cells <- unique(c(mef.clones$cell.barcode, d3.clones$cell.barcode, d13.clones$cell.barcode))

celltag_data <- data.frame(clone.cells, row.names = clone.cells)

celltag_data$CellTagV1 <- NA

celltag_data$CellTagV2 <- NA

celltag_data$CellTagV3 <- NA

celltag_data[mef.clones$cell.barcode, "CellTagV1"] <- mef.clones$CellTagV1

celltag_data[d3.clones$cell.barcode, "CellTagV2"] <- d3.clones$CellTagV2

celltag_data[d13.clones$cell.barcode, "CellTagV3"] <- d13.clones$CellTagV3

celltag_data <- celltag_data[, -1]

row.names(celltag_data) <- paste0(rownames(celltag_data), "-1")


colnames(celltag_data) <- c("CellTagV1", "CellTagV2", "CellTagV3")

linkList <- convertCellTagMatrix2LinkList(celltag_data)


Nodes <- getNodesfromLinkList(linkList)


additional_data <- data.frame(sample(1:10, size = length(row.names(celltag_data)), replace = TRUE), row.names = row.names(celltag_data))
Nodes <- addData2Nodes(Nodes, additional_data)
colnames(Nodes)[4] <- "Cluster"
head(Nodes)


drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", linkList = linkList, Nodes = Nodes )

bar.data <- celltag_data

bar.data$Cell.BC <- row.names(bar.data)

bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)


ggplot(data = bar.data) + geom_bar(mapping = aes(x = CellTag, fill = factor(Clone)), position = "fill", show.legend = FALSE) + scale_y_continuous(labels = scales::percent_format())
