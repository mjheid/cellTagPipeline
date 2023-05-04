library(igraph)
library(proxy)
library(corrplot)
library(data.table)

source("./scripts/CellTagCloneCalling_Function.R")

#load cell matrices and binarize, name columns after corresponding cellTag

mef.mat <- as.data.frame(readRDS("./hf1.d15.v1.celltag.matrix.Rds"))

d3.mat <- as.data.frame(readRDS("./hf1.d15.v2.celltag.matrix.Rds"))

d13.mat <- as.data.frame(readRDS("./hf1.d15.v3.celltag.matrix.Rds"))

rownames(mef.mat) <- mef.mat$Cell.BC

rownames(d3.mat) <- d3.mat$Cell.BC

rownames(d13.mat) <- d13.mat$Cell.BC

mef.mat <- mef.mat[,-1]

d3.mat <- d3.mat[,-1]

d13.mat <- d13.mat[,-1]

mef.bin <- SingleCellDataBinarization(celltag.dat = mef.mat, 2)

d3.bin <- SingleCellDataBinarization(celltag.dat = d3.mat, 2)

d13.bin <- SingleCellDataBinarization(celltag.dat = d13.mat, 2)

#Fillter cells based on whitelist
mef.filt <- SingleCellDataWhitelist(celltag.dat = mef.bin, whitels.cell.tag.file = "whitelist/V1.CellTag.Whitelist.csv")

d3.filt <- SingleCellDataWhitelist(celltag.dat = d3.bin, whitels.cell.tag.file = "whitelist/V2.CellTag.Whitelist.csv")

d13.filt <- SingleCellDataWhitelist(celltag.dat = d13.bin, whitels.cell.tag.file = "whitelist/V3.CellTag.Whitelist.csv")


mef.filt <- MetricBasedFiltering(whitelisted.celltag.data = mef.filt, cutoff = 20, comparison = "less")

d3.filt <- MetricBasedFiltering(whitelisted.celltag.data = d3.filt, cutoff = 20, comparison = "less")

d13.filt <- MetricBasedFiltering(whitelisted.celltag.data = d13.filt, cutoff = 20, comparison = "less")


mef.filt <- MetricBasedFiltering(whitelisted.celltag.data = mef.filt, cutoff = 2, comparison = "greater")

d3.filt <- MetricBasedFiltering(whitelisted.celltag.data = d3.filt, cutoff = 2, comparison = "greater")

d13.filt <- MetricBasedFiltering(whitelisted.celltag.data = d13.filt, cutoff = 2, comparison = "greater")

# Clone calling and saving outputs
mef.sim <- JaccardAnalysis(whitelisted.celltag.data = mef.filt, plot.corr = FALSE, id = "mef")

d3.sim <- JaccardAnalysis(whitelisted.celltag.data = d3.filt, plot.corr = FALSE, id = "d3")

d13.sim <- JaccardAnalysis(whitelisted.celltag.data = d13.filt, plot.corr = FALSE, id = "d13")


mef.clones <- CloneCalling(Jaccard.Matrix = mef.sim, output.dir = "./", output.filename = "hf1.d15.v1.clones.csv", correlation.cutoff = 0.7)

d3.clones <- CloneCalling(Jaccard.Matrix = d3.sim, output.dir = "./", output.filename = "hf1.d15.v2.clones.csv", correlation.cutoff = 0.7)

d13.clones <- CloneCalling(Jaccard.Matrix = d13.sim, output.dir = "./", output.filename = "hf1.d15.v3.clones.csv", correlation.cutoff = 0.7)