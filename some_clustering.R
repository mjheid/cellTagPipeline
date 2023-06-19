suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
bam.test.obj <- readRDS(paste0("data/out/big_b1_l1_h20", ".RDS"))
bam.test.obj
q()
suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
dt.mtx.path <- system.file("extdata", "bam_v123_obj.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)
'
)
dt.mtx.path <- system.file("extdata", "bam_v123_obj.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)
bam.test.obj@version
bam.test.obj@celltag.version
bam.test.obj@bam.parse.rslt
x= CellTagObject(object.name = "test", fastq.bam.directory = "data/bam/possorted_genome_bam.filtered.bam")
x=CellTagExtraction(x, celltag.version = "v1")
x@bam.parse.rslt
bam.test.obj
bam.test.obj@network
bam.test.obj <- drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", celltag.obj = bam.test.obj)
bam.test.obj <- convertCellTagMatrix2LinkList(bam.test.obj)
bam.test.obj <- getNodesfromLinkList(bam.test.obj)
additional_data <- data.frame(sample(1:10, size = length(rownames(bam.test.obj@celltag.aggr.final)), replace = TRUE), row.names = rownames(bam.test.obj@celltag.aggr.final))
colnames(additional_data) <- "Cluster"
bam.test.obj <- addData2Nodes(bam.test.obj, additional_data)
bam.test.obj <- drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", celltag.obj = bam.test.obj)
bam.test.obj@network
pdf("net.pdf")
bam.test.obj@network
saveNetwork(bam.test.obj@network, "net.html")
bam.test.obj@celltag.aggr.final
dim(bam.test.obj@celltag.aggr.final)
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)
bar.data
dim(bar.data)
gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)
1123*3
bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)
dim(bar.data)
additional_data
additional_datadt.mtx.path <- system.file("extdata", "bam_v123_obj.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)
dt.mtx.path <- system.file("extdata", "bam_v123_obj.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)
bam.test.obj <- convertCellTagMatrix2LinkList(bam.test.obj)
bam.test.obj <- getNodesfromLinkList(bam.test.obj)
bam.test.obj@clones
bam.test.obj@clone.composition
ls
bam.test.obj@clone.size.info
bam.test.obj@network.link.list
bam.test.obj@nodes
bam.test.obj
additional_data <- data.frame(sample(1:10, size = length(rownames(bam.test.obj@celltag.aggr.final)), replace = TRUE), row.names = rownames(bam.test.obj@celltag.aggr.final))
colnames(additional_data) <- "Cluster"
additional_data
bam.test.obj@network.link.list
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)
bar.data
bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)
bar.data
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)
bar.data@CellTag
bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)
bar.data
bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)
bar.data
bam.test.obj@clone.composition
bam.test.obj@clone.composition$v1
bam.test.obj@clone.composition$v1[rownames(bam.test.obj@clone.composition$v1), "clone.id"]
bam.test.obj@celltag.aggr.final
bam.test.obj@metric.filtered.count
bam.test.obj@celltag.aggr.final
bar.data <- bam.test.obj@celltag.aggr.final
bar.data$Cell.BC <- rownames(bar.data)
bar.data
q()
suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
kek1  <- readRDS("data/out/v1insert.RDS")
kek2  <- readRDS("data/out/v2insert.RDS")
kek3  <- readRDS("data/out/v3insert.RDS")
kek  <- readRDS("data/out/insert.RDS")
kek
kek1
kek@celltag.version
kek1@celltag.version
kek2@celltag.version
kek3@celltag.version
kek3@clone.composition
kek3@clone.composition$v3
kek$clone.composition$v1 = kek1$clone.composition$v1
kek1@clone.composition
kek1@clone.composition$v1
kek@clone.composition$v1 = kek1@clone.composition$v1
kek@clone.composition$v2 = kek1@clone.composition$v2
kek@clone.composition
kek@clone.composition$v2 = kek2@clone.composition$v2
kek@clone.composition
kek <- convertCellTagMatrix2LinkList(kek)
kek <- getNodesfromLinkList(kek)
additional_data <- data.frame(sample(1:10, size = length(rownames(bam.test.obj@celltag.aggr.final)), replace = TRUE), row.names = rownames(bam.test.obj@celltag.aggr.final))
colnames(additional_data) <- "Cluster"
kek <- addData2Nodes(kek, additional_data)
bam.test.obj <- drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", celltag.obj = kek)
saveNetwork(kek@network, "test_kek.html")
kek<- drawSubnet(tag = "CellTagV1_2", overlay = "Cluster", celltag.obj = kek)
saveNetwork(kek@network, "test_kek.html")
kek@celltag.aggr.final
kek@nodes
kek<- drawSubnet(tag = "CellTagV2_3", overlay = "Cluster", celltag.obj = kek)
saveNetwork(kek@network, "test_kek.html")
quit()
suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
bam.test.obj <- readRDS("data/out/insert.RDS")
bam.test.obj@celltag.version
bam.test.obj@metric.filtered.count$v3
bam.test.obj@metric.filtered.count
bam.test.obj@clone.composition
bam.test.obj@clone.composition$v1
bam.test.obj@clone.composition$v3
bam.test.obj@clone.composition$v3$cell.barcode
q()
rm(list = ls()),
rm(list = ls())
q()
suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
There were 20 warnings (use warnings() to see them)
x = readRDS("data/out/anika/anika.RDS")
pbmc = readRDS("data/out/anika/anika_reduction.RDS")
subset = readRDS("data/out/anika/anika_subset_reduction.RDS")
filtered_bc =rownames(x@jaccard.mtx)
Error: unexpected symbol in:
"suppressMessages(library("Seurat"))
x = readRDS("data/out/anika/anika.RDS")
pbmc = readRDS("data/out/anika/anika_reduction.RDS")
subset = readRDS("data/out/anika/anika_subset_reduction.RDS")
filtered_bc =rownames(x@jaccard.mtx)
Idents(pbmc) = "Without quality CellTag"

x = readRDS("data/out/anika/anika.RDS")
pbmc = readRDS("data/out/anika/anika_reduction.RDS")
subset = readRDS("data/out/anika/anika_subset_reduction.RDS")
filtered_bc =rownames(x@jaccard.mtx)
Idents(pbmc) = "Without quality CellTag"
cluster9_bc = x@clone.composition$v1$cell.barcode[x@clone.composition$v1$clone.id == 9]
cluster13_bc = x@clone.composition$v1$cell.barcode[x@clone.composition$v1$clone.id == 13]
Idents(pbmc) = "Without quality CellTag"
pbmc = SetIdent(object = pbmc, cells = cluster9_bc, value = "clone9")
pbmc = SetIdent(object = pbmc, cells = cluster13_bc, value = "clone13")
cluster5.markers <- FindMarkers(pbmc, ident.1 = "clone9", ident.2 = "clone13", min.pct = 0.25)
cluster5.markers
cluster5.markershead(cluster5.markers, n = 5)
head(cluster5.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = "clone9", ident.2 = "clone13", min.pct = 0.25, logfc.threshold = 0.25)
head(cluster5.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = "clone9", ident.2 = "clone13", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
head(cluster5.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = "clone13", ident.2 = "clone9", min.pct = 0.25, logfc.threshold = 0.25)
head(cluster5.markers, n = 5)
c9_13__bc = append(cluster9_bc, cluster13_bc)
c9_13 = subset(pbmc, cells = c9_13)
c9_13__bc = append(cluster9_bc, cluster13_bc)
c9_13 = subset(pbmc, cells = c9_13__bc)
pbmc = c9_13
pbmc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 26)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 1.5)
pbmc <- FindClusters(pbmc, resolution = 1)
length(Indents(pbmc) == 0)
length(Idents(pbmc) == 0)
length(Idents(pbmc) == 1)
Idents(pbmc)
pdf("9_13_clustered.pdf")
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
dev.off()
c9_13__bc = append(cluster9_bc, cluster13_bc)
c9_13 = subset(pbmc, cells = cluster9_bc)
pbmc = c9_13
pbmc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 26)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- FindClusters(pbmc, resolution = 1.2)
pbmc <- FindClusters(pbmc, resolution = 1.1)
pbmc <- FindClusters(pbmc, resolution = 1.05)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10, n_neighbors=3)
pbmc <- FindClusters(pbmc, resolution = 1.05)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:5)
pbmc <- RunUMAP(pbmc, dims = 1:5, verbose=F)
pbmc <- RunTSNE(pbmc)
pbmc <- RunTSNE(pbmc, dims=1:10)
pbmc <- RunTSNE(pbmc, dims=1:5)
pbmc <- RunTSNE(pbmc, dims=1:3)
pbmc <- RunTSNE(pbmc, dims=1:2)
pbmc <- RunTSNE(pbmc, dims=1)
pbmc <- RunTSNE(pbmc, dims=0)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
pbmc.markers
savehistory(file = "some_clustering.R")
