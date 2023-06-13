> suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
There were 20 warnings (use warnings() to see them)
Warning message:
package ‘Seurat’ was built under R version 4.2.2 
> bam.test.obj = readRDS("Patrick_Anika_22-10/scRNAseq/07_extCellTag_Patrick_Anika/output/SI-TT-H4/outs/Anika_out/bam.test.obj.RDS")
> bam.test.obj
Object name:  bam.cell.tag.obj 
Raw CellTag Counts =  2102 
Raw Number of Cells with CellTag =  7000 
Collapsed CellTag Counts =  2083 
Whitelisted CellTag Counts =  915 
Whitelisted Number of Cells with CellTag =  1919 
> x = readRDS("data/out/anika/anika.RDS")
> x
Object name:  bam.cell.tag.obj 
Raw CellTag Counts =  1970 
Raw Number of Cells with CellTag =  7000 
Collapsed CellTag Counts =  1955 
Whitelisted CellTag Counts =  792 
Whitelisted Number of Cells with CellTag =  1804 
> length(rownames(bam.test.obj@jaccard.mtx))
[1] 1016
> length(rownames(bam.test.obj@metric.filtered.count))
[1] 1919
> curr.mtx <- slot(bam.test.obj, "raw.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1919
> curr.mtx <- slot(bam.test.obj, "whitelisted.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1016
> curr.mtx <- slot(bam.test.obj, "metric.filtered.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1016
> bam.test.obj@clone.composition
$v1
     clone.id       cell.barcode
  1:        1 ACGTCCTAGAGATTCA-1
  2:        1 AAAGTGAGTTTCACAG-1
  3:        2 CATACAGCATGGTGGA-1
  4:        2 CTCACTGTCCGCTAGG-1
  5:        2 AACAAAGCAAAGTGTA-1
 ---                            
270:      100 TCGGGTGTCGGAGCAA-1
271:      101 TTAGGCAAGTCAGGGT-1
272:      101 TCTATACAGGGTGAGG-1
273:      102 TTGGTTTAGCAGTAAT-1
274:      102 TGCACGGTCATAAGGA-1

> length(rownames(bam.test.obj@clone.composition$v1$cell.barcode))
[1] 0
> length(bam.test.obj@clone.composition$v1$cell.barcode)
[1] 274
> curr.mtx <- slot(bam.test.obj, "whitelisted.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1919
> x
Object name:  bam.cell.tag.obj 
Raw CellTag Counts =  1970 
Raw Number of Cells with CellTag =  7000 
Collapsed CellTag Counts =  1955 
Whitelisted CellTag Counts =  792 
Whitelisted Number of Cells with CellTag =  1804 
> curr.mtx <- slot(x, "whitelisted.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1804
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 879
> curr.mtx <- slot(x, "metric.filtered.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 879
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub)),]                                                           
length(rownames(full.mtx.sub))
[1] 879
> savehistory(file = "session_history.R")
> curr.mtx <- slot(x, "collapsed.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub)),]                                                           
length(rownames(full.mtx.sub))
[1] 1804
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1804

> curr.mtx <- slot(x, "collapsed.count")
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub)),]                                                           
length(rownames(full.mtx.sub))
[1] 1804
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 1804
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub <= 1) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 35
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub <= 2) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 3
> curr.version <- bam.test.obj@curr.version
curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub <= 3) == ncol(curr.mtx.sub)),]
length(rownames(full.mtx.sub))
[1] 0
> 

