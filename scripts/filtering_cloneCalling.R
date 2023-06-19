suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
 
option_list = list(
    make_option(c("--visualize"), action="store_true", default=FALSE, 
              help="Visualize filtering"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--collapsing_name"), type="character", default="collapsing.txt", 
              help="Name of collapsing file", metavar="character"),
    make_option(c("--save_progress_name"), type="character", default="test", 
              help="Name appended to outputs to be saved", metavar="character"),
    make_option(c("--whitelist_version"), type="character", default="v1", 
                help="Version of whitelist", metavar="character"),
    make_option(c("--whitelist_path"), type="character", default="data/whitelist/V1.CellTag.Whitelist.csv", 
              help="Path of whitelist", metavar="character"),
    make_option(c("--high_filter"),  default=20, 
              help="Upper bound for filtering celltags", metavar="int"),
    make_option(c("--low_filter"),  default=1, 
              help="Lower bound for filtering celltags", metavar="int"),
    make_option(c("--tagged"),  default=1, 
              help="When binarizing, the number of cellTags needed to consider cell as tagged", metavar="int")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tryCatch(
{

print_cellTag = function(object) {
            cat("Object name: ", object@obj.name, "\n")
            cat("Library version: ", object@curr.version, "\n")
            
            curr.mtx <- slot(object, "raw.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols <- !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub <- curr.mtx.sub[, valid_cols]
            cat("Unique Raw CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Raw CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            curr.mtx <- slot(object, "raw.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Raw Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            curr.mtx <- slot(object, "collapsed.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols <- !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub <- curr.mtx.sub[, valid_cols]
            cat("Unique Collapsed CellTag Counts = ", ncol(full.mtx.sub), "\n")
            cat("Collapsed CellTag Counts = ", sum(full.mtx.sub, na.rm = TRUE), "\n")
            
            curr.mtx <- slot(object, "whitelisted.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols <- !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub <- curr.mtx.sub[, valid_cols]
            cat("Unique Whitelisted CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Whitelisted CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            curr.mtx <- slot(object, "whitelisted.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Whitelisted Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            curr.mtx <- slot(object, "metric.filtered.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols <- !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub <- curr.mtx.sub[, valid_cols]
            cat("Unique Filtered CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Filtered CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            curr.mtx <- slot(object, "metric.filtered.count")
            curr.version <- object@curr.version
            curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Filtered Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            cat("Number of identified clones = ", length(object@clone.composition$v1$cell.barcode), "\n")
            cat("Number of unique clones = ", length(unique(object@clone.composition$v1$clone.id)), "\n")
            print("Number of clones per unique clone = ")
            print(table(object@clone.composition$v1$clone.id))
}


DIR <- opt$out

bam.test.obj <- readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

bam.test.obj <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj, collapsed.rslt.file = paste0(DIR, opt$whitelist_version, opt$collapsing_name))

bam.test.obj <- SingleCellDataBinarization(bam.test.obj, opt$tagged)


if (opt$visualize) {
pdf(paste0(DIR, opt$whitelist_version, opt$save_progress_name, "_Metric_plots.pdf"))
MetricPlots(bam.test.obj)
dev.off()    
}



# 4. Apply the whitelisted CellTags generated from assessment
bam.test.obj <- SingleCellDataWhitelist(bam.test.obj, opt$whitelist_path)


# 5. Check metric plots after whitelist filtering
# Recheck the metric similar to Step 3
if (opt$visualize) {
pdf(paste0(DIR, opt$whitelist_version, opt$save_progress_name, "_Metric_plots_after_whitelist_filtering.pdf"))
MetricPlots(bam.test.obj)
dev.off()  
}


# 6. Additional filtering
# Filter out cells with more than 20 CellTags
bam.test.obj@metric.filtered.count <- as(matrix(NA, 0, 0), "dgCMatrix")
bam.test.obj <- MetricBasedFiltering(bam.test.obj, opt$high_filter, comparison = "less")

# Filter out cells with less than 2 CellTags
bam.test.obj <- MetricBasedFiltering(bam.test.obj, opt$low_filter, comparison = "greater")


# 7. Last check of metric plots
if (opt$visualize) {
pdf(paste0(DIR, opt$whitelist_version, opt$save_progress_name, "_Metric_plots_after_whitelist_filtering_last_check.pdf"))
MetricPlots(bam.test.obj)
dev.off()
}

# 8. Clone Calling
# I. Jaccard Analysis
bam.test.obj <- JaccardAnalysis(bam.test.obj, fast = TRUE)

'''
doesnt work bc has to be not fast version of jaccard, but in that version there is a bug :)
if (opt$visualize) {
pdf(paste0(DIR, opt$whitelist_version, opt$save_progress_name, "_jaccard.pdf"))
Jac = bam.test.obj@jaccard.mtx
diag(Jac) <- 1
corrplot(Jac, method="color", order="hclust", hclust.method ="ward.D2", cl.lim=c(0,1), tl.cex=0.1)
dev.off()
}
'''

# Call clones
bam.test.obj <- CloneCalling(celltag.obj = bam.test.obj, correlation.cutoff=0.7)

print_cellTag(bam.test.obj)

saveRDS(bam.test.obj, paste0(DIR, opt$save_progress_name, ".RDS"))
saveRDS(bam.test.obj, paste0(DIR, opt$whitelist_version, opt$save_progress_name, ".RDS"))

},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)
