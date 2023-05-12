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

DIR <- opt$out

bam.test.obj <- readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

bam.test.obj <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj, collapsed.rslt.file = paste0(DIR, opt$whitelist_version, opt$collapsing_name))

bam.test.obj <- SingleCellDataBinarization(bam.test.obj, opt$tagged)


if (opt$visualize) {
pdf(paste0(DIR, opt$save_progress_name, "_Metric_plots.pdf"))
MetricPlots(bam.test.obj)
dev.off()    
}



# 4. Apply the whitelisted CellTags generated from assessment
bam.test.obj <- SingleCellDataWhitelist(bam.test.obj, opt$whitelist_path)


# 5. Check metric plots after whitelist filtering
# Recheck the metric similar to Step 3
if (opt$visualize) {
pdf(paste0(DIR, opt$save_progress_name, "_Metric_plots_after_whitelist_filtering.pdf"))
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
pdf(paste0(DIR, opt$save_progress_name, "_Metric_plots_after_whitelist_filtering_last_check.pdf"))
MetricPlots(bam.test.obj)
dev.off()
}

# 8. Clone Calling
# I. Jaccard Analysis
bam.test.obj <- JaccardAnalysis(bam.test.obj, fast = T)

# Call clones
bam.test.obj <- CloneCalling(celltag.obj = bam.test.obj, correlation.cutoff=0.7)

saveRDS(bam.test.obj, paste0(DIR, opt$save_progress_name, ".RDS"))
saveRDS(bam.test.obj, paste0(OUTDIR, opt$whitelist_version, opt$save_progress_name, ".RDS"))
write.csv(as.matrix(bam.test.obj@whitelisted.count), file=paste0(DIR, opt$save_progress_name, ".csv"))

