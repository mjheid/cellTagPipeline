suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
source("scripts/helperFunctions.R")

# command line options
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
    make_option(c("--jac_cut"),  default=0.7, 
              help="Cutoff for Jaccard similarity to be considerered for clone calling", metavar="double"),
    make_option(c("--tagged"),  default=1, 
              help="When binarizing, the number of cellTags needed to consider cell as tagged", metavar="int")
)

# read arguments and put into data.frame opt
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# catch errors occuring. When catched, print error and exit programm.
tryCatch(
{

# Output directory. Assumes that data from previous steps is in here
DIR = opt$out

# Read celltag obj of previous step.
obj = readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))

# Collapse celtags based on stardust output
obj = CellTagDataPostCollapsing(celltag.obj = obj, collapsed.rslt.file = paste0(DIR, opt$whitelist_version, opt$collapsing_name))

# Binarize celltag matrix
obj = SingleCellDataBinarization(obj, opt$tagged)

# Visualize celltag metrics at this time in the process
if (opt$visualize) {
MyMetricPlots(obj, "Metric_plots", opt)  
}

# Apply the whitelisted CellTags generated from assessment
obj = SingleCellDataWhitelist(obj, opt$whitelist_path)

# Visualize celltag metrics at this time in the process
if (opt$visualize) {
MyMetricPlots(obj, "Metric_plots_after_whitelist_filtering", opt)
}

# Filter out cells with more than 20 CellTags
obj@metric.filtered.count = as(matrix(NA, 0, 0), "dgCMatrix")
obj = MetricBasedFiltering(obj, opt$high_filter, comparison = "less")

# Filter out cells with less than 2 CellTags
obj = MetricBasedFiltering(obj, opt$low_filter, comparison = "greater")


# Visualize celltag metrics at this time in the process
if (opt$visualize) {
MyMetricPlots(obj, "Metric_plots_after_all_filtering", opt)
}

# Jaccard Analysis, if fast = FALSE there is an error and execution stops?
obj = JaccardAnalysis(obj, fast = TRUE)


# Visualize Jaccard simiarity
if (opt$visualize) {
pdf(paste0(DIR, opt$whitelist_version, opt$save_progress_name, "_jaccard.pdf"))
Jac = as.matrix(obj@jaccard.mtx)
corrplot(Jac, method="color", order="hclust", hclust.method ="ward.D2", cl.lim=c(0,1), tl.cex=0.1)
dev.off()
}


# Call clones
obj = CloneCalling(celltag.obj = obj, correlation.cutoff=opt$jac_cut)

# Print data stats of cellTag obj
print_cellTag(obj)

# Save celltag obj, once for the executed library version, once for all
saveRDS(obj, paste0(DIR, opt$save_progress_name, ".RDS"))
saveRDS(obj, paste0(DIR, opt$whitelist_version, opt$save_progress_name, ".RDS"))

},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)
