suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))

# command line options
option_list = list(
    make_option(c("--out"), type="character", default="data/out/", 
                help="output file dir", metavar="character"),
    make_option(c("--collapsing_name"), type="character", default="collapsing.txt", 
                help="Name of collapsing file", metavar="character"),
    make_option(c("--save_progress_name"), type="character", default="test", 
                help="Name appended to outputs to be saved", metavar="character"),
    make_option(c("--sample_name"), type="character", default="SI-TT-H4", 
                help="Name of experiment sample", metavar="character"),
    make_option(c("--whitelist_version"), type="character", default="v1", 
                help="Version of whitelist", metavar="character"),
    make_option(c("--bam_data"), type="character", default="data/bam/possorted_genome_bam.filtered.bam", 
                help="Path from current dir to .bam file", metavar="character")
)

# read arguments and put into data.frame opt
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# catch errors occuring. When catched, print error and exit programm.
tryCatch(
{

# Set up the CellTag Object
if (file.exists(paste0(opt$out, opt$save_prograss_name, ".RDS"))) {
    # If the object already exists, load it.
    obj = readRDS(paste0(opt$out, opt$save_progress_name, ".RDS"))
} else {
    # If the object does not already exist, create it.
    obj = CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = opt$bam_data)
}

# Extract all celltags. Does simple pattern matching of beginning, tag, and end.
obj = CellTagExtraction(obj, celltag.version = opt$whitelist_version)

# Create celltag matrix cound with the cell barcodes
obj = CellTagMatrixCount(celltag.obj = obj,
	barcodes.file = paste0("data/samples/", opt$sample_name, "/outs/filtered_feature_bc_matrix/barcodes.tsv"))

# Prepare celltag data for collapsing, creates file of celltags, similar ones get collapsed.
obj = CellTagDataForCollapsing(celltag.obj = obj, output.file = paste0(opt$out, opt$whitelist_version, opt$save_progress_name, "_", opt$collapsing_name))

# Save working object, once for total object over multiple versions, once for single version.
saveRDS(obj, paste0(opt$out, opt$whitelist_version, opt$save_progress_name, ".RDS"))
saveRDS(obj, paste0(opt$out, opt$save_progress_name, ".RDS"))

},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)