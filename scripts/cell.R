suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
 
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
    make_option(c("--bam_data"), type="character", default="data/bam/", 
                help="Path from current dir to output dir of filtered .bam files", metavar="character"),
    make_option(c("--bamfilter"), action="store_true", default=FALSE, 
                help="Path from current dir to output dir of filtered .bam files")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

DIR <- opt$bam_data
OUTDIR <- opt$out

bamfilter=""

if (opt$bamfilter){
    bamfilter = "_merged"
    print("HEEEEEEEEEEEELLLOOOOOOO????")
}

# Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = paste0(DIR, "possorted_genome_bam.filtered", bamfilter, ".bam"))


bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = opt$whitelist_version)


bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj,
	barcodes.file = paste0("data/", opt$sample_name, "/outs/filtered_feature_bc_matrix/barcodes.tsv"))


bam.test.obj <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = paste0(OUTDIR, opt$save_progress_name, "_", opt$collapsing_name))
saveRDS(bam.test.obj, paste0(OUTDIR, opt$save_progress_name, ".RDS"))