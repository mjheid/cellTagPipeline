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
    make_option(c("--bam_data"), type="character", default="data/bam/possorted_genome_bam.filtered.bam", 
                help="Path from current dir to .bam file", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

OUTDIR <- opt$out

# Set up the CellTag Object
if (file.exists(paste0(opt$save_prograss_name, ".RDS"))) {
    bam.test.obj <- readRDS(paste0(DIR, opt$save_progress_name, ".RDS"))
} else {
    bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = opt$bam_data)
}

bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = opt$whitelist_version)

bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj,
	barcodes.file = paste0("data/samples/", opt$sample_name, "/outs/filtered_feature_bc_matrix/barcodes.tsv"))

bam.test.obj <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = paste0(OUTDIR, opt$whitelist_version, opt$save_progress_name, "_", opt$collapsing_name))

saveRDS(bam.test.obj, paste0(OUTDIR, opt$whitelist_version, opt$save_progress_name, ".RDS"))
