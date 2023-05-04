library("CellTagR")

# Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "/home/kaies/Downloads/wgEncodeUwRepliSeqBjG1bAlnRep1.bam")

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file = "./barcodes.tsv")
# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)

# Read the RDS file and get the object
dt.mtx.path <- system.file("extdata", "Demo_V1.Rds", package = "CellTagR")
bam.test.obj <- readRDS(dt.mtx.path)

# Generating the collapsing file
bam.test.obj <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "~/Desktop/collapsing.txt")
