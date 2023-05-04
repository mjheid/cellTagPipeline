library("CellTagR")

# Read in the data file that come with the package
fpath <- system.file("extdata", "V2-1_R1.zip", package = "CellTagR")
extract.dir <- "."
# Extract the dataset
unzip(fpath, overwrite = FALSE, exdir = ".")
full.fpath <- paste0(extract.dir, "/", "V2-1_S2_L001_R1_001.fastq")
# Set up the CellTag Object
test.obj <- CellTagObject(object.name = "v2.whitelist.test", fastq.bam.directory = full.fpath)
# Extract the CellTags
test.obj <- CellTagExtraction(celltag.obj = test.obj, celltag.version = "v2")

# Count and Sort the CellTags in descending order of occurrence
test.obj <- AddCellTagFreqSort(test.obj)
# Check the stats
test.obj@celltag.freq.stats

# Generate the whitelist
test.obj <- CellTagWhitelistFiltering(celltag.obj = test.obj, percentile = 0.9, output.dir = NULL)
