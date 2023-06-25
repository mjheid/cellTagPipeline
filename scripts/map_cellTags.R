suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
source("scripts/helperFunctions.R")

# Command line options
option_list = list(
    make_option(c("--visualize_umap"), action="store_true", default=FALSE, 
              help="Visualize filtering using umap"),
    make_option(c("--visualize_tsne"), action="store_true", default=FALSE, 
              help="Visualize filtering using tsne"),
    make_option(c("--runUMAP"), action="store_true", default=FALSE, 
              help="Run UMAP"),
    make_option(c("--runTSNE"), action="store_true", default=FALSE, 
              help="Run TSNE"),
    make_option(c("--filter"), action="store_true", default=FALSE, 
              help="Filter data"),
    make_option(c("--jackstraw"), action="store_true", default=FALSE, 
              help="Use jackstraw to confirm PCs p-value, otherwise elbow plot"),
    make_option(c("--visualize_ver"), type="character", default="all", 
              help="visualize specified version, either v1, v2, v3 or all.", metavar="character"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--gene_list"), type="character", default="empty", 
              help="gene list file location and name", metavar="character"),
    make_option(c("--vis_clone"),  default=1, 
              help="Clone to visualize.", metavar="int"),
    make_option(c("--min_cells"),  default=3, 
              help="Minimum cell count to load data.", metavar="int"),
    make_option(c("--min_features"),  default=200, 
              help="Minimum features count to load data.", metavar="int"),
    make_option(c("--scale_factor"),  default=10000, 
              help="Factor to scale cells with.", metavar="int"),
    make_option(c("--npcs"),  default=100, 
              help="Number of PCs to calculate.", metavar="int"),
    make_option(c("--n_var_features"),  default=2000, 
              help="Number of varible features to work with.", metavar="int"),
    make_option(c("--neighbours_dims"),  default=10, 
              help="Number of PCAs for FindNeighbours to look at.", metavar="int"),
    make_option(c("--cluster_resolution"),  default=0.5, 
              help="Resolution of clusters to calculate. Between 0 and 1.", metavar="float"),
    make_option(c("--save_progress_name"), type="character", default="insert", 
              help="Name appended to outputs to be saved", metavar="character"),
    make_option(c("--sample_name"), type="character", default="SI-TT-H4", 
                help="Name of experiment sample", metavar="character")
)

# read arguments and put into data.frame opt
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# catch errors occuring. When catched, print error and exit programm.
tryCatch(
{

# Read complete celltag obj of previous steps.
obj = readRDS(paste0(opt$out, opt$save_progress_name, ".RDS"))

# If the pipeline has been run with library v1, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))) {
    v1.obj = readRDS(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))
    obj@clone.composition$v1 = v1.obj@clone.composition$v1
    v1 = obj@clone.composition$v1$cell.barcode
    v1_clones = obj@clone.composition$v1$clone.id
    v1_c = rownames(v1.obj@jaccard.mtx)
}
# If the pipeline has been run with library v2, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))) {
    v2.obj = readRDS(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))
    obj@clone.composition$v2 = v2.obj@clone.composition$v2
    v2 = obj@clone.composition$v2$cell.barcode
    v2_clones = obj@clone.composition$v2$clone.id
    v2_c = rownames(v2.obj@jaccard.mtx)
}
# If the pipeline has been run with library v3, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))) {
    v3.obj = readRDS(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))
    obj@clone.composition$v3 = v3.obj@clone.composition$v3
    v3 = obj@clone.composition$v3$cell.barcode
    v3_clones = obj@clone.composition$v3$clone.id
    v3_c = rownames(v3.obj@jaccard.mtx)
}

# If the object already exists, and a reduction therefore has been calculated,  load the object.
if (file.exists(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))) {
    ct = readRDS(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}
    
# If we choose to perform filtering and prepareing data for dimensionalty reduction
if (opt$filter) {
    # Read RNA-seq data and create seurat object
    ct = Read10X(data.dir = paste0("data/samples/", opt$sample_name, "/outs/filtered_feature_bc_matrix/"))
    ct = CreateSeuratObject(counts = ct, project = "PVNct", min.cells = opt$min_cells, min.features = opt$min_features)
    
    # Find mitochrondrial genes
    ct[["percent.mt"]] = PercentageFeatureSet(ct, pattern = "^MT-")
    
    # QC plots
    pdf(paste0(opt$out, opt$save_progress_name, "_QC.pdf"))
    g = VlnPlot(ct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(g)
    dev.off()
    pdf(paste0(opt$out, opt$save_progress_name, "_QC_FeatureScatter.pdf"))
    plot1 = FeatureScatter(ct, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 = FeatureScatter(ct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    g = plot1 + plot2
    print(g)
    dev.off()
    
    
    # Normalize data
    ct = NormalizeData(ct, normalization.method="LogNormalize", scale.factor=opt$scale_factor)
    
    # Find varible features
    ct = FindVariableFeatures(ct, nfeatures=opt$n_var_features)
    
    # Plot itttt
    top10 = head(VariableFeatures(ct), 10)
    # plot variable features with and without labels
    pdf(paste0(opt$out, opt$save_progress_name, "_variableFeatures.pdf"))
    plot1 = VariableFeaturePlot(ct)
    plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
    g = plot1 + plot2
    print(g)
    dev.off()
    
    # Score genes for cell cycleing
    s.genes = cc.genes$s.genes
    g2m.genes = cc.genes$g2m.genes
    ct = CellCycleScoring(ct, s.features=s.genes, g2m.features=g2m.genes)
    
    # Recommended as only regressing genes can negatively impact downstream analysis, particularly in differentiating processes 
    ct$CC.Difference = ct$S.Score - ct$G2M.Score
    
    # Scale data
    ct = ScaleData(ct, features=rownames(ct), vars.to.regress=c("CC.Difference", "percent.mt"))

    # Run PCA on varible features
    ct = RunPCA(ct, features=VariableFeatures(object=ct), npcs=opt$npcs)
    
    # Some Plots
    pdf(paste0(opt$out, opt$save_progress_name, "_pca_dimheatmap.pdf"))
    g = DimHeatmap(ct, dims = 1:15, cells = 500, balanced = TRUE)
    print(g)
    dev.off
    
    # Test significance of PCs
    if (opt$jackstraw) {
        ct <- JackStraw(ct, num.replicate = 100)
        ct <- ScoreJackStraw(ct, dims = 1:90)
        pdf(paste0(opt$out, opt$save_progress_name, "_jackstraw.pdf"))
        g = JackStrawPlot(ct, dims = 1:90)
        print(g)
        dev.off()
    } else {
        pdf(paste0(opt$out, opt$save_progress_name, "_elbow.pdf"))
        g = ElbowPlot(ct)
        print(g)
        dev.off()
    }

    # Find neighbors and clusters in data
    ct = FindNeighbors(ct, dims = 1:opt$neighbours_dims)
    ct = FindClusters(ct, resolution = opt$cluster_resolution)
    
    saveRDS(ct, file = paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}

#  Run UMAP and save seurat object
if (opt$runUMAP) {
    ct = RunUMAP(ct, dims = 1:opt$neighbours_dims)
    saveRDS(ct, file = paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}

# Run TSNE and save seurat object
if (opt$runUMAP) {
    ct = RunTSNE(ct)
    saveRDS(ct, file = paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}


# Visualize UMAP
if (opt$visualize_umap | opt$visualize_tsne) {
    plotMap(opt, ct, v1=v1, v2=v2, v3=v3, v1_clones=v1_clones, v2_clones=v2_clones, 
            v3_clones=v3_clones, v1_c=v1_c, v2_c=v2_c, v3_c=v3_c)
}

    
},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)
