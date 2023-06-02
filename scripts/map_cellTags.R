suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))

option_list = list(
    make_option(c("--visualize_umap"), action="store_true", default=FALSE, 
              help="Visualize filtering using umap"),
    make_option(c("--visualize_tsne"), action="store_true", default=FALSE, 
              help="Visualize filtering using tsne"),
    make_option(c("--runUMAP"), action="store_true", default=FALSE, 
              help="Run UMAP"),
    make_option(c("--runTSNE"), action="store_true", default=FALSE, 
              help="Run TSNE"),
    make_option(c("--visualize_ver"), type="character", default="all", 
              help="visualize specified version, either v1, v2, v3 or all.", metavar="character"),
    make_option(c("--out"), type="character", default="data/out/", 
              help="output file dir", metavar="character"),
    make_option(c("--vis_clone"),  default=1, 
              help="Clone to visualize.", metavar="int"),
    make_option(c("--min_cells"),  default=3, 
              help="Minimum cell count to load data.", metavar="int"),
    make_option(c("--min_features"),  default=200, 
              help="Minimum features count to load data.", metavar="int"),
    make_option(c("--scale_factor"),  default=10000, 
              help="Factor to scale cells with.", metavar="int"),
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

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


tryCatch(
{


plot_umap_with_specific_clone <- function(pbmc, version, method, subset, subset_clust) {
    # Create a data frame with reduction values from Seurat object
    if  (method=="umap") {
        red_df <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
    } else if (method=="tsne") {
        red_df <- as.data.frame(pbmc@reductions$tsne@cell.embeddings)
    }
    
    red_df$subset <- "Not specified Cells"
    
    subset = subset[subset_clust==opt$vis_clone]
    
    # Assign colors to subsets in the data frame
    red_df$subset[rownames(red_df) %in% subset] <- paste0("Clone ", opt$vis_clone)
    red_df$alpha <- 1
    red_df$alpha[!rownames(red_df) %in% subset] <- NA


    # Plot the UMAP values with colored subsets
    if  (method=="umap") {
        plot <- ggplot(red_df, aes(UMAP_1, UMAP_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c("red", "grey")) +
          labs(color = "CellTags")
    } else if (method=="tsne") {
        plot <- ggplot(red_df, aes(TSNE_1, TSNE_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c("red", "grey")) +
          labs(color = "CellTags")
    }

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_", method,  "_", version,  
                  "_", opt$vis_clone ,".pdf"), plot, width = 6, height = 6, dpi = 1)
}


plot_umap_clones <- function(pbmc, version, method, subset, subset_clust) {
    # Create a data frame with reduction values from Seurat object
    if  (method=="umap") {
        red_df <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
    } else if (method=="tsne") {
        red_df <- as.data.frame(pbmc@reductions$tsne@cell.embeddings)
    }
    
    red_df$subset <- "Not specified Cells"
    red_df$col <- 0
    indices_sub = which(rownames(red_df) %in% subset)
    red_df$col[indices_sub] = subset_clust
    red_df$col = as.factor(red_df$col)
    
    # Assign colors to subsets in the data frame
    red_df$subset[rownames(red_df) %in% subset] <- paste0("Clone ", opt$vis_clone)
    red_df$alpha <- 1
    red_df$alpha[!rownames(red_df) %in% subset] <- NA


    # Plot the UMAP values with colored subsets
    if  (method=="umap") {
        plot <- ggplot(red_df, aes(UMAP_1, UMAP_2, color = col)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1), show.legend = FALSE) +
          scale_color_discrete() +
          labs(color = "CellTags")
    } else if (method=="tsne") {
        plot <- ggplot(red_df, aes(TSNE_1, TSNE_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_discrete() +
          labs(color = "CellTags")
    }

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_", method,  "_", version,  
            "_allclones.pdf"), plot, width = 6, height = 6, dpi = 1)
}



bam.test.obj <- readRDS(paste0(opt$out, opt$save_progress_name, ".RDS"))

if (file.exists(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))) {
    bam.v1.obj <- readRDS(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v1 = bam.v1.obj@clone.composition$v1
    v1 = bam.test.obj@clone.composition$v1$cell.barcode
    v1_clones = bam.test.obj@clone.composition$v1$clone.id
}
if (file.exists(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))) {
    bam.v2.obj <- readRDS(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v2 = bam.v2.obj@clone.composition$v2
    v2 = bam.test.obj@clone.composition$v2$cell.barcode
    v2_clones = bam.test.obj@clone.composition$v2$clone.id
}
if (file.exists(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))) {
    bam.v3.obj <- readRDS(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v3 = bam.v3.obj@clone.composition$v3
    v3 = bam.test.obj@clone.composition$v3$cell.barcode
    v3_clones = bam.test.obj@clone.composition$v3$clone.id
}

if (file.exists(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))) {
    pbmc = readRDS(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
} else {
    pbmc.data <- Read10X(data.dir = paste0("data/samples/", opt$sample_name, "/outs/filtered_feature_bc_matrix/"))
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "test", min.cells = opt$min_cells, min.features = opt$min_features)

    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = opt$scale_factor)

    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = opt$n_var_features)

    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)

    pbmc <- RunPCA(pbmc)

    pbmc <- FindNeighbors(pbmc, dims = 1:opt$neighbours_dims)
    pbmc <- FindClusters(pbmc, resolution = opt$cluster_resolution)
}

if (opt$runUMAP) {
    pbmc <- RunUMAP(pbmc, dims = 1:opt$neighbours_dims)
    saveRDS(pbmc, file = paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}
if (opt$runUMAP) {
    pbmc <- RunTSNE(pbmc)
    saveRDS(pbmc, file = paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}



if (opt$visualize_umap) {
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v1.pdf"))
    g = DimPlot(pbmc, reduction = "umap", cells.highlight=v1, cols.highlight= "red", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v2.pdf"))
    g = DimPlot(pbmc, reduction = "umap", cells.highlight=v2, cols.highlight= "blue", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v3.pdf"))
    g = DimPlot(pbmc, reduction = "umap", cells.highlight=v2, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v1", "umap", v1, v1_clones)
    }
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v2", "umap", v2, v2_clones)
    }
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v3", "umap", v3, v3_clones)
    }
    
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v1", "umap", v1, v1_clones)
    }
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v2", "umap", v2, v2_clones)
    }
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v3", "umap", v3, v3_clones)
    }

    # Create a data frame with UMAP values from Seurat object
    if (opt$visualize_ver=="all") {
    umap_df <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
    umap_df$subset <- "Without cellTag"

    # Define subsets and corresponding colors
    subsets <- list(v1, v2, v3)
    colors <- c("red", "blue", "green")

    # Assign colors to subsets in the data frame
    umap_df$alpha <- NA
    for (i in seq_along(subsets)) {
        umap_df$subset[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
        umap_df$alpha[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
    }

    # Plot the UMAP values with colored subsets
    plot <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subset)) +
      geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
      scale_color_manual(values = c(colors, "black")) +
      labs(color = "CellTags")

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_umap.pdf"), plot, width = 6, height = 6, dpi = 1)
    }

}

if (opt$visualize_tsne) {
    f (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v1.pdf"))
    g = DimPlot(pbmc, reduction = "tsne", cells.highlight=v1, cols.highlight= "red", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v2.pdf"))
    g = DimPlot(pbmc, reduction = "tsne", cells.highlight=v2, cols.highlight= "blue", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v3.pdf"))
    g = DimPlot(pbmc, reduction = "tsne", cells.highlight=v2, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v1", "tsne", v1, v1_clones)
    }
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v2", "tsne", v2, v2_clones)
    }
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(pbmc, "v3", "tsne", v3, v3_clones)
    }
    
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v1", "tsne", v1, v1_clones)
    }
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v2", "tsne", v2, v2_clones)
    }
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_clones(pbmc, "v3", "tsne", v3, v3_clones)
    }

    # Create a data frame with UMAP values from Seurat object
    if (opt$visualize_ver=="all") {
    umap_df <- as.data.frame(pbmc@reductions$tsne@cell.embeddings)
    umap_df$subset <- "Without cellTag"

    # Define subsets and corresponding colors
    subsets <- list(v1, v2, v3)
    colors <- c("red", "blue", "green")

    # Assign colors to subsets in the data frame
    umap_df$alpha <- NA
    for (i in seq_along(subsets)) {
        umap_df$subset[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
        umap_df$alpha[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
    }

    # Plot the UMAP values with colored subsets
    plot <- ggplot(umap_df, aes(TSNE_1, TSNE_2, color = subset)) +
      geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
      scale_color_manual(values = c(colors, "black")) +
      labs(color = "CellTags")

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_umap.pdf"), plot, width = 6, height = 6, dpi = 1)
    }
}

},
error = function(err) {
    # Print the error message
    print(paste("Error:", conditionMessage(err)))
    # Stop the execution of the program
    stop("Exiting the program.")
}
)
