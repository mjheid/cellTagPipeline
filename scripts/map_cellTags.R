suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))

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

    
"""
Plot  existing clone in dimension reduction plot.

ct: Seurat obj; Contains dimension reduction data and more
version: str; Library version of clone
method: str, either umap or tsne, data reduction method
subset: vector; contains cell barcodes of clones
subset_clust: vector; contains cluster id of clones
"""
plot_umap_with_specific_clone = function(ct, version, method, subset, subset_clust) {
    # Create a data frame with reduction values from Seurat object
    if  (method=="umap") {
        red_df = as.data.frame(ct@reductions$umap@cell.embeddings)
    } else if (method=="tsne") {
        red_df = as.data.frame(ct@reductions$tsne@cell.embeddings)
    }
    
    # Creat subset column of dataframe and set it equal to Not specified Cells
    red_df$subset = "Not specified Cells"
    
    # Set subset to clone cluster id for cell barcodes od corresponding cells
    subset = subset[subset_clust==opt$vis_clone]
    
    # Assign colors to subsets in the data frame
    red_df$subset[rownames(red_df) %in% subset] = paste0("Clone ", opt$vis_clone)
    red_df$alpha <- 1
    red_df$alpha[!rownames(red_df) %in% subset] = NA


    # Plot the UMAP values with colored subsets
    if  (method=="umap") {
        plot = ggplot(red_df, aes(UMAP_1, UMAP_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c("red", "grey")) +
          labs(color = "CellTags")
    } else if (method=="tsne") {
        plot = ggplot(red_df, aes(TSNE_1, TSNE_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c("red", "grey")) +
          labs(color = "CellTags")
    }

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_", method,  "_", version,  
                  "_", opt$vis_clone ,".pdf"), plot, width = 6, height = 6, dpi = 1)
}


"""
Plot  all existing clones in dimension reduction plot.

ct: Seurat obj; Contains dimension reduction data and more
version: str; Library version of clone
method: str, either umap or tsne, data reduction method
subset: vector; contains cell barcodes of clones
subset_clust: vector; contains cluster id of clones
"""
plot_umap_clones <- function(ct, version, method, subset, subset_clust) {
    # Create a data frame with reduction values from Seurat object
    if  (method=="umap") {
        red_df = as.data.frame(ct@reductions$umap@cell.embeddings)
    } else if (method=="tsne") {
        red_df = as.data.frame(ct@reductions$tsne@cell.embeddings)
    }
    
    # Creat subset column of dataframe and set it equal to Not specified Cells
    red_df$subset = "Not specified Cells"
    
    # Set colour of all cells to 0. Cells subset get assigned their cluster id as colour.
    red_df$col = 0
    indices_sub = which(rownames(red_df) %in% subset)
    red_df$col[indices_sub] = subset_clust
    red_df$col = as.factor(red_df$col)
    
    # Set column subset eual to clone ID. TODO: Not needed here no?
    red_df$subset[rownames(red_df) %in% subset] = paste0("Clone ", opt$vis_clone)
    
    #Set all cells who dont have a celltag to column alpha=1, otherwise NA
    red_df$alpha = 1
    red_df$alpha[!rownames(red_df) %in% subset] = NA


    # Plot the UMAP values with colored subsets
    if  (method=="umap") {
        plot = ggplot(red_df, aes(UMAP_1, UMAP_2, color = col)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1), show.legend = FALSE) +
          scale_color_discrete() +
          labs(color = "CellTags")
    } else if (method=="tsne") {
        plot = ggplot(red_df, aes(TSNE_1, TSNE_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(red_df$alpha), 1, 0.1)) +
          scale_color_discrete() +
          labs(color = "CellTags")
    }

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_", method,  "_", version,  
            "_allclones.pdf"), plot, width = 6, height = 6, dpi = 1)
}


"""
Visualize gene count of specified gene in cell containing cellTags on dimension reduction plot.

ct: Seurat obj; Contains dimension reduction data and more
version: str; Library version of clone
method: str, either umap or tsne, data reduction method
subset: vector; contains cell barcodes of clones
gene_sp: str; gene name to visualize abundance in clones
"""
plot_gene_plus_cellTag = function(ct, version, method, subset, gene_sp) {
    # Create a data frame with reduction values from Seurat object
    if  (method=="umap") {
        df = as.data.frame(ct@reductions$umap@cell.embeddings)
    } else if (method=="tsne") {
        df = as.data.frame(ct@reductions$tsne@cell.embeddings)
    }
    
    # Set barcode column equal to rownames of dataframe
    df$barcode = rownames(df)
    
    # Fetch gene count of gene_sp, filter out all cell barcodes who are not in subset
    gene = FetchData(ct, vars=gene_sp)
    gene$barcode = rownames(gene)
    filtered_gene_count = gene[gene$barcode %in% subset, ]
    
    # Partition df into dataframe that contain subset cell barcodes, and do not.
    pos_filtered_gene_count = df[rownames(df) %in% rownames(filtered_gene_count), ]
    df = df[!rownames(df) %in% rownames(filtered_gene_count), ]
    
    # Merge dataframe containing subset cell barcodes and gene count by cell barcode
    # Name gene count column of resulting dataframe gene_count
    pos_filtered_gene_count = merge(pos_filtered_gene_count, filtered_gene_count, by="barcode")
    pos_filtered_gene_count$gene_count = pos_filtered_gene_count[, 4]

    # Plot the UMAP values with colored subsets
    if  (method=="umap") {
        p = ggplot(data = df, aes(UMAP_1, UMAP_2)) +
                geom_point(color = "gray")
        p = p + geom_point(data = pos_filtered_gene_count, aes(UMAP_1, UMAP_2, color = gene_count)) +
                  scale_color_gradient(low = "lightblue", high = "darkblue")+
                  labs(color = paste0("CellTags with ", gene_sp))
        print(p)
    } else if (method=="tsne") {
        p = ggplot(data = df, aes(TSNE_1, TSNE_2)) +
                geom_point(color = "gray")
        p = p + geom_point(data = pos_filtered_gene_count, aes(TSNE_1, TSNE_2, color = gene_count)) +
                  scale_color_gradient(low = "lightblue", high = "darkblue")+
                  labs(color = paste0("CellTags with ", gene_sp))
        print(p)
    }

    # Save the plot as an image file
    ggsave(paste0(opt$out, opt$save_progress_name, "_", method,  "_", version,  
            "_", gene_sp, ".pdf"), p, width = 6, height = 6, dpi = 0.1)
}
    
    
# Read complete celltag obj of previous steps.
bam.test.obj = readRDS(paste0(opt$out, opt$save_progress_name, ".RDS"))

# If the pipeline has been run with library v1, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))) {
    bam.v1.obj = readRDS(paste0(opt$out, "v1", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v1 = bam.v1.obj@clone.composition$v1
    v1 = bam.test.obj@clone.composition$v1$cell.barcode
    v1_clones = bam.test.obj@clone.composition$v1$clone.id
}
# If the pipeline has been run with library v2, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))) {
    bam.v2.obj = readRDS(paste0(opt$out, "v2", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v2 = bam.v2.obj@clone.composition$v2
    v2 = bam.test.obj@clone.composition$v2$cell.barcode
    v2_clones = bam.test.obj@clone.composition$v2$clone.id
}
# If the pipeline has been run with library v3, load the version object and get clone cell barcodes and clone ids.
if (file.exists(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))) {
    bam.v3.obj = readRDS(paste0(opt$out, "v3", opt$save_progress_name, ".RDS"))
    bam.test.obj@clone.composition$v3 = bam.v3.obj@clone.composition$v3
    v3 = bam.test.obj@clone.composition$v3$cell.barcode
    v3_clones = bam.test.obj@clone.composition$v3$clone.id
}

# If the object already exists, and a reduction therefore has been calculated,  load the object.
if (file.exists(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))) {
    ct = readRDS(paste0(opt$out, opt$save_progress_name, "_reduction.RDS"))
}
    
# If we choose to perform filtering and prepareing data for dimensionalty reduction
if (opt$filter) {
    # Read RNA-seq data and create seurat object
    ct = Read10X(data.dir = paste0("data/samples/", opt$sample_name, "/outs/filtered_feature_bc_matrix/"))
    ct = CreateSeuratObject(counts = ct, project = "test", min.cells = opt$min_cells, min.features = opt$min_features)
    
    # Normalize data
    ct = NormalizeData(ct, normalization.method = "LogNormalize", scale.factor = opt$scale_factor)
    
    # Find varible features
    ct = FindVariableFeatures(ct, selection.method = "vst", nfeatures = opt$n_var_features)

    # Scale data
    all.genes = rownames(ct)
    ct = ScaleData(ct, features = all.genes)

    # Run PCA
    ct = RunPCA(ct)

    # Find neighbors and clusters in data
    ct = FindNeighbors(ct, dims = 1:opt$neighbours_dims)
    ct = FindClusters(ct, resolution = opt$cluster_resolution)
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
if (opt$visualize_umap) {
    # If library version is v1 or all libraries, visualize clones of v1 in UMAP
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v1.pdf"))
    g = DimPlot(ct, reduction = "umap", cells.highlight=v1, cols.highlight= "red", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v2 or all libraries, visualize clones of v2 in UMAP
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v2.pdf"))
    g = DimPlot(ct, reduction = "umap", cells.highlight=v2, cols.highlight= "blue", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v3 or all libraries, visualize clones of v3 in UMAP
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v3.pdf"))
    g = DimPlot(ct, reduction = "umap", cells.highlight=v2, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    
    # If library version is v1 or all libraries, visualize specific clone of v1 in UMAP
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v1", "umap", v1, v1_clones)
    }
    # If library version is v2 or all libraries, visualize specific clone of v2 in UMAP
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v2", "umap", v2, v2_clones)
    }
    # If library version is v3 or all libraries, visualize specific clone of v3 in UMAP
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v3", "umap", v3, v3_clones)
    }
    
    # If library version is v1 or all libraries, colour all clones of v1 in UMAP
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v1", "umap", v1, v1_clones)
    }
    # If library version is v12or all libraries, colour all clones of v2 in UMAP
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v2", "umap", v2, v2_clones)
    }
    # If library version is v3 or all libraries, colour all clones of v3 in UMAP
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v3", "umap", v3, v3_clones)
    }

    # If we have data of all libraries, visualize clones of the libraries in UMAP
    if (opt$visualize_ver=="all") {
        # Get UMAP embedings of entire data and set subset column to Without CellTag
        umap_df = as.data.frame(ct@reductions$umap@cell.embeddings)
        umap_df$subset <- "Without cellTag"

        # Define subsets and corresponding colors
        subsets = list(v1, v2, v3)
        colors = c("red", "blue", "green")

        # Assign colors to subsets in the data frame
        umap_df$alpha = NA
        for (i in seq_along(subsets)) {
            umap_df$subset[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
            umap_df$alpha[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
        }

        # Plot the UMAP values with colored subsets
        plot = ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c(colors, "black")) +
          labs(color = "CellTags")

        # Save the plot as an image file
        ggsave(paste0(opt$out, opt$save_progress_name, "_umap.pdf"), plot, width = 6, height = 6, dpi = 1)
    }
    
    #Plot specific gene abundance
    if (file.exists(opt$gene_list)) {
        # Read the file and store the gene data into a data frame
        gene_data = read.table(opt$gene_list, header = FALSE, sep = "\t")

        # Combine the gene columns into a vector
        genes = unlist(gene_data)

        # Iterate over the genes in a for loop
        for (gene in genes) {
            # Catch error if gene is not found in data and print it, do not exit
            tryCatch ({
                # Create gene Feature plot of entire dataset
                pdf(paste0(opt$out, opt$save_progress_name, "_", gene, "_umap.pdf"))
                g = FeaturePlot(ct, features = c(gene), reduction = "umap")
                print(g)
                dev.off()

                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v1", "umap", v1, gene)
                }
                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v2", "umap", v2, gene)
                }
                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v3", "umap", v3, gene)
                }
            
            },
            error = function(err) {
                # Print the error message
                print(paste("Error:", conditionMessage(err)))
            })

        }

    }

}

# Visualize UMAP
if (opt$visualize_tsne) {
    # If library version is v1 or all libraries, visualize clones of v1 in TSNE
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v1.pdf"))
    g = DimPlot(ct, reduction = "tsne", cells.highlight=v1, cols.highlight= "red", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v2 or all libraries, visualize clones of v2 in TSNE
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v2.pdf"))
    g = DimPlot(ct, reduction = "tsne", cells.highlight=v2, cols.highlight= "blue", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v2 or all libraries, visualize clones of v2 in TSNE
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_umap_v3.pdf"))
    g = DimPlot(ct, reduction = "tsne", cells.highlight=v2, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    
    # If library version is v1 or all libraries, visualize specific clone of v1 in TSNE
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v1", "tsne", v1, v1_clones)
    }
    # If library version is v2 or all libraries, visualize specific clone of v2 in TSNE
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v2", "tsne", v2, v2_clones)
    }
    # If library version is v3 or all libraries, visualize specific clone of v3 in TSNE
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v3", "tsne", v3, v3_clones)
    }
    
    # If library version is v1 or all libraries, colour all clones of v1 in TSNE
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v1", "tsne", v1, v1_clones)
    }
    # If library version is v2 or all libraries, colour all clones of v2 in TSNE
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v2", "tsne", v2, v2_clones)
    }
    # If library version is v3 or all libraries, colour all clones of v3 in TSNE
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v3", "tsne", v3, v3_clones)
    }

    # If we have data of all libraries, visualize clones of the libraries in TSNE
    if (opt$visualize_ver=="all") {
        # Get TSNE embedings of entire data and set subset column to Without CellTag
        umap_df = as.data.frame(ct@reductions$tsne@cell.embeddings)
        umap_df$subset = "Without cellTag"

        # Define subsets and corresponding colors
        subsets = list(v1, v2, v3)
        colors = c("red", "blue", "green")

        # Assign colors to subsets in the data frame
        umap_df$alpha = NA
        for (i in seq_along(subsets)) {
            umap_df$subset[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
            umap_df$alpha[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
        }

        # Plot the TSNE values with colored subsets
        plot = ggplot(umap_df, aes(TSNE_1, TSNE_2, color = subset)) +
          geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
          scale_color_manual(values = c(colors, "black")) +
          labs(color = "CellTags")

        # Save the plot as an image file
        ggsave(paste0(opt$out, opt$save_progress_name, "_tsne.pdf"), plot, width = 6, height = 6, dpi = 1)
    }
    
    #Plot specific gene abundance
    if (file.exists(opt$gene_list)) {
        # Read the file and store the gene data into a data frame
        gene_data = read.table(opt$gene_list, header = FALSE, sep = "\t")

        # Combine the gene columns into a vector
        genes = unlist(gene_data)

        # Iterate over the genes in a for loop
        for (gene_sp in genes) {
            # Catch error if gene is not found in data and print it, do not exit
            tryCatch ({
                # Create gene Feature plot of entire dataset
                pdf(paste0(opt$out, opt$save_progress_name, "_", gene_sp, "_tsne.pdf"))
                g = FeaturePlot(ct, features = c(gene_sp), reduction = "tsne")
                print(g)
                dev.off()

                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v1", "tsne", v1, gene_sp)
                }
                # If library version is v2 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v2", "tsne", v2, gene_sp)
                }
                # If library version is v3 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v3", "tsne", v3, gene_sp)
                }
            
            },
            error = function(err) {
                # Print the error message
                print(paste("Error:", conditionMessage(err)))
            })
        }

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
