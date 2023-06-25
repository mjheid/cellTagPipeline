suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
suppressMessages(library("ggplot2"))


# Plot  existing clone in dimension reduction plot.

# ct: Seurat obj; Contains dimension reduction data and more
# version: str; Library version of clone
# method: str, either umap or tsne, data reduction method
# subset: vector; contains cell barcodes of clones
# subset_clust: vector; contains cluster id of clones
# opt: data.frame; contains arguments of script call
plot_umap_with_specific_clone = function(ct, version, method, subset, subset_clust, opt) {
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



# Plot  all existing clones in dimension reduction plot.

# ct: Seurat obj; Contains dimension reduction data and more
# version: str; Library version of clone
# method: str, either umap or tsne, data reduction method
# subset: vector; contains cell barcodes of clones
# subset_clust: vector; contains cluster id of clones
# opt: data.frame; contains arguments of script call
plot_umap_clones <- function(ct, version, method, subset, subset_clust, opt) {
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



# Visualize gene count of specified gene in cell containing cellTags on dimension reduction plot.

# ct: Seurat obj; Contains dimension reduction data and more
# version: str; Library version of clone
# method: str, either umap or tsne, data reduction method
# subset: vector; contains cell barcodes of clones
# gene_sp: str; gene name to visualize abundance in clones
# opt: data.frame; contains arguments of script call
# name: str; contains information about subset, default clones
plot_gene_plus_cellTag = function(ct, version, method, subset, gene_sp, opt, name="clones") {
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
            "_", gene_sp, "_", name, ".pdf"), p, width = 6, height = 6, dpi = 0.1)
}
    


#Visualize gene count of specified gene in specific clone on dimension reduction plot.

#ct: Seurat obj; Contains dimension reduction data and more
#version: str; Library version of clone
#method: str, either umap or tsne, data reduction method
#subset: vector; contains cell barcodes of clones
#subset_clust: vector; contains cluster id of clones, should be same length as subset
#gene_sp: str; gene name to visualize abundance in clones
#opt: data.frame; contains arguments of script call
plot_gene_plus_clone = function(ct, version, method, subset, subset_clust, gene_sp, opt) {
    # Create a data frame with reduction values from Seurat object
    tryCatch({
        if  (method=="umap") {
            df = as.data.frame(ct@reductions$umap@cell.embeddings)
        } else if (method=="tsne") {
            df = as.data.frame(ct@reductions$tsne@cell.embeddings)
        }

        # Set subset to only specific clone
        subset = subset[subset_clust==opt$vis_clone]

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
                "_", gene_sp, "_", opt$vis_clone, ".pdf"), p, width = 6, height = 6, dpi = 0.1)
    },
    error = function(err) {
        # Print the error message
        print(paste("Error:", conditionMessage(err)))
    })
    
}



# Plots all the stuff

# opt: data.frame; contains args
# v1: vector; contains clones cell barcodes
# v2: vector; contains clones cell barcodes
# v3: vector; contains clones cell barcodes
# v1_clones: vector; contains clones clone id
# v2_clones: vector; contains clones clone id
# v3_clones: vector; contains clones clone id
# v1_c: vector; contains filtered cell barcodes
# v1_c: vector; contains filtered cell barcodes
# v1_c: vector; contains filtered cell barcodes
plotMap = function(opt, ct, v1=NULL, v2=NULL, v3=NULL, v1_clones=NULL, v2_clones=NULL, v3_clones=NULL,
                  v1_c=NULL, v2_c=NULL, v3_c=NULL) {
    if (opt$visualize_umap) {
        method = "umap"
    } else {
        method = "tsne"
    }
    
    # If library version is v1 or all libraries, visualize clones of v1 in UMAP
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v1_clones.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v1, cols.highlight= "red", cols = "gray")
    print(g)
    dev.off()
    
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v1_ct.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v1_c, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v2 or all libraries, visualize clones of v2 in dim reduction
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v2_clones.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v2, cols.highlight= "blue", cols = "gray")
    print(g)
    dev.off()
    
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v2_ct.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v2_c, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    # If library version is v3 or all libraries, visualize clones of v3 in dim reduction
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v3_clones.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v3, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    
    pdf(paste0(opt$out, opt$save_progress_name, "_", method, "_v3_ct.pdf"))
    g = DimPlot(ct, reduction = method, cells.highlight=v3_c, cols.highlight= "green", cols = "gray")
    print(g)
    dev.off()
    }
    
    # If library version is v1 or all libraries, visualize specific clone of v1 in dim reduction
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v1", method, v1, v1_clones, opt)
    }
    # If library version is v2 or all libraries, visualize specific clone of v2 in dim reduction
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v2", method, v2, v2_clones, opt)
    }
    # If library version is v3 or all libraries, visualize specific clone of v3 in dim reduction
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_with_specific_clone(ct, "v3", method, v3, v3_clones, opt)
    }
    
    # If library version is v1 or all libraries, colour all clones of v1 in dim reduction
    if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v1", method, v1, v1_clones, opt)
    }
    # If library version is v12or all libraries, colour all clones of v2 in dim reduction
    if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v2", method, v2, v2_clones, opt)
    }
    # If library version is v3 or all libraries, colour all clones of v3 in dim reduction
    if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
    plot_umap_clones(ct, "v3", method, v3, v3_clones, opt)
    }

    # If we have data of all libraries, visualize clones of the libraries in dim reduction
    if (opt$visualize_ver=="all") {
        # Get embedings of entire data and set subset column to Without CellTag
        if (opt$visualize_umap) {
            umap_df = as.data.frame(ct@reductions$umap@cell.embeddings)
        } else {
            umap_df = as.data.frame(ct@reductions$tsne@cell.embeddings)
        }
        
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

        # Plot the dim reduction values with colored subsets
        if (opt$visualize_umap) {
            plot = ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subset)) +
              geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
              scale_color_manual(values = c(colors, "black")) +
              labs(color = "CellTags")
        } else {
            plot = ggplot(umap_df, aes(TSNE_1, TSNE_2, color = subset)) +
              geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
              scale_color_manual(values = c(colors, "black")) +
              labs(color = "CellTags")
        }

        # Save the plot as an image file
        ggsave(paste0(opt$out, opt$save_progress_name, "_", method, "_clones.pdf"), plot, width = 6, height = 6, dpi = 1)
        
        # Get embedings of entire data and set subset column to Without CellTag
        if (opt$visualize_umap) {
            umap_df = as.data.frame(ct@reductions$umap@cell.embeddings)
        } else {
            umap_df = as.data.frame(ct@reductions$tsne@cell.embeddings)
        }
        
        umap_df$subset <- "Without cellTag"

        # Define subsets and corresponding colors
        subsets = list(v1_c, v2_c, v3_c)
        colors = c("red", "blue", "green")

        # Assign colors to subsets in the data frame
        umap_df$alpha = NA
        for (i in seq_along(subsets)) {
            umap_df$subset[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
            umap_df$alpha[rownames(umap_df) %in% subsets[[i]]] <- paste0("v", i)
        }

        # Plot the dim reduction values with colored subsets
        if (opt$visualize_umap) {
            plot = ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subset)) +
              geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
              scale_color_manual(values = c(colors, "black")) +
              labs(color = "CellTags")
        } else {
            plot = ggplot(umap_df, aes(TSNE_1, TSNE_2, color = subset)) +
              geom_point(alpha = ifelse(!is.na(umap_df$alpha), 1, 0.1)) +
              scale_color_manual(values = c(colors, "black")) +
              labs(color = "CellTags")
        }

        # Save the plot as an image file
        ggsave(paste0(opt$out, opt$save_progress_name, "_", method, "_celltags.pdf"), plot, width = 6, height = 6, dpi = 1)
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
                pdf(paste0(opt$out, opt$save_progress_name, "_", gene, "_", method, ".pdf"))
                g = FeaturePlot(ct, features = c(gene), reduction = method)
                print(g)
                dev.off()

                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v1" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v1", method, v1, gene, opt)
                plot_gene_plus_cellTag(ct, "v1", method, v1_c, gene, opt, name="cellTags")
                plot_gene_plus_clone(ct, "v1", method, v1, v1_clones, gene, opt)
                }
                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v2" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v2", method, v2, gene, opt)
                plot_gene_plus_cellTag(ct, "v2", method, v2_c, gene, opt, name="cellTags")
                plot_gene_plus_clone(ct, "v2", method, v2, v2_clones, gene, opt)
                }
                # If library version is v1 or all libraries, create gene feature plot of cells with celltag
                if (opt$visualize_ver=="v3" || opt$visualize_ver=="all") {
                plot_gene_plus_cellTag(ct, "v3", method, v3, gene, opt)
                plot_gene_plus_cellTag(ct, "v3", method, v3_c, gene, opt, name="cellTags")
                plot_gene_plus_clone(ct, "v3", method, v3, v3_clones, gene, opt)
                }
            
            },
            error = function(err) {
                # Print the error message
                print(paste("Error:", conditionMessage(err)))
            })

        }

    }

}


#Print celltag.R object with all relevent stats of each step in the pipeline.

#object: cellTagR object;
print_cellTag = function(object) {
            cat("Object name: ", object@obj.name, "\n")
            cat("Library version: ", object@curr.version, "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "raw.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols = !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub <- curr.mtx.sub[, valid_cols]
            cat("Unique Raw CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Raw CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "raw.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub = curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Raw Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "collapsed.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols = !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub = curr.mtx.sub[, valid_cols]
            cat("Unique Collapsed CellTag Counts = ", ncol(full.mtx.sub), "\n")
            cat("Collapsed CellTag Counts = ", sum(full.mtx.sub, na.rm = TRUE), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "whitelisted.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols = !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub = curr.mtx.sub[, valid_cols]
            cat("Unique Whitelisted CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Whitelisted CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "whitelisted.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub = curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Whitelisted Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "metric.filtered.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            valid_cols = !(Matrix::colSums(curr.mtx.sub != 0, na.rm = TRUE) == 0 & Matrix::colSums(!is.na(curr.mtx.sub)) == 0)
            full.mtx.sub = curr.mtx.sub[, valid_cols]
            cat("Unique Filtered CellTag Counts = ", (ncol(full.mtx.sub)), "\n")
            cat("Filtered CellTag Counts = ", (sum(full.mtx.sub, na.rm = TRUE)), "\n")
            
            # Get non emtpy matrix
            curr.mtx = slot(object, "metric.filtered.count")
            curr.version = object@curr.version
            curr.mtx.sub = curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
            colnames(curr.mtx.sub) = gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
            full.mtx.sub = curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
            cat("Filtered Number of Cells with CellTag = ", nrow(full.mtx.sub), "\n")
            
            curr.version = object@curr.version
            if (curr.version == "v1") {
                cat("Number of identified clones = ", length(object@clone.composition$v1$cell.barcode), "\n")
                cat("Number of unique clones = ", length(unique(object@clone.composition$v1$clone.id)), "\n")
                cat("Number of clones per unique clone = ")
                print(table(object@clone.composition$v1$clone.id))
            } else if (curr.version == "v2") {
                cat("Number of identified clones = ", length(object@clone.composition$v2$cell.barcode), "\n")
                cat("Number of unique clones = ", length(unique(object@clone.composition$v2$clone.id)), "\n")
                cat("Number of clones per unique clone = ")
                print(table(object@clone.composition$v2$clone.id))
            } else {
                cat("Number of identified clones = ", length(object@clone.composition$v3$cell.barcode), "\n")
                cat("Number of unique clones = ", length(unique(object@clone.composition$v3$clone.id)), "\n")
                cat("Number of clones per unique clone = ")
                print(table(object@clone.composition$v3$clone.id))
            }
            
}


#Get working matrix of slot from cellTagR object, filters out empty rows.

#celltag.obj:  cellTagR object;
#slot.to.select: str; name of slot to select
GetCellTagCurrentVersionWorkingMatrix <- function(celltag.obj, slot.to.select) {
  curr.mtx <- slot(celltag.obj, slot.to.select)
  if (nrow(curr.mtx) <= 0) {
    return(curr.mtx)
  } else {
    curr.version <- celltag.obj@curr.version
    curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
full.mtx.sub <- curr.mtx.sub[!(Matrix::rowSums(is.na(curr.mtx.sub)) == ncol(curr.mtx.sub) | Matrix::rowSums(curr.mtx.sub == 0) == ncol(curr.mtx.sub)),]
    
    return(full.mtx.sub)
  }
}


# Plot Metrics of cellTag obj

# celltag.obj: cellTagR object;
# name: str; name to concatonate to resulting output files.
# opt: data.frame; contains arguments of script call
MyMetricPlots <- function(celltag.obj, name, opt) {
  
      obj.metric.filtered.count <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "metric.filtered.count")
      obj.whitelisted.count <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "whitelisted.count")

      if (ncol(obj.metric.filtered.count) <= 0) {
        if (ncol(obj.whitelisted.count) <= 0) {
          celltag.data <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "binary.mtx")
        } else {
          celltag.data <- obj.whitelisted.count
        }
      } else {
        celltag.data <- obj.metric.filtered.count
      }

      CellTags.per.cell.whitelisted.pf <- Matrix::rowSums(celltag.data)
      CellTags.per.cell.avg <- mean(CellTags.per.cell.whitelisted.pf)
      CellTags.frequency.whitelisted.pf <- Matrix::colSums(celltag.data)
      CellTags.freq.avg <- mean(CellTags.frequency.whitelisted.pf)

      pdf(paste0(opt$out, opt$whitelist_version, opt$save_progress_name, "_", name, "_CCIC.pdf"))
      g = plot(CellTags.per.cell.whitelisted.pf, main = "CellTag Counts of Individual Cells", xlab = "Cell Index", ylab = "CellTag Counts")
      print(g)
      dev.off()
      pdf(paste0(opt$out, opt$whitelist_version, opt$save_progress_name, "_", name, "_COFAAC.pdf"))
      g = plot(CellTags.frequency.whitelisted.pf, main = "CellTag Occurrence Frequency Across All Cells", xlab = "Cell Index", ylab = "CellTag Frequency")
      print(g)
      dev.off()

      hist_data <- hist(CellTags.per.cell.whitelisted.pf, plot = FALSE)
      counts = hist_data$counts
      midpoints = seq(1, length(counts))
      midpoints = midpoints[!counts==0]
      counts = counts[!counts==0]
      hist_df <- data.frame(counts = counts, midpoints=midpoints)
      g <- ggplot(hist_df, aes(x = factor(midpoints), y = counts)) +
              geom_bar(stat = "identity") +
              labs(x = "CellTag Counts", y = "Count") +
              ggtitle("Histogram of CellTag Counts of Individual Cells")
      ggsave(paste0(opt$out, opt$whitelist_version, opt$save_progress_name, "_", name, "_HCCIC.pdf"),
             g, width = 6, height = 6, dpi = 1)
      
      
      hist_data <- hist(CellTags.frequency.whitelisted.pf, plot = FALSE)
      counts <- hist_data$counts
      midpoints <- hist_data$mids
      midpoints <- midpoints[counts != 0]
      counts <- counts[counts != 0]
      hist_df <- data.frame(counts = counts, midpoints = midpoints)

      # Create the histogram plot using ggplot
      g <- ggplot(hist_df, aes(x = factor(midpoints), y = counts)) +
        geom_bar(stat = "identity") +
        labs(x = "CellTag Occurrence Frequency", y = "Count") +
        ggtitle("Histogram of CellTag Occurrence Frequency Across All Cells")

      # Print the plot
      ggsave(paste0(opt$out, opt$whitelist_version, opt$save_progress_name, "_", name, "_HCOFAAC.pdf"),
             g, width = 6, height = 6, dpi = 1)
      
      cat("Average: ", CellTags.per.cell.avg, "\n")
      cat("Frequency: ", CellTags.freq.avg, "\n")
}
