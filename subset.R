> suppressMessages(library("CellTagR"))
suppressMessages(library("Seurat"))
There were 20 warnings (use warnings() to see them)
Warning message:
package ‘Seurat’ was built under R version 4.2.2 
> pbmc = readRDS("data/out/anika/anika_subset_reduction.RDS")
> pbmc
An object of class Seurat 
24354 features across 230 samples within 1 assay 
Active assay: RNA (24354 features, 2000 variable features)
 3 dimensional reductions calculated: pca, umap, tsne
> pdf("subset.pdf")
> DimPlot(pbmc, reduction = "umap")
> dev.off()
null device 
          1 
> pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
Calculating cluster 0
For a more efficient implementation of the Wilcoxon Rank Sum Test,
(default method for FindMarkers) please install the limma package
--------------------------------------------
install.packages('BiocManager')
BiocManager::install('limma')
--------------------------------------------
After installation of limma, Seurat will automatically use the more 
efficient implementation (no further action necessary).
This message will be shown once per session
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
Calculating cluster 1
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
Calculating cluster 2
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
> pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
# A tibble: 6 × 7
# Groups:   cluster [3]
     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
     <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
1 5.50e-22       1.68 0.992 0.737  1.34e-17 0       PDE10A
2 1.73e-17       1.39 0.977 0.606  4.21e-13 0       TRPM3 
3 3.73e- 7       1.94 0.708 0.673  9.09e- 3 1       CENPF 
4 1.16e- 4       1.81 0.354 0.17   1   e+ 0 1       ANKRD1
5 2.35e-22       2.34 0.971 0.24   5.73e-18 2       SRRM4 
6 2.82e-25       2.16 0.912 0.117  6.88e-21 2       SRRM3 
> table(x@clone.composition$v1$clone.id)

 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
 3  3  2  2  3  2  3  3 27  2  3  3 13  2  2  3  2  2  2  2  2  2  2  2  3  2 
27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 
 3  2  4  3  3  5  2  2  2  2  2  2  2  2  2  7  2  2  2  3  2  2  2  2  2  2 
53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 
 2  3  2  2  2  2  2  4  2  2  3  2  2  2  3  2  2  2  2  2  2  2  2  2  2  2 
79 80 81 82 83 
 2  2  2  2  2 
> cluster9_bc = x@clone.composition$v1$cell.barcode[x@clone.composition$v1$clone.id == 9]
> length(cluster9_bc)
[1] 27
> cluster13_bc = x@clone.composition$v1$cell.barcode[x@clone.composition$v1$clone.id == 13]
> pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
# A tibble: 15 × 7
# Groups:   cluster [3]
      p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
      <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
 1 5.50e-22       1.68 0.992 0.737  1.34e-17 0       PDE10A  
 2 1.73e-17       1.39 0.977 0.606  4.21e-13 0       TRPM3   
 3 3.45e-10       1.33 0.443 0.071  8.40e- 6 0       EFS     
 4 3.85e- 4       1.31 0.496 0.253  1   e+ 0 0       MYH7B   
 5 3.13e- 7       1.26 0.672 0.354  7.61e- 3 0       ADAMTS18
 6 3.73e- 7       1.94 0.708 0.673  9.09e- 3 1       CENPF   
 7 1.16e- 4       1.81 0.354 0.17   1   e+ 0 1       ANKRD1  
 8 5.86e-18       1.73 0.954 0.964  1.43e-13 1       TPM1    
 9 6.84e- 3       1.53 0.277 0.158  1   e+ 0 1       CLSPN   
10 5.51e- 4       1.39 0.569 0.527  1   e+ 0 1       COL2A1  
11 2.35e-22       2.34 0.971 0.24   5.73e-18 2       SRRM4   
12 2.82e-25       2.16 0.912 0.117  6.88e-21 2       SRRM3   
13 4.97e-18       2.10 0.853 0.194  1.21e-13 2       GRIA2   
14 3.26e-17       1.89 0.618 0.071  7.94e-13 2       GRIK3   
15 3.09e-13       1.86 0.647 0.128  7.54e- 9 2       TH      
> Idents(pbmc) = "Not interested"
> pbmc = SetIdent(object = pbmc, cells = cluster13_bc, value = "cluster13")
> pbmc = SetIdent(object = pbmc, cells = cluster9_bc, value = "cluster9")
> Idents(pbmc)
CATACAGCATGGTGGA-1 CTCACTGTCCGCTAGG-1 AACAAAGCAAAGTGTA-1 GTGTCCTTCGTCCTTG-1 
    Not interested     Not interested     Not interested     Not interested 
TCCGGGACAACATCGT-1 AACCTTTTCTGTGCAA-1 AATCACGTCGGAATTC-1 AACGAAAGTAGCCCTG-1 
    Not interested     Not interested     Not interested     Not interested 
TAATCTCGTATCGTGT-1 AACGGGAAGTCTGGAG-1 GTTACCCGTACCCGAC-1 TCTCACGTCGTTCCTG-1 
    Not interested     Not interested     Not interested     Not interested 
AAGACTCAGGGATGTC-1 CCGGTGACAGTCCCGA-1 AAGTCGTAGCATAGGC-1 CTTACCGGTAAGCTCT-1 
    Not interested     Not interested     Not interested     Not interested 
GACCCTTAGCCGGAAT-1 AATGGAACAGCACACC-1 ATTACCTGTATGGGAC-1 TAGGGTTGTTACCCTC-1 
    Not interested     Not interested     Not interested     Not interested 
AATTCCTCATCCGAAT-1 ACTGTGAGTGGTATGG-1 ACTTATCAGATGCAGC-1 AGACTCATCTCTATGT-1 
    Not interested           cluster9           cluster9           cluster9 
AGATCCATCCGATTAG-1 ATGAGTCCAGGCTTGC-1 CAACCAAGTTCCAGGC-1 CAGGTATCACCAATTG-1 
          cluster9           cluster9           cluster9           cluster9 
CAGTTCCAGTCACGCC-1 CATCAAGAGGCTCCCA-1 CATTGTTTCGGAGATG-1 CCTCAGTAGCTGTTAC-1 
          cluster9           cluster9           cluster9           cluster9 
CCTCTAGCAGTATGAA-1 CGTTCTGCAACCCTCT-1 GATGATCCATTCCTAT-1 GCTGGGTAGACTGGGT-1 
          cluster9           cluster9           cluster9           cluster9 
GCTGGGTAGTCTAGCT-1 GGAACCCAGTGTTGTC-1 GGTGTTAGTTACCTGA-1 GGTTCTCAGGACAACC-1 
          cluster9           cluster9           cluster9           cluster9 
GTTGTAGCACTGGAAG-1 TATCGCCTCGTCTCAC-1 TGAATGCGTTACCCTC-1 TGACTCCGTCGGCACT-1 
          cluster9           cluster9           cluster9           cluster9 
TGATTTCTCAATCTCT-1 TTAATCCTCAAACGAA-1 TTACTGTTCCAGCCTT-1 ACAAGCTTCATGGAGG-1 
          cluster9           cluster9           cluster9           cluster9 
TAGTGCACAGACCCGT-1 ACACGCGCAGGTAGTG-1 ATACCTTGTGCTCTTC-1 TCTACCGTCAGGTAAA-1 
    Not interested     Not interested     Not interested     Not interested 
ACATCGAAGTGACCTT-1 CGAAGGATCTAGTGTG-1 TTGACCCTCGTAACTG-1 ACATGCAAGTAACCGG-1 
    Not interested     Not interested     Not interested     Not interested 
AGTACTGCATTCTTCA-1 AGTCAACGTATGGAGC-1 AGTGTTGCAGTCGTTA-1 CACCGTTGTTGGGCCT-1 
         cluster13          cluster13          cluster13          cluster13 
CGGGTGTAGTGGATAT-1 GACTTCCGTGCCCAGT-1 GGCTTGGGTAAGACCG-1 GGGCTACAGGTAGTAT-1 
         cluster13          cluster13          cluster13          cluster13 
GTCTACCTCCTTTAGT-1 GTGGTTATCCCGAGTG-1 TCAGTTTGTCAACACT-1 TGAGGTTAGACATACA-1 
         cluster13          cluster13          cluster13          cluster13 
ACCAACAAGCGTCTCG-1 GCCAGCACACATGACT-1 ACCAACAGTCGTACAT-1 GCCAGTGAGGATTCCT-1 
         cluster13     Not interested     Not interested     Not interested 
ACCCAAAGTTCTCTAT-1 CACGTGGGTTACGTAC-1 CATGCCTAGGAGTATT-1 ACCTGAACAGGAGGTT-1 
    Not interested     Not interested     Not interested     Not interested 
CTTCAATTCCCGAGGT-1 ACGATCAAGAGCTTTC-1 GGGAGATTCAAGCCTA-1 ACGGAAGCACGGTGAA-1 
    Not interested     Not interested     Not interested     Not interested 
ATCATTCTCGCTACAA-1 ACGTCCTCATGTACGT-1 TCAGGTACAGGGAATC-1 ACGTCCTGTTCTTGTT-1 
    Not interested     Not interested     Not interested     Not interested 
GAAGAATAGTTACGTC-1 ACTCCCATCATCGACA-1 TTGACCCCAAGCGCTC-1 ACTTAGGCAGAGAGGG-1 
    Not interested     Not interested     Not interested     Not interested 
TGAACGTTCACCATAG-1 AGAAGCGAGCCGTCGT-1 TCATGCCCAAGAAATC-1 AGACAGGAGAGATTCA-1 
    Not interested     Not interested     Not interested     Not interested 
CACTTCGTCGAGTGAG-1 GAAATGATCGTGGACC-1 AGACCATAGGGCAAGG-1 AGTGATCCAACCCGCA-1 
    Not interested     Not interested     Not interested     Not interested 
AGCCAATCATCCGGCA-1 GGGTTTATCGCAGTGC-1 TCCCACAAGGCTGTAG-1 AGCTTCCTCATTGCGA-1 
    Not interested     Not interested     Not interested     Not interested 
GCGGAAACATCACCAA-1 AGGACGATCCGTGGGT-1 CATCCACGTGTATTGC-1 GTAATGCAGTATGAGT-1 
    Not interested     Not interested     Not interested     Not interested 
TACCCGTAGCTAAATG-1 AGGACTTGTTGTGGAG-1 GAGTTACTCTGTGCTC-1 TTACGTTGTTGTGTAC-1 
    Not interested     Not interested     Not interested     Not interested 
AGGCTGCGTAGCGTTT-1 CAGCAGCTCAATCTTC-1 GTGTGATTCTGATGGT-1 AGGGTCCTCCATCTGC-1 
    Not interested     Not interested     Not interested     Not interested 
GTGAGCCGTAGGCAGT-1 TCCTCTTGTGCGGATA-1 TCCTGCACACTCAGAT-1 TGTGAGTCAAACCATC-1 
    Not interested     Not interested     Not interested     Not interested 
AGTAGTCAGAGAATCT-1 CATGCAAGTTGTAGCT-1 AGTCATGCACAAATGA-1 TACCTCGGTCTTGTCC-1 
    Not interested     Not interested     Not interested     Not interested 
AGTGCCGCAGCTATTG-1 GAATAGAAGCACCAGA-1 ATACTTCAGTACCCTA-1 GTCTACCCAAGCTCTA-1 
    Not interested     Not interested     Not interested     Not interested 
ATCAGGTCAGCTCTGG-1 ATTCCTACACACGCCA-1 ATGGTTGAGAGTTGTA-1 GTGGTTAGTCTTCATT-1 
    Not interested     Not interested     Not interested     Not interested 
ATTACCTAGTGTTCCA-1 CATTGAGCAGCAGGAT-1 ATTACCTCAGGGTCTC-1 TGTGGCGTCGTTTACT-1 
    Not interested     Not interested     Not interested     Not interested 
ATTACTCCACAAGCCC-1 TCGACGGTCGACGATT-1 ATTCATCAGACGGATC-1 CATCGGGAGTGGATTA-1 
    Not interested     Not interested     Not interested     Not interested 
GACCCAGTCAGACCTA-1 GGGACAACAGAGCTAG-1 GTAACCAGTAGTGATA-1 TGTCCACAGGGCAAGG-1 
    Not interested     Not interested     Not interested     Not interested 
TTAGGCAGTCAGACTT-1 ATTCGTTCAGCTGTGC-1 TGTACAGAGATCCGAG-1 CAACAACCATAGAATG-1 
    Not interested     Not interested     Not interested     Not interested 
TTTATGCAGGAATCGC-1 CACAACACACTGCTTC-1 GTAGGTTCAGTTGCGC-1 CACACAACACCGCTGA-1 
    Not interested     Not interested     Not interested     Not interested 
GGTGGCTGTGGCAGAT-1 TCGGATACAGTAGATA-1 CAGGTATGTATCCTCC-1 TATACCTGTACTCGTA-1 
    Not interested     Not interested     Not interested     Not interested 
CATCGGGAGACTCATC-1 TGTAACGAGCCGATTT-1 CATTCTAGTAGAGTTA-1 TACGCTCGTGCCTATA-1 
    Not interested     Not interested     Not interested     Not interested 
CCCAACTCAAGATGGC-1 TGATTCTCAGCAGTAG-1 CCGAACGTCCAAGGGA-1 GCATGATTCCGTAGTA-1 
    Not interested     Not interested     Not interested     Not interested 
CCTCATGTCGGATTAC-1 GATCAGTTCCAACCAA-1 CGATCGGGTCCAACGC-1 CTACAGAGTGCCTAAT-1 
    Not interested     Not interested     Not interested     Not interested 
CGCAGGTTCAGACCCG-1 GGGACAACAGTCAGCC-1 GGGCGTTCAACGATCT-1 CGCATGGCAACGTAAA-1 
    Not interested     Not interested     Not interested     Not interested 
TCGACCTTCATTTGGG-1 CGTGATACATGGGCAA-1 CTGAGGCTCCATCAGA-1 CGTGCTTCACTTGAAC-1 
    Not interested     Not interested     Not interested     Not interested 
TCATCATGTTAGTCGT-1 CGTGTCTGTCGCTTAA-1 CTGTACCTCTTCGGAA-1 CTAACCCCAAAGGATT-1 
    Not interested     Not interested     Not interested     Not interested 
TACGGTATCCGTGGCA-1 CTACGGGAGTTGCCTA-1 GGTGTCGCAAACGAGC-1 TAGTGCAGTGTCTTCC-1 
    Not interested     Not interested     Not interested     Not interested 
TCCCACATCTTGAACG-1 CTCAACCCATTCACAG-1 GTAGAGGAGCCATATC-1 CTCAAGACAGCAAGAC-1 
    Not interested     Not interested     Not interested     Not interested 
GTGGGAATCCAAGCTA-1 CTGAATGAGGTCCAGA-1 CTTAGGATCTTTCCGG-1 GAAACCTTCAAGGAGC-1 
    Not interested     Not interested     Not interested     Not interested 
CTGCTCAAGAGTGACC-1 GCATCGGGTTCCTTGC-1 CTGTACCTCCAGTACA-1 TGGATCATCGACCTAA-1 
    Not interested     Not interested     Not interested     Not interested 
CTGTAGACAAGTTCCA-1 TTTAGTCGTTCTTAGG-1 CTTCAATGTGTGTGGA-1 TCGCACTGTATTCCGA-1 
    Not interested     Not interested     Not interested     Not interested 
TTGGGCGCACTTGAAC-1 CTTCGGTAGTGAGTTA-1 TTTACGTAGAGTGTGC-1 GAGTTTGCACCTGTCT-1 
    Not interested     Not interested     Not interested     Not interested 
TTTGACTTCTTAATCC-1 GCAGCCATCGACGTCG-1 TCCTCCCCAGGCACTC-1 GCATTAGGTGTCGCTG-1 
    Not interested     Not interested     Not interested     Not interested 
GCTACAAGTCCAGAAG-1 GCATTAGTCCGAGGCT-1 TCAGTCCCAGTTGGTT-1 GCTTCACGTTCGTGCG-1 
    Not interested     Not interested     Not interested     Not interested 
GGGTGTCCACGGTGAA-1 GGGACCTGTCTACGTA-1 TGGTTAGCAAGGTTGG-1 GGGACCTTCAGCAATC-1 
    Not interested     Not interested     Not interested     Not interested 
GTGATGTTCGAGTCTA-1 GGTCACGAGGTAGCCA-1 TGTCAGATCTGCTAGA-1 GGTTCTCTCATAAGGA-1 
    Not interested     Not interested     Not interested     Not interested 
GTGTTCCGTGTGGTCC-1 GTAACACAGTGTACAA-1 TAGACCATCGCGTTTC-1 GTCCTCACACAGCTGC-1 
    Not interested     Not interested     Not interested     Not interested 
TCTTGCGTCCTTCAGC-1 GTGCAGCGTTCCTTGC-1 TTGAGTGTCAACACCA-1 GTGTTCCGTGGCAACA-1 
    Not interested     Not interested     Not interested     Not interested 
TGTGGCGGTAACCCTA-1 TCAGCAAAGTAACAGT-1 TTAGGCAAGTCAGGGT-1 TCTATACAGGGTGAGG-1 
    Not interested     Not interested     Not interested     Not interested 
TTGGTTTAGCAGTAAT-1 TGCACGGTCATAAGGA-1 
    Not interested     Not interested 
Levels: cluster9 cluster13 Not interested
> pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
Calculating cluster cluster9
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
Calculating cluster cluster13
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
Calculating cluster Not interested
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02s  
> pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
# A tibble: 15 × 7
# Groups:   cluster [3]
       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster        gene          
       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>          <chr>         
 1 0.0000318      1.21  0.333 0.074     0.774 cluster9       CHRND         
 2 0.00441        1.14  0.444 0.222     1     cluster9       PARP14        
 3 0.00438        0.942 0.815 0.567     1     cluster9       AUXG01000058.1
 4 0.000652       0.924 0.333 0.103     1     cluster9       RBL1          
 5 0.0000493      0.881 0.407 0.118     1     cluster9       TATDN2        
 6 0.00238        1.84  0.308 0.069     1     cluster13      RAB13         
 7 0.000109       1.45  0.615 0.194     1     cluster13      USP7          
 8 0.00279        1.40  0.846 0.429     1     cluster13      FLNB          
 9 0.000345       1.34  0.462 0.115     1     cluster13      SUPV3L1       
10 0.000393       1.34  0.385 0.083     1     cluster13      DNAH11        
11 0.00794        0.821 0.532 0.375     1     Not interested RPLP0         
12 0.00921       -0.268 0.121 0.3       1     Not interested NYNRIN        
13 0.00948       -0.288 0.126 0.3       1     Not interested NAE1          
14 0.00838       -0.329 0.105 0.275     1     Not interested DOC2B         
15 0.00783       -0.334 0.189 0.4       1     Not interested WDR19         
> cluster2.markers <- FindMarkers(pbmc, ident.1 = "cluster9", min.pct = 0.25)
head(cluster2.markers, n = 5)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
              p_val avg_log2FC pct.1 pct.2    p_val_adj
FSIP1  1.044809e-08  0.7667346 0.296 0.025 0.0002544529
HDHD5  4.052633e-06  0.6962610 0.333 0.059 0.0986978177
CHRND  3.178771e-05  1.2060144 0.333 0.074 0.7741578419
KAT5   4.304879e-05  0.6281316 0.259 0.044 1.0000000000
TATDN2 4.926712e-05  0.8808845 0.407 0.118 1.0000000000
> cluster2.markers <- FindMarkers(pbmc, ident.1 = "cluster13", min.pct = 0.25)
head(cluster2.markers, n = 5)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04s  
                 p_val avg_log2FC pct.1 pct.2   p_val_adj
MPND      2.821349e-07  0.9634384 0.308 0.023 0.006871112
WDR1      5.285807e-07  0.7723548 0.462 0.060 0.012873055
NDUFB8    5.035966e-06  0.7393046 0.385 0.051 0.122645917
AZIN1-AS1 5.274712e-06  1.0295656 0.385 0.051 0.128460340
COPRS     2.625170e-05  0.5394031 0.308 0.037 0.639333813
> cluster2.markers <- FindMarkers(pbmc, ident.1 = "cluster9",indent.2 = "cluster13", min.pct = 0.25)
head(cluster2.markers, n = 5)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
              p_val avg_log2FC pct.1 pct.2    p_val_adj
FSIP1  1.044809e-08  0.7667346 0.296 0.025 0.0002544529
HDHD5  4.052633e-06  0.6962610 0.333 0.059 0.0986978177
CHRND  3.178771e-05  1.2060144 0.333 0.074 0.7741578419
KAT5   4.304879e-05  0.6281316 0.259 0.044 1.0000000000
TATDN2 4.926712e-05  0.8808845 0.407 0.118 1.0000000000
> cluster2.markers <- FindMarkers(pbmc, ident.1 = "cluster9",indent.2 = "cluster13", min.pct = 0.25, logfc.threshold = 0.25)
head(cluster2.markers, n = 5)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
              p_val avg_log2FC pct.1 pct.2    p_val_adj
FSIP1  1.044809e-08  0.7667346 0.296 0.025 0.0002544529
HDHD5  4.052633e-06  0.6962610 0.333 0.059 0.0986978177
CHRND  3.178771e-05  1.2060144 0.333 0.074 0.7741578419
KAT5   4.304879e-05  0.6281316 0.259 0.044 1.0000000000
TATDN2 4.926712e-05  0.8808845 0.407 0.118 1.0000000000
> cluster2.markers <- FindMarkers(pbmc, ident.1 = "cluster13",indent.2 = "cluster9", min.pct = 0.25, logfc.threshold = 0.25)
head(cluster2.markers, n = 5)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04s  
                 p_val avg_log2FC pct.1 pct.2   p_val_adj
MPND      2.821349e-07  0.9634384 0.308 0.023 0.006871112
WDR1      5.285807e-07  0.7723548 0.462 0.060 0.012873055
NDUFB8    5.035966e-06  0.7393046 0.385 0.051 0.122645917
AZIN1-AS1 5.274712e-06  1.0295656 0.385 0.051 0.128460340
COPRS     2.625170e-05  0.5394031 0.308 0.037 0.639333813