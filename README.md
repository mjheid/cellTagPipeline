# cellTag Pipeline
In this Project we implement a command line tool for the CellTag pipline using 
[CellTagR](https://github.com/morris-lab/CellTagR). Most necessary steps for running 
the bioinformatics celltag pipline can be done by this tool. We also provide a simple 
way to create data visualisations of the pipeline, to create a first look into the 
data to analyse.


## Instalation


It is recommended to create a conda envrionment to install all needed packages. Do 
as follows:
```sh
conda create -n celltag  
conda activate celltag  
conda install -c bioconda r-seurat r-optparse r-devtools starcode r-ggplot2  
```
It might take some time for conda to find compatible verions to install everything. 
To complete the installation open the R command line of the celltag environment and 
do the following:
```R
library("devtools")  
devtools::install_github("morris-lab/CellTagR")  
```
Again, this might take quite some time to install.


## Structure

```
.gitignore  
README.md  
analysis.Rmd  
analysis.html  
control.sh  
data  
   |-- CellTag.UTR  
   |-- GFP.CDS  
   |-- gene_list.txt
   |-- bam  
   |   |-- .gitinclude  
   |-- fastq  
   |   |-- .gitinclude  
   |-- out  
   |   |-- .gitinclude  
   |-- refgenomes  
   |   |-- .gitinclude  
   |-- samples  
   |   |-- .gitinclude  
   |-- sourceFile  
   |   |-- .gitinclude  
   |-- whitelist  
   |   |-- V1.CellTag.Whitelist.csv  
   |   |-- V2.CellTag.Whitelist.csv  
   |   |-- V3.CellTag.Whitelist.csv  
p2.sh  
pipe.sh  
scripts  
   |-- cell.R  
   |-- cellTag.sh  
   |-- filtering_cloneCalling.R  
   |-- helperFunctions.R  
   |-- insert_data.py  
   |-- map.py  
   |-- map_cellTags.R  
   |-- network.R  
   |-- preprocessing.sh  
```

The folder `/scripts` contains all scripts that are used by the pipeline and are called over 
`control.sh`. `scripts/helperFunctions.R` contains functions for plotting and printing of information 
used by other scripts. `data/` contains data needed and produced by the pipeline. `data/CellTag.UTR` is a fasta 
file of the CellTag.UTR transposon, `data/GFP.CDS` is a fasta file of the GFP.CDS transposon. `data/gene_list.txt`
is a file containing row seperated gene names. These are used when visualising gene expression. `data/bam/`
 contains `.bam` files needed to create CellTag UMI count matrixes. `data/fastq/` contains fastq
files , which contains the raw sequenced data. `data/outt` contains output files produced during the 
execution of the pipeline. `data/refgenome/` contains the reference genome output created from 
cellranger mkref. `data/samples/` contains the output from cellranger count. It is assumed that 
gene count data is found here, as if cellranger count had this dir as an output. `data/sourceFile/`
contains `.fa` and `.gtf` files used to build the reference genome. They should enclude CellTag.UTR 
and GFP.CDS. `data/whitelist/` contains the three whitelists of the three CellTag librarys. 
`scripts/preprocessing.sh` contains the preprocessing part of the pipeline, expects user to have 
`samtools/1.16.1` and `cellranger/5.0.1 `loadable as folows: `module load [tool]`. `scripts/cellTag.sh` 
contains the CellTag UMI count matrix generation and filtering, as well as clone calling. This gets 
executed by firstly calling `scripts/cell.R`, loading `.bam`data and preparing CellTags for collapsing, 
followed by `scripts/filtering_cloneCalling.R` in which CellTags are filtered and clones called. 
`scripts/network.R` creates a clone network for data in which all libraries have been clone called. 
`scripts/map_cellTags.R` runs the seurat [workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) 
while [regressing out cell cycle genes](https://satijalab.org/seurat/articles/cell_cycle_vignette.html) as 
well as [mitochrondial genes](https://github.com/satijalab/seurat/issues/6613). Afterwards several visualisation 
options exist to visualize distribution of CellTags, clones and genes. `analysis.Rmd` creates `analysis.html`,
showcaseing some gene expression analysis of CellTag data. `control.sh` is the script over which the created 
tool can be excuted and controlled.

## control of CellTag Pipeline

The control center of the Pipeline. Executed as follows:

```sh
./control.sh [command] [arguments]
```

Type `./control.sh --help` for the following help message on how to use the control script:

```
Usage: ./control.sh [command] [arguments]  
Options:  
  preprocess [arguments]    Run the preprocessing script with the provided arguments.  
                            Build a reference genome, count stuff, filter stuff with samtools.  
  cellTag [arguments]    Run the cellTagR pipeline with the provided arguments  
  network [arguments]    Visualize clone lineages through networks and more.  
  map [arguments]    Visualize cellTags in dimension reduction plots, TSNE and UMAP supported.  
  --help    Print this help. For more specific help on commands type [command] --help  
```

## Preprocessing

Executed via `./control.sh preprocess [arguments`. List of arguments can be seen via `./control.sh preprocess --help`:

```
Usage: scripts/preprocessing.sh [arguments]
Arguments:
  --mkref    Run cellranger mkref. Default: false
  --count    Run cellranger count. Default: false
  --filter_sam    Run samtools to filter .bam provided by cellranger count.
                  Default: false
  --genome [Path]   Path from current dir to fasta file of genome. Required 
                    for --mkref. Default /sourceFile/genome_new.fa
  --genes [Path]   Path from current dir to .gtf file of genes of genome. 
                   Required for --mkref. Default: /sourceFile/genes_new.gtf
  --refgenome [String]   Name of reference genome. Needed for --mkref as output 
                         dir name and --count as input. Default: refdata-gex-GRCh38-2020-A-cellTag-barcode
  --sample_name [String]   Name of Sample, output from ilumina. Used for --count
                           as well as --filter_sam. Default: SI-TT-H4
  --fastqs [Path]   Path from current dir to fastqs of ilumina run. Contains 
                    --sample_name samples. Used for --count. Default: None
  --expect_cells [int]   Number of expected cells for --count. Default: 5000
  --GFP [Path]   Path from current dir to GFP.CDS file, for filtering for 
                 reads mapped to this sequence. Used with --filter_sam. Default: data/GFP.CDS
  --UTR [Path]   Path from current dir to CellTag.UTR file, for filtering for
  --GFP [Path]   reads mapped to this sequence. Used with --filter_sam. Default data/CellTag.UTR
  --bam_data [Path]   Path from current dir to output dir of filtered .bam 
                      files when using --filter_sam. Default: data/bam
  --help    Print this help
```

Output when using '--mkref' will be in 'data/refgenomes/[--refgenome]/'.
Output when using '--count' will be in 'data/samples/[--sample_name]/'.
Example execution:

```sh
./control.sh preprocess --mkref --count --filter_sam --fastqs "data/fastq/"
```

## CellTag Pipline

Executed via `./control.sh cellTag [arguments]`. List of arguments can be seen via `./control.sh cellTag --help`:

```
Usage: scripts/cellTag.sh [arguments]  
Arguments:  
  --visualize    Visualize filtering steps. Default: false  
  --out [Path]   Path from current dir to dir where produced files  
                 should get saved. Default: data/out/  
  --collapsing_name [String]   Name of collapsing file, used for  
                               starcode. Default: collapsing.txt  
  --save_progress_name [String]   Name of outputs. Default: test  
  --sample_name [String]   Name of Sample, output from ilumina.  
                           Used to find path barcodes.tsv. Default: SI-TT-H4  
  --whitelist_version [String]   Version of the whitelist. Either v1, v2 or v3.  
                                 Default: v1  
  --whitelist_path [Path]   Path from current dir to the whitelist csv, gets adujusted.  
                            based on specified version. Default: data/whitelist/V1.CellTag.Whitelist.csv  
  --high_filter [int]   Upper bound for filtering celltags. Default: 20  
  --low_filter [int]   Lower bound for filtering celltags. Default: 2  
  --bam_data [Path]   Path from current dir to output dir of filtered .bam  
                      files when using --filter_sam. Default: data/bam/possorted_genome_bam.filtered.bam  
  --tagged [int]   When binarizing, the number of cellTags needed to  
                   consider cell as tagged. Default: 2  
  --jac_cut [float]   When calculating simiarity of cells voa Jaccard coefficiants,  
                      this is the cut-off used to consider cells for clone calling. Default: 0.7  
  --help    Print this help  
```

This command executes two other R script, `script/cell.R` and `scripts/filtering_cloneCalling.R`.
Between their execution starcode gets executed. If users want to add arguments, they should 
add arguments in `scripts/cellTag.sh`, the script in which that argument is used, as well as in 
the call of the script inside of 'scripts/cellTag.sh'. Also, make sure that all paths specified exist.  
Several outputfiles will be produced: a [whitelist_version][save_progress_name].RDS, containing 
current working matrixes, [save_progress_name].RDS which will accumalte data of different 
library versions as long as [save_progress_name] is the same during runs. Each run will produce 
multiple plots as well.
A possible way to execute this part of the script for all libaries could be as follows:  
```sh
./control.sh cellTag --visualize --save_progress_name "ct" --tagged 1 --low_filter 1 --bam_data "data/bam/possorted_genome_bam.filtered_merged.bam"  --whitelist_version "v1" --sample_name "SI-TT-H4" --out "data/out/cellTag/"  
#./control.sh cellTag --visualize --save_progress_name "ct" --tagged 1 --low_filter 1 --bam_data "data/bam/possorted_genome_bam.filtered_merged.bam"  --whitelist_version "v2" --sample_name "SI-TT-H4" --out "data/out/cellTag/"  
#./control.sh cellTag --visualize --save_progress_name "ct" --tagged 1 --low_filter 1 --bam_data "data/bam/possorted_genome_bam.filtered_merged.bam"  --whitelist_version "v3" --sample_name "SI-TT-H4" --out "data/out/cellTag/"  
```

## Network Construction

Executed via `./control.sh network [arguments]`. List of arguments can be seen via `./control.sh network --help`:  

```
Usage: scripts/network.R [options]  
Options:  
    --visualize  
        Visualize filtering. Default: false  
    --out=CHARACTER  
        output file dir. Default: data/out/  
    --start_node=CHARACTER  
        Start node to build Graph from. Default: BiggestNode  
    --save_progress_name=CHARACTER  
        Name appended to outputs to be saved. Default: test  
    -h, --help  
        Show this help message and exit  
```
`--start_node` if not specified finds the node with the highest amount of connections and uses 
that one to start building the network. If wanted you can take a look at nodes in the saved `.RDS` file as follows: `object@nodes` and `object@network.link.list`.  
This produces an `.html` file in which a clone network is visualized, as well as a plot which 
show how many cells contain CellTags from each library.  
A possible way to execute this would be as follows:  
```sh
./control.sh network --visualize --save_progress_name "ct"
```

## Mapping 

```
Usage: scripts/map_cellTags.R [options]  
Options:  
    --visualize_umap  
        Visualize filtering using umap. Default: false  
    --visualize_tsne  
        Visualize filtering using tsne. Default: false  
    --runUMAP  
        Run UMAP. Default: false  
    --runTSNE  
        Run TSNE. Default: false  
    --filter  
        Seurat workflow, normalizing, scaling, calculating PCs, calculating clusters,  
        regressing cell cycle genes and mitochrondrial genes out. Default: false  
    --jackstraw  
        Use jackstraw to confirm PCs p-value, otherwise elbow plot. Default: false  
    --visualize_ver=CHARACTER  
        visualize specified version, either v1, v2, v3 or all. Default: all  
    --out=CHARACTER  
        output file dir. Default: data/out/  
    --gene_list=CHARACTER  
        Gene list file location and name. Contains list of genes. Default: empty  
    --vis_clone=INT  
        Clone to visualize. Default: 1  
    --min_cells=INT  
        Minimum cell count to load data. Default: 3  
    --min_features=INT  
        Minimum features count to load data. Default: 200  
    --scale_factor=INT  
        Factor to scale cells with. Default: 10000  
    --npcs=INT  
        Number of PCs to calculate. Default: 50  
    --n_var_features=INT  
        Number of varible features to work with. Default: 2000  
    --neighbours_dims=INT  
        Number of PCAs for FindNeighbours to look at. Default: 10  
    --cluster_resolution=FLOAT  
        Resolution of clusters to calculate. Higher values find more clusters. Default: 0.5  
    --save_progress_name=CHARACTER  
        Name appended to outputs to be saved. Default: test  
    --sample_name=CHARACTER  
        Name of experiment sample. Default: SI-TT-H4  
    -h, --help  
        Show this help message and exit  
```

It is recommended to first filter, utalising the seurat workflow. Several quality control plots 
will be plotted, it is advised to take them into considertion and rerun with `--filter` if deemed 
necessary. Whenever TSNE or UMAP gets calculated, the seurat object `[save_progress_name]_reduction.RDS` 

gets saved. Whenever filter gets used  this object is saved as well. When choosing to visualize several plots 
get plotted: Location of CellTags, identified clones and specified clone on the 2D embedding. If gene list 
is specified, gene expression in all cells, as well as cells with ellTags, identified clones and specified clone 
get visualized. This is done for the specified library version. If "all" is specified as library version, 
all librarys get plotted.

A possible execution might look as follows:

```sh
./control.sh map --save_progress_name "ct" --sample_name "SI-TT-H4" --min_cells 3 --min_features 200 --filter --out "data/out/cellTag/"  
./control.sh map --visualize_umap --runUMAP --visualize_ver "v1" --save_progress_name "ct" --gene_list "data/gene_list.txt" --out "data/out/cellTag/" --vis_clone 9  
./control.sh map --visualize_umap --visualize_ver "v1" --save_progress_name "ct" --gene_list "data/gene_list.txt" --out "data/out/cellTag/" --vis_clone 13  
```
