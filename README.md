# cellTag Pipeline
In this Project we implement a command line tool for the CellTag pipline using 
[CellTagR](https://github.com/morris-lab/CellTagR). Most necessary steps for running 
the bioinformatics celltag pipline can be done by this tool. We also provide a simple 
way to create data visualisations of the pipeline, to create a first look into the 
data to analyse.


## Instalation


It is recommended to create a conda envrionment to install all needed packages. Do 
as follows:
sh'''
conda create -n celltag
conda activate celltag
conda install -c bioconda r-seurat r-optparse r-devtools starcode r-ggplot2
'''
It might take some time for conda to find compatible verions to install everything. 
To complete the installation open the R command line of the celltag environment and 
do the following:
R'''
library("devtools")
devtools::install_github("morris-lab/CellTagR")
'''
Again, this might take quite some time to install.


## Structure

'''
cellTagPipeline/
.gitignore
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
pipe.sh
scripts
   |-- cell.R
   |-- cellTag.sh
   |-- filtering_cloneCalling.R
   |-- insert_data.py
   |-- map.py
   |-- map_cellTags.R
   |-- network.R
   |-- preprocessing.sh
   |-- test.R

'''

## control of CellTag Pipeline

The control center of the Pipeline. Executed as follows:

sh'''
./control.sh [command] [arguments]
'''

Type './control.sh --help' for the following help message on how to use the control script:

sh'''
Usage: ./control.sh [command] [arguments]
Options:
  preprocess [arguments]    Run the preprocessing script with the provided arguments.
                            Build a reference genome, count stuff, filter stuff with samtools.
  cellTag [arguments]    Run the cellTagR pipeline with the provided arguments
  network [arguments]    Visualize clone lineages through networks and more.
  map [arguments]    Visualize cellTags in dimension reduction plots, TSNE and UMAP supported.
  --help    Print this help. For more specific help on commands type [command] --help
'''

## Preprocessing

Executed via './control.sh preprocess [arguments]'. List of arguments can be seen via 
'./control.sh preprocess --help':

sh'''
Usage: scripts/preprocessing.sh [arguments]
Arguments:
  --mkref    Run cellranger mkref. Default: false
  --count    Run cellranger count. Default: false
  --filter_sam    Run samtools to filter .bam provided by cellranger count. 
                  Default: false
  --genome [Path]   Path from current dir to fasta file of genome. 
                    Required for --mkref. Default: "/sourceFile/genome_new.fa"
  --genes [Path]   Path from current dir to .gtf file of genes of genome. 
                   Required for --mkref. Default: "/sourceFile/genes_new.gtf"
  --refgenome [String]   Name of reference genome. Needed for --mkref as 
                         output dir name and --count as input
  --sample_name [String]   Name of Sample, output from ilumina. Used for 
                           --count as well as --filter_sam
  --fastqs [Path]   Path from current dir to fastqs of ilumina run. Contains 
                    --sample_name samples. Used for --count
  --expect_cells [int]   Number of expected cells for --count. Default: 5000.
  --GFP [Path]   Path from current dir to GFP.CDS file, for filtering for 
                 reads mapped to this sequence. Used with --filter_sam. Default:
                 "data/GFP.CDS"
  --UTR [Path]   Path from current dir to CellTag.UTR file, for filtering 
                 for reads mapped to this sequence. Used with --filter_sam.
                 Default: "data/CellTag.UTR"
  --bam_data [Path]   Path from current dir to output dir of filtered .bam 
                      files when using --filter_sam. Default: "data/bam"
  --help    Print this help
'''

Output when using '--mkref' will be in 'data/refgenomes/[--refgenome]/'.
Output when using '--count' will be in 'data/samples/[--sample_name]/'.
Example execution:

sh'''
./control.sh preprocess --mkref --count --filter_sam
'''

## CellTag Pipline

Executed via './control.sh cellTag [arguments]'. List of arguments can be seen via 
'./control.sh cellTag --help':

sh'''
Usage: scripts/cellTag.sh [arguments]
Arguments:
  --visualize    Visualize filtering steps. Default: false
  --out [Path]   Path from current dir to dir where produced files 
                 should get saved. Default: "data/out/"
  --collapsing_name [String]   Name of collapsing file, used for starcode.
                               Default: "collapsing.txt"
  --save_progress_name [String]   Name of outputs
  --sample_name [String]   Name of Sample, output from ilumina. Used to 
                           find path to barcodes.tsv
  --whitelist_version [String]   Version of the whitelist. Either "v1", "v2" or "v3".
  --whitelist_path [Path]   Path from current dir to the whitelist csv.
                            Only specifiy when using self created whitelist.
  --high_filter [int]   Upper bound for filtering celltags. Default: 20
  --low_filter [int]   Lower bound for filtering celltags. Default: 2
  --bam_data [Path]   Path from current dir to output dir of filtered .bam 
                      files when using --filter_sam
  --tagged [int]   When binarizing, the number of cellTags needed to consider 
                   cell as tagged. Default: 2
  --help    Print this help
'''

This command executes two other R script, 'script/cell.R' and 'scripts/filtering_cloneCalling.R'.
Between their execution starcode gets executed.
