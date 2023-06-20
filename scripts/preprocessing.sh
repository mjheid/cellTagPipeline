#!/bin/bash
# cat genes.gtf Anika.genes.gtf >> genes_new.gtf 
# cat genome.fa Anika.genome.fa >> genome_new.fa

# Define usage function
function usage {
  echo "Usage: $0 [arguments]"
  echo "Arguments:"
  echo "  --mkref    Run cellranger mkref"
  echo "  --count    Run cellranger count"
  echo "  --filter_sam    Run samtools to filter .bam provided by cellranger count"
  echo "  --genome [Path]   Path from current dir to fasta file of genome. Required for --mkref"
  echo "  --genes [Path]   Path from current dir to .gtf file of genes of genome. Required for --mkref"
  echo "  --refgenome [String]   Name of reference genome. Needed for --mkref as output dir name and --count as input"
  echo "  --sample_name [String]   Name of Sample, output from ilumina. Used for --count as well as --filter_sam"
  echo "  --fastqs [Path]   Path from current dir to fastqs of ilumina run. Contains --sample_name samples. Used for --count"
  echo "  --expect_cells [int]   Number of expected cells for --count"
  echo "  --GFP [Path]   Path from current dir to GFP.CDS file, for filtering for reads mapped to this sequence. Used with --filter_sam"
  echo "  --UTR [Path]   Path from current dir to CellTag.UTR file, for filtering for reads mapped to this sequence. Used with --filter_sam"
  echo "  --bam_data [Path]   Path from current dir to output dir of filtered .bam files when using --filter_sam"
  echo "  --help    Print this help"
  exit 1
}

# Initialize variables
mkref=false
count=false
filter_sam=false
refgenome="refdata-gex-GRCh38-2020-A-cellTag-barcode"
genome="/sourceFile/genome_new.fa"
genes="/sourceFile/genes_new.gtf"
sample_name="SI-TT-H4"
expect_cells=5000
GFP="data/GFP.CDS"
UTR="data/CellTag.UTR"
bam_data="data/bam"

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "No Arguments provided"
    echo "Execute with --help for a list of arguments"
    exit 1
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --mkref)
        mkref=true
        shift
        ;;
    --count)
        count=true
        shift
        ;;
    --filter_sam)
        filter_sam=true
        shift
        ;;
    --genome)
        genome=$2
        shift
        shift
        ;;
    --genes)
        genes=$2
        shift
        shift
        ;;
    --refgenome)
        ref_name=$2
        shift
        shift
        ;;
    --sample_name)
        sample_name=$2
        shift
        shift
        ;;
    --fastqs)
        fastqs=$2
        shift
        shift
        ;;
    --expect_cells)
        expect_cells=$2
        shift
        shift
        ;;
    --GFP)
        GFP=$2
        shift
        shift
        ;;
    --UTR)
        UTR=$2
        shift
        shift
        ;;
    --bam_data)
        bam_data=$2
        shift
        shift
        ;;
    --help)
        usage
        ;;
    *)
        echo "Invalid option: $key"
        exit 1
        ;;
  esac
done


# Executing cellranger mkref, building a reference genome
if $mkref; then
    conda_env=$(conda info --envs | grep \* | cut -d " "  -f 1)
    conda deactivate
    conda deactivate
    
    genome="$(pwd)$genome"
    genes="$(pwd)$genes"
    
    cd data/refgenomes
    
    module load cellranger/5.0.1 
    
    echo "cellranger mkref --genome=$refgenome --fasta=$genome --genes=$genes"  >> data/Log.log
    cellranger mkref --genome=$refgenome --fasta=$genome --genes=$genes
    
    cd ../..
    conda activate $conda_env
    
fi


# Executing cellranger count
if $count; then

    conda_env=$(conda info --envs | grep \* | cut -d " "  -f 1)
    conda deactivate
    conda deactivate
    
    fastqs="$(pwd)$fastqs"
    refgenome="$(pwd)refgenomes/$refgenome"
    
    cd data/samples
    
    echo "cellranger count --id=${sample_name}  \
                    --fastqs=${fastqs}  \
                    --sample=${sample_name}  \
                    --transcriptome=${refgenome}  \
                    --localcores=40  \
                    --expect-cells=${expect_cells}  " >> data/Log.log
    cellranger count --id=${sample_name}  \
                    --fastqs=${fastqs}  \
                    --sample=${sample_name}  \
                    --transcriptome=${refgenome}  \
                    --localcores=40  \
                    --expect-cells=${expect_cells}  \
    
    cd ../..
    conda activate $conda_env
fi

# Filtering with samtools
if $filter_sam; then

    module load samtools/1.16.1
    
    echo "Filtering..."
    
    echo "samtools view -b -f 4 data/samples/$sample_name/outs/possorted_genome_bam.bam > $bam_data/possorted_genome_bam.filtered_umapped.bam &" >> data/Log.log
    samtools view -b -f 4 data/samples/$sample_name/outs/possorted_genome_bam.bam > $bam_data/possorted_genome_bam.filtered_umapped.bam &
    pidid1=$!
    
    echo "samtools view -b data/samples/$sample_name/outs/possorted_genome_bam.bam $GFP  > $bam_data/possorted_genome_bam.filtered_GFP.CDS.bam &" >> data/Log.log
    samtools view -b data/samples/$sample_name/outs/possorted_genome_bam.bam $GFP  > $bam_data/possorted_genome_bam.filtered_GFP.CDS.bam &
    pidid2=$!

    echo "samtools view -b data/samples/$sample_name/outs/possorted_genome_bam.bam $UTR  > $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam &" >> data/Log.log
    samtools view -b data/samples/$sample_name/outs/possorted_genome_bam.bam $UTR  > $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam &
    pidid3=$!

    wait $pidid1 $pidid2 $pidid3
    
    echo "Combining..."

    echo "cat $filter_out/possorted_genome_bam.filtered_umapped.bam $filter_out/possorted_genome_bam.filtered_GFP.CDS.bam $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam > $filter_out/possorted_genome_bam.filtered.bam" >> data/Log.log
    cat $filter_out/possorted_genome_bam.filtered_umapped.bam $filter_out/possorted_genome_bam.filtered_GFP.CDS.bam $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam > $filter_out/possorted_genome_bam.filtered.bam

    echo "time samtools merge  $bam_data/possorted_genome_bam.filtered_merged.bam $bam_data/possorted_genome_bam.filtered_umapped.bam  $bam_data/possorted_genome_bam.filtered_GFP.CDS.bam $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam" >> data/Log.log
    time samtools merge  $bam_data/possorted_genome_bam.filtered_merged.bam $bam_data/possorted_genome_bam.filtered_umapped.bam  $bam_data/possorted_genome_bam.filtered_GFP.CDS.bam $bam_data/possorted_genome_bam.filtered_CellTag.UTR.bam
    
    echo "Finished"

fi


