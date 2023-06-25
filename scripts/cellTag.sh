#!/bin/bash

# Define usage function
function usage {
  echo "Usage: $0 [arguments]"
  echo "Arguments:"
  echo "  --visualize    Visualize filtering steps"
  echo "  --out [Path]   Path from current dir to dir where produced files should get saved"
  echo "  --collapsing_name [String]   Name of collapsing file, used for starcode"
  echo "  --save_progress_name [String]   Name of outputs"
  echo "  --sample_name [String]   Name of Sample, output from ilumina. Used to find path barcodes.tsv"
  echo "  --whitelist_version [String]   Version of the whitelist. Either v1, v2 or v3."
  echo "  --whitelist_path [Path]   Path from current dir to the whitelist csv"
  echo "  --high_filter [int]   Upper bound for filtering celltags."
  echo "  --low_filter [int]   Lower bound for filtering celltags"
  echo "  --bam_data [Path]   Path from current dir to output dir of filtered .bam files when using --filter_sam"
  echo "  --tagged [int]   When binarizing, the number of cellTags needed to consider cell as tagged"
  echo "  --help    Print this help"
  exit 1
}

# Initialize variables
sample_name="SI-TT-H4"
bam_data="data/bam/possorted_genome_bam.filtered.bam"
out="data/out/"
collapsing_name="collapsing.txt"
save_progress_name="test"
whitelist_version="v1"
whitelist_path="data/whitelist/V1.CellTag.Whitelist.csv"
whitelist_path1="data/whitelist/V1.CellTag.Whitelist.csv"
whitelist_path2="data/whitelist/V2.CellTag.Whitelist.csv"
whitelist_path3="data/whitelist/V3.CellTag.Whitelist.csv"
high_filter=20
low_filter=2
tagged=2
jac_cut=0.7
visualize=false

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
    --bam_data)
        bam_data=$2
        shift
        shift
        ;;
    --out)
        out=$2
        shift
        shift
        ;;
    --sample_name)
        sample_name=$2
        shift
        shift
        ;;
    --collapsing_name)
        collapsing_name=$2
        shift
        shift
        ;;
    --save_progress_name)
        save_progress_name=$2
        shift
        shift
        ;;
    --whitelist_version)
        whitelist_version=$2
        shift
        shift
        ;;
    --whitelist_path)
        whitelist_path=$2
        shift
        shift
        ;;
    --high_filter)
        high_filter=$2
        shift
        shift
        ;;
    --low_filter)
        low_filter=$2
        shift
        shift
        ;;
    --tagged)
        tagged=$2
        shift
        shift
        ;;
    --jac_cut)
        jac_cut=$2
        shift
        shift
        ;;
    --bamfilter)
        bamfilter=true
        shift
        ;;
    --visualize)
        visualize=true
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

# check if whitelist path got changed, if not, set whitelist 
# path in accordance with specified version
if [ "$whitelist_path" = "$whitelist_path1" ]; then
    if [ "$whitelist_version" = "v2" ]; then
        whitelist_path="$whitelist_path2"
    elif [ "$whitelist_version" = "v3" ]; then
        whitelist_path="$whitelist_path3"
    fi
fi

# gunzip barcodes if .gz
echo "gunzip data/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" >> data/Log.log
gunzip data/samples/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# Execute cell.R
echo "Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version --sample_name $sample_name" >> data/Log.log
Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version --sample_name $sample_name

# Annalyse which celltags to collapse
echo "starcode -s --print-clusters $out${whitelist_version}${save_progress_name}_$collapsing_name > $out${whitelist_version}${save_progress_name}_${collapsing_name%.*}_result.txt" >> data/Log.log
starcode -s --print-clusters $out${whitelist_version}${save_progress_name}_$collapsing_name > $out${whitelist_version}${save_progress_name}_${collapsing_name%.*}_result.txt

# Execute filtering_cloneCalling.R
if $visualize; then
echo "Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name '${save_progress_name}_${collapsing_name%.*}_result.txt' --save_progress_name $save_progress_name --whitelist_path $whitelist_path --whitelist_version $whitelist_version --high_filter $high_filter --low_filter $low_filter --tagged $tagged --visualize" >> data/Log.log
Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --whitelist_version $whitelist_version --high_filter $high_filter --low_filter $low_filter --tagged $tagged --visualize --jac_cut $jac_cut
else
echo "Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name '${save_progress_name}_${collapsing_name%.*}_result.txt' --save_progress_name $save_progress_name --whitelist_path $whitelist_path --whitelist_version $whitelist_version --high_filter $high_filter --low_filter $low_filter --tagged $tagged" >> data/Log.log
Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --whitelist_version $whitelist_version --high_filter $high_filter --low_filter $low_filter --tagged $tagged --jac_cut $jac_cut
fi

# gzip barcodes if .tsv
echo "gzip data/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv" >> data/Log.log
gzip data/samples/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv