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
bam_data="data/bam/"
out="data/out/"
collapsing_name="collapsing.txt"
save_progress_name="test"
whitelist_version="v1"
whitelist_path="data/whitelist/V1.CellTag.Whitelist.csv"
high_filter=20
low_filter=1
tagged=1
bamfilter=false
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
        ref_name=$2
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
    --whitelist_version)
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

echo "gunzip data/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" >> data/Log.log
gunzip data/$sample_name/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

if $bamfilter; then
echo 'Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version --bamfilter' >> data/Log.log
Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version --bamfilter
else
echo 'Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version' >> data/Log.log
Rscript --vanilla scripts/cell.R --out $out --bam_data $bam_data --collapsing_name $collapsing_name --save_progress_name $save_progress_name --whitelist_version $whitelist_version
fi

echo 'starcode -s --print-clusters $out${save_progress_name}_$collapsing_name > $out/${save_progress_name}_${collapsing_name%.*}_result.txt' >> data/Log.log
starcode -s --print-clusters $out${save_progress_name}_$collapsing_name > $out/${save_progress_name}_${collapsing_name%.*}_result.txt

if $visualize; then
echo 'Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --high_filter $high_filter --low_filter $low_filter --tagged $tagged --visualize' >> data/Log.log
Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --high_filter $high_filter --low_filter $low_filter --tagged $tagged --visualize
else
echo 'Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --high_filter $high_filter --low_filter $low_filter --tagged $tagged' >> data/Log.log
Rscript --vanilla scripts/filtering_cloneCalling.R --out $out --collapsing_name "${save_progress_name}_${collapsing_name%.*}_result.txt" --save_progress_name $save_progress_name --whitelist_path $whitelist_path --high_filter $high_filter --low_filter $low_filter --tagged $tagged
fi