Rscript ./extract_cellTags_fq.R

./filter_bam.sh

Rscript ./extract_cellTags_fq

./starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

Rscript ./name.R