
#./control.sh preprocess --filter_sam
#./control.sh cellTag --visualize --save_progress_name "test"
#./control.sh cellTag --visualize --save_progress_name "test_low2" --tagged 2 --low_filter 2
#./control.sh cellTag --visualize --save_progress_name "test_merged" --bamfilter
#./control.sh cellTag --visualize --save_progress_name "test_merged_low2" --tagged 2 --low_filter 2 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag1_merged_low2" --tagged 1 --low_filter 2 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag2_merged_low1" --tagged 2 --low_filter 1 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag1_low2" --tagged 1 --low_filter 2
#./control.sh cellTag --visualize --save_progress_name "tag2_low1" --tagged 2 --low_filter 1
./control.sh cellTag --visualize --save_progress_name "og_b1_l1_h20" --tagged 1 --low_filter 1 --bam_data "data/bam/hf1.d15.filtered.bam" --sample_name "og"
./control.sh cellTag --visualize --save_progress_name "og_b2_l1_h20" --tagged 2 --low_filter 1 --bam_data "data/bam/hf1.d15.filtered.bam" --sample_name "og"
./control.sh cellTag --visualize --save_progress_name "og_b1_l2_h20" --tagged 1 --low_filter 2 --bam_data "data/bam/hf1.d15.filtered.bam" --sample_name "og"
./control.sh cellTag --visualize --save_progress_name "og_b2_l2_h20" --tagged 2 --low_filter 2 --bam_data "data/bam/hf1.d15.filtered.bam" --sample_name "og"