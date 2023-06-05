#./control.sh preprocess --filter_sam
#./control.sh cellTag --visualize --save_progress_name "test"
#./control.sh cellTag --visualize --save_progress_name "test_low2" --tagged 2 --low_filter 2
#./control.sh cellTag --visualize --save_progress_name "test_merged" --bamfilter
#./control.sh cellTag --visualize --save_progress_name "test_merged_low2" --tagged 2 --low_filter 2 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag1_merged_low2" --tagged 1 --low_filter 2 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag2_merged_low1" --tagged 2 --low_filter 1 --bamfilter
#./control.sh cellTag --visualize --save_progress_name "tag1_low2" --tagged 1 --low_filter 2
#./control.sh cellTag --visualize --save_progress_name "tag2_low1" --tagged 2 --low_filter 1
#./control.sh cellTag --visualize --save_progress_name "big_b1_l1_h20" --tagged 1 --low_filter 1 --bam_data "data/test_data/hf1.d15.possorted_genome_bam.bam" --sample_name "og" --whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "big_b2_l1_h20" --tagged 2 --low_filter 1 --bam_data "data/test_data/hf1.d15.possorted_genome_bam.bam" --sample_name "og" #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "big_b1_l2_h20" --tagged 1 --low_filter 2 --bam_data "data/test_data/hf1.d15.possorted_genome_bam.bam" --sample_name "og" #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "big_b2_l2_h20" --tagged 2 --low_filter 2 --bam_data "data/test_data/hf1.d15.possorted_genome_bam.bam" --sample_name "og" #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "bigo_b1_l1_h20" --tagged 1 --low_filter 1 --bam_data "data/samples/SI-TT-H4/outs/possorted_genome_bam.bam"  --whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "bigo_b2_l1_h20" --tagged 2 --low_filter 1 --bam_data "data/samples/SI-TT-H4/outs/possorted_genome_bam.bam"  #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "bigo_b1_l2_h20" --tagged 1 --low_filter 2 --bam_data "data/samples/SI-TT-H4/outs/possorted_genome_bam.bam"  #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
#./control.sh cellTag --visualize --save_progress_name "bigo_b2_l2_h20" --tagged 2 --low_filter 2 --bam_data "data/samples/SI-TT-H4/outs/possorted_genome_bam.bam"  #--whitelist_version "v2" --whitelist_path "data/whitelist/V2.CellTag.Whitelist.csv"
./control.sh cellTag --visualize --save_progress_name "original" --tagged 2 --low_filter 2 --bam_data "data/tutorial/hf1.d15.possorted_genome_bam.bam"  --whitelist_version "v1" --sample_name "og"
./control.sh cellTag --visualize --save_progress_name "original" --tagged 2 --low_filter 2 --bam_data "data/tutorial/hf1.d15.possorted_genome_bam.bam"  --whitelist_version "v2" --sample_name "og"
./control.sh cellTag --visualize --save_progress_name "original" --tagged 2 --low_filter 2 --bam_data "data/tutorial/hf1.d15.possorted_genome_bam.bam"  --whitelist_version "v3" --sample_name "og"
./control.sh network --visualize --save_progress_name "original"