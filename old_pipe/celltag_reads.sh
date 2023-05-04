#bash

#filtering for celtags
#MEF
samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'GGT[ACTG]{8}GAATTC' > v1.celltag.reads.out

#D3
samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'GTGATG[ACTG]{8}GAATTC' > v2.celltag.reads.out

#D13
samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'TGTACG[ACTG]{8}GAATTC' > v3.celltag.reads.out

#parsing relevent information to .tsv
#MEF
./scripts/celltag.parse.reads.10x.sh -v tagregex="CCGGT([ACTG]{8})GAATTC" v1.celltag.reads.out > v1.celltag.parsed.tsv

#D3
./scripts/celltag.parse.reads.10x.sh -v tagregex="GTGATG([ACTG]{8})GAATTC" v2.celltag.reads.out > v2.celltag.parsed.tsv

#D13
./scripts/celltag.parse.reads.10x.sh -v tagregex="TGTACG([ACTG]{8})GAATTC" v3.celltag.reads.out > v3.celltag.parsed.tsv

#quantify the CellTag for each cell and create a Cell Barcode x CellTag matrix
#MEF
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v1.celltag.parsed.tsv hf1.d15.v1

#D3
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v2.celltag.parsed.tsv hf1.d15.v2

#D13
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v3.celltag.parsed.tsv hf1.d15.v3