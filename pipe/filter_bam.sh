# filter unmapped reads
samtools view -b -f 4 ./hf1.d15.bam > ./hf1.d15.filtered.bam
# filter transgene reads
samtools view -b  ./hf1.d15.bam GFP >> ./hf1.d15.filtered.bam
# dont execute if CellTag UTR not included in ref
samtools view -b ./hf1.d15.bam CellTag.UTR >> ./hf1.d15.filtered.bam
