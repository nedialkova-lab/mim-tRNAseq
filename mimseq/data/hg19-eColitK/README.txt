# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# Extract sequences
bedtools getfasta -fi /data/genome_misc/hg19/hg19.fa -fo hg19-tRNAs-all.tmp.fa -bed hg19-tRNAs.bed -name -s -fullHeader

# Adjust header naming convention for mimseq using hg38-tRNAs_name_map.txt
./FastaHeadersforMimseq.py hg19-tRNAs_name_map.txt hg19-tRNAs-detailed.out hg19-tRNAs-all.tmp.fa hg19-tRNAs-all.fa
rm hg19-tRNAs-all.tmp.fa

# Add E. coli Lys-TTT reference sequence to the top of hg38-tRNAs-all.fa
