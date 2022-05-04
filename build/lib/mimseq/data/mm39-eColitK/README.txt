# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# Extract sequences
bedtools getfasta -fi /data/genome_misc/mouse_mm39/mm39.fa -fo mm39-tRNAs-all.tmp.fa -bed mm39-tRNAs.bed -name -s -fullHeader

# Adjust header naming convention for mimseq using hg38-tRNAs_name_map.txt
./FastaHeadersforMimseq.py mm39-tRNAs_name_map.txt mm39-tRNAs-detailed.out mm39-tRNAs-all.tmp.fa mm39-tRNAs-all.fa
rm mm39-tRNAs-all.tmp.fa

# Add E. coli Lys-TTT reference sequence to the top of mm39-tRNAs-all.fa
# cat out file for mm39 and E. coli and save
cat mm39-tRNAs-detailed.out ../eschColi-K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.out > mm39-tRNAs-detailed2.out
mv mm39-tRNAs-detailed2.out mm39-tRNAs-detailed.out