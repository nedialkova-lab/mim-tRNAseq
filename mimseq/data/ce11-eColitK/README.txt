# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# Extract sequences
bedtools getfasta -fi /data/genome_data/Caenorhabditis_elegans_ce11/ce11.fa -fo ce11-tRNAs-all.tmp.fa -bed ce11-tRNAs.bed -name -s -fullHeader

# Adjust header naming convention for mimseq using hg38-tRNAs_name_map.txt
./FastaHeadersforMimseq.py ce11-tRNAs_name_map.txt ce11-tRNAs-detailed.out ce11-tRNAs-all.tmp.fa ce11-tRNAs-all.fa

# Add E. coli Lys-TTT reference sequence to the top of ce11-tRNAs-all.fa

###### NB ######
Some tRNAs incorrectly named!


### Manually removed sequences due to strange sequence, poor score and issues with infernal alignment with mimseq
### See ce11-tRNAs-filtered.fa


#### mitochondrial Tyr-GTA missing a 3'-A which leads to inaccuracies in alignment and shows as mostly CCA-less reads in CCA plots ####
#### Gene should end on CCA and gain an additional CCA after processing
#### manually added this A
#### see https://www.nature.com/articles/s41467-020-18068-6#Sec24 (supplementary data 4)