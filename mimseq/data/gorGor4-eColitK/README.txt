# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# Edit bed file chromosome names to match genome
sed -e 's/chr//g' gorGor4-tRNAs.bed > gorGor4-tRNAs_noChr.bed
sed -i -e 's/^Un_//g' gorGor4-tRNAs_noChr.bed
sed -i -e 's/^1_//g' gorGor4-tRNAs_noChr.bed
sed -i -e 's/_random//g' gorGor4-tRNAs_noChr.bed
sed -i -e 's/v1/.1/g' gorGor4-tRNAs_noChr.bed

# Extract sequences
bedtools getfasta -fi /data/genome_misc/gorilla_gorGor4/Gorilla_gorilla.gorGor4.dna.toplevel.fa -fo gorGor4-tRNAs-all.tmp.fa -bed gorGor4-tRNAs_noChr.bed -name -s -fullHeader

# Adjust header naming convention for mimseq using hg38-tRNAs_name_map.txt
./FastaHeadersforMimseq.py gorGor4-tRNAs_name_map.txt gorGor4-tRNAs-detailed.out gorGor4-tRNAs-all.tmp.fa gorGor4-tRNAs-all.fa

# Add E. coli Lys-TTT reference sequence to the top of gorGor4-tRNAs-all.fa

##### Manually rename problem reference names #######

# >Gorilla_gorilla_tRNA-Gly-CCC-3-2 does not match sequence to Gly-CCC-3-1
# renamed to >Gorilla_gorilla_tRNA-Gly-CCC-6-1   

# Gorilla_gorilla_tRNA-Glu-CTC-2-1 has same sequence as Glu-CTC-1 isodecoder
# renamed to >Gorilla_gorilla_tRNA-Glu-CTC-1-6

# >Gorilla_gorilla_tRNA-Glu-CTC-1-1 does not match to other Glu-CTC-1 genes
# renamed to >Gorilla_gorilla_tRNA-Glu-CTC-2-1

# >Gorilla_gorilla_tRNA-Ser-AGA-1-2 does not match Ser-AGA-1-1
# renamed to >Gorilla_gorilla_tRNA-Ser-AGA-4-1

#### mitochondrial Tyr-GTA missing a 3'-A which leads to inaccuracies in alignment and shows as mostly CCA-less reads in CCA plots ####
#### Gene should end on CCA and gain an additional CCA after processing
#### manually added this A

#### NB! ####
####-----####

Gorilla_gorilla_tRNA-Pro-AGG-6-1 has been removed as NNN repeat causes problems with clustering
This is a high-confidence tRNA so removing it might not be the best fix, but the AGG anticodon family has very very few reads aligned anyway