# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# Extract sequences
bedtools getfasta -fi /data/genome_misc/hg38/Homo_sapiens.GRCh38.primary_assembly.genome.fa -fo hg38-tRNAs-all.tmp.fa -bed hg38-tRNAs.bed -name -s -fullHeader

# Manually add two sequences from a reference scaffold not in the genome assembly (chr1_KI270713v1_random)

>tRNA-Glu-TTC-14-1::chr1_KI270713v1_random:31635-31706(-)
TCCCTGGTGGTCTAGTGGCtAGGATTCGGCGCTTTCACCGCTGCGGCCCGGGTTCGATTCCCGGTCAGGGAA
>tRNA-Asn-GTT-26-1::chr1_KI270713v1_random:28755-28828(-)
GTCTCTGTGGCGCAGTCGGTtAGCGCGTTCGGCTGTTAACCGGAAGGtTGGTGGTTCGAGACCACCCAAGGACG

# Adjust header naming convention for mimseq using hg38-tRNAs_name_map.txt
./FastaHeadersforMimseq.py hg38-tRNAs_name_map.txt hg38-tRNAs-detailed.out hg38-tRNAs-all.tmp.fa hg38-tRNAs-all.fa
rm hg38-tRNAs-all.tmp.fa

# Add E. coli Lys-TTT reference sequence to the top of hg38-tRNAs-all.fa

###### NB ######
Some tRNAs incorrectly named!
Homo_sapiens_tRNA-Asn-GTT-16-5 does not share exact sequence with other GTT-16 genes. Rename to Homo_sapiens_tRNA-Asn-GTT-28-1.
Homo_sapiens_tRNA-Asn-GTT-9-2 does not share sequence with GTT-9-1. Rename to Homo_sapiens_tRNA-Asn-GTT-14-1
Homo_sapiens_tRX-Asn-GTT-3-1 shares sequence with Homo_sapiens_tRX-Asn-GTT-2-1. Rename to Homo_sapiens_tRX-Asn-GTT-2-2
Homo_sapiens_tRX-Phe-GAA-1-1 shares sequence with Homo_sapiens_tRNA-Phe-GAA-10-1. Rename to Homo_sapiens_tRNA-Phe-GAA-10-2
Homo_sapiens_tRNA-Ala-AGC-7-1 shares sequence with Homo_sapiens_tRNA-Ala-AGC-2. Renamed to Homo_sapiens_tRNA-Ala-AGC-2-3


#### mitochondrial Tyr-GTA missing a 3'-A which leads to inaccuracies in alignment and shows as mostly CCA-less reads in CCA plots ####
#### Gene should end on CCA and gain an additional CCA after processing
#### manually added this A
#### see https://www.nature.com/articles/s41467-020-18068-6#Sec24 (supplementary data 4)