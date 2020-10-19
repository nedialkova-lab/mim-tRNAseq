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

### Manually removed sequences due to strange sequence, poor score and issues with infernal alignment with mimseq
### See hg38-tRNAs-filtered.fa
# Homo_sapiens_tRX-Phe-GAA-2-1 (tRNAscan-SE ID: chr1.trna134) Phe (GAA) 80 bp Sc: 26.6 chr1:121009543-121009623(-)
# Homo_sapiens_tRX-Lys-CTT-5-1 (tRNAscan-SE ID: chr7.trna1) Lys (CTT) 84 bp Sc: 31.1 chr7:12601963-12602047(+)
# Homo_sapiens_tRX-Lys-CTT-4-1 (tRNAscan-SE ID: chr7.trna31) Lys (CTT) 73 bp Sc: 34.5 chr7:97141572-97141645(-)
# Homo_sapiens_tRNA-Phe-GAA-9-1 (tRNAscan-SE ID: chr1.trna124) Phe (GAA) 80 bp Sc: 32.2 chr1:143792723-143792803(-)
# Homo_sapiens_tRX-Arg-CCT-2-1 (tRNAscan-SE ID: chr11.trna9) Arg (CCT) 79 bp Sc: 21.2 chr11:118241371-118241450(+)
# Homo_sapiens_tRX-Ile-AAT-4-1 (tRNAscan-SE ID: chr17.trna38) Ile (AAT) 87 bp Sc: 20.5 chr17:8206072-8206159(-)
# Homo_sapiens_tRNA-Glu-TTC-13-1 (tRNAscan-SE ID: chr2.trna6) Glu (TCA) 73 bp Sc: 21.3 chr2:74896918-74896991(+)
# Homo_sapiens_tRNA-Tyr-GTA-11-1 (tRNAscan-SE ID: chr7.trna29) Tyr (GTA) 77 bp Sc: 37.1 chr7:149356652-149356729(-)
# Homo_sapiens_tRNA-Asn-GTT-20-1 (tRNAscan-SE ID: chr1.trna116) Asn (GTT) 76 bp Sc: 41.0 chr1:146049197-146049273(-)
# Homo_sapiens_tRX-Cys-GCA-4-1 (tRNAscan-SE ID: chr7.trna16) Cys (GCA) 83 bp Sc: 26.6 chr7:149628152-149628235(+)
# Homo_sapiens_tRNA-Cys-GCA-25-1 (tRNAscan-SE ID: chr3.trna11) Cys (GCA) 80 bp Sc: 40.9 chr3:17699892-17699972(-)