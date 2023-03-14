# Arabidopsis TAIR10 reference from GtRNAdb
# Plastid and mito reference sequences from PtRNAdb - http://14.139.61.8/PtRNAdb/index.php
# convertPtRNAdbSearch.py to create araTha1-plastidtRNAs.fa and araTha1-mitotRNAs.fa in correct format as in other species mitotRNAs files

# Generate fasta file of full tRNA set using bed file and custom script to adjust tRNA naming in final fasta

# rename chromosome names in bed file. E.g., chr1 -> Chr1
sed -i -e 's/chr/Chr/g' araTha1-tRNAs.bed

# Extract sequences
bedtools getfasta -fi /data/genome_misc/arabidopsis_TAIR10/TAIR10_chr_all.fas -fo araTha1-tRNAs-all.tmp.fa -bed araTha1-tRNAs.bed -name -s -fullHeader

# Adjust header naming convention for mimseq using araTha1-tRNAs_name_map.txt
./FastaHeadersforMimseq.py araTha1-tRNAs_name_map.txt araTha1-tRNAs-detailed.out araTha1-tRNAs-all.tmp.fa araTha1-tRNAs-all.fa
rm araTha1-tRNAs-all.tmp.fa

# Add E. coli Lys-TTT reference sequence to the top of araTha1-tRNAs-all.fa

###### NB ######
Some tRNAs incorrectly named based on duplicates sequence

Arabidopsis_thaliana_tRNA-Arg-TCT-4-1 shares sequence with Arabidopsis_thaliana_tRNA-Arg-TCT-3-1. Renamed to Arabidopsis_thaliana_tRNA-Arg-TCT-3-2