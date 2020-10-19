# Filter isotype mismatch, unexpected anticodon and undertermined isotype tRNAs from D. rerio reference

# generate list of genes to be filtered using "out" file
grep -E "isotype mismatch|unexpected anticodon|undetermined isotype" danRer11_eschColi-tRNAs.out | awk 'OFS=".trna" {print $1,$2}' > filterList.txt

# run filtertRNAs.py to output filtered tRNA fasta