# Filter isotype mismatch, unexpected anticodon and undertermined isotype tRNAs from D. rerio reference

# generate list of genes to be filtered using "out" file
grep -E "isotype mismatch|unexpected anticodon|undetermined isotype" danRer11_eschColi-tRNAs.out | awk 'OFS=".trna" {print $1,$2}' > filterList.txt

# match these chrx.trnaxxx names to fasta header and reverse match to get sequences to keep
grep -wf filterList.txt danRer11_eColitK.fa | grep -vf - <(grep ">" danRer11_eColitK.fa) > filteredMatches.txt

# Extract these (plus 2 lines for sequence) from the fasta file
grep --no-group-separator -A2 -f filteredMatches.txt danRer11_eColitK.fa > danRer11_eColitK_filtered.fa 