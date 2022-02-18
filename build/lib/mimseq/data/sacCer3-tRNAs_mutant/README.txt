##--- Custom sacCer3 tRNA references with edited tRNA sequences for mutant strains ---##

##- General procedure is to copy the original reference, change the sequence where edited, and edit the name so that new reference forms a unique isodecoder (as long as that is the case), and is handled correctly with respect to the anticodon if that is the site of mutation.
##- Note also the "Mut" added directly after the isotype, e.g. ArgMut, so that this Mutated reference is distinguishable easily in the data.

# Note, in all cases the tRNAScanSE ID info for additional reference sequences is kept the same as the original reference that has been edited so that the same original entry in the "out" file is found and used for the new reference.

1. Leu-GAG to Leu-TAG
# Copy Saccharomyces_cerevisiae_tRNA-Leu-GAG-1-1
# Change name to Saccharomyces_cerevisiae_tRNA-LeuMut-TAG-2-1 (TAG-1 already exists)
# Chnage anticodon to TAG

2. Arg-CCT to Arg-TCT
# Copy Saccharomyces_cerevisiae_tRNA-Arg-CCT-1-1
# Change name to Saccharomyces_cerevisiae_tRNA-ArgMut-TCT-2-1 (TCT-1 already exists)
# Chnage anticodon to TCT