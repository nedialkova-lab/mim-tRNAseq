##--- Custom Hsap38 tRNA references with edited tRNA sequences for Gly-TCC-4 edited B-box strain ---##

##- General procedure is to copy the original reference, change the sequence where edited, and edit the name so that new reference forms a unique isodecoder (as long as that is the case), and is handled correctly with respect to the anticodon if that is the site of mutation.
##- Note also the "Mut" added directly after the isotype, e.g. GlyMut, so that this Mutated reference is distinguishable easily in the data.

# Note, in all cases the tRNAScanSE ID info for additional reference sequences is kept the same as the original reference that has been edited so that the same original entry in the "out" file is found and used for the new reference.

1. Gly-TCC-4 to GlyMut-TCC-4
# Copy Homo_sapiens_tRNA-Gly-TCC-4-1
# Change name to Homo_sapiens_tRNA-GlyMut-TCC-4-1
# Change B-box sequence "...GGCTCGA..." to "...GGTTCGA..."