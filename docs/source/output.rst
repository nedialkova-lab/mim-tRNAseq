Outputs
=======


mim-tRNAseq automatically generates many data and figure output files, some are dependent on paramter options specified when running mim-tRNAseq.
The output are split into various subdirectories described below:

**indices**

This folder contains the indices required for GSNAP alignment.

* `tRNAgenome\` contains the index of mature, processed tRNA transcripts
* `snp_index\` contains the SNP index generated from modified positions, needed by GSNAP for SNP-tolerant alignment. Only generated if --snp-tolerance is enabled.
* `\*.log` files contain log info from index generation

**annotation**

This folder contains various files describing the tRNA trasncriptome of the genome of interest.

* 