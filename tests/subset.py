#!/usr/bin/env python
import Bio.SeqIO
import gzip
import os

def create_subset(in_fastq_path, out_fastq_path, n_seqs = 1000):
    with gzip.open(in_fastq_path, 'rt') as file_handle:
        records = list()
        for seq_num, seq_record in enumerate(Bio.SeqIO.parse(file_handle, 'fastq')):
            if seq_num >= n_seqs:
                break
            records.append(seq_record)
    with gzip.open(out_fastq_path, 'wt') as file_handle:
        Bio.SeqIO.write(records, file_handle, 'fastq')

def main():
    for fastq_filename in "mimseq_hek_1.fastq.gz  mimseq_hek_2.fastq.gz  mimseq_k562_1.fastq.gz  mimseq_k562_2.fastq.gz".split():
        subset_filename = os.path.join('tests', 'data', fastq_filename.replace(".fastq.gz", ".subset.fastq.gz"))
        create_subset(fastq_filename, subset_filename)

if __name__ == "__main__":
    main()