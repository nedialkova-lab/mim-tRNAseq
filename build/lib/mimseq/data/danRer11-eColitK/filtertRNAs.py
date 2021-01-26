#! /usr/bin/env python3

from Bio import SeqIO
import re

filterList = list()
with open("filterList.txt", "r") as filterFile:
	for line in filterFile:
		line = line.strip()
		filterList.append(line)

seq_dict = SeqIO.to_dict(SeqIO.parse("danRer11_eColitK.fa","fasta"))

with open("danRer11_eColitK_filtered.fa", "w") as filteredSeqs:
	for seq in seq_dict:
		id = re.search("tRNAscan-SE ID: (.*?)\).|\((chr.*?)-",seq_dict[seq].description).groups()
		id = list(filter(None, id))[0]
		if not id in filterList:
			filteredSeqs.write(">" + seq_dict[seq].description + "\n" + str(seq_dict[seq].seq) + "\n")

