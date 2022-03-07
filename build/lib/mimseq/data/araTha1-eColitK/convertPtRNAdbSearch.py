#! /usr/bin/env python3

db = open("araTha1-PtRNAdb-searchResults.txt", "r")

with open("araTha1-plastidtRNAs.fa", "w") as out:
	for line in db:
		if not line.startswith("PtRNAdb"):
			line = line.strip()
			ID = line.split("\t")[0]
			species = line.split("\t")[3]
			isotype = line.split("\t")[9]
			isodecoder = line.split("\t")[10]
			seq = line.split("\t")[22]
			out.write(">" + ID + "|" + species + "|XXX|" + isotype + "|" + isodecoder + "\n")
			out.write(seq + "\n")