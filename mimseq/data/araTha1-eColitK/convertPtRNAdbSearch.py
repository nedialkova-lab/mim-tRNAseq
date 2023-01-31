#! /usr/bin/env python3

files = ["araTha1-PtRNAdb-plastid-searchResults.txt", "araTha1-PtRNAdb-mito-searchResults.txt"]

for fn in files:
	org = fn.split("-")[2]
	db = open(fn, "r")

	with open("araTha1-" + org + "tRNAs.fa", "w") as out:
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