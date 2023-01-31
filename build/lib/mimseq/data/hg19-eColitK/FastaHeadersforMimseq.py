#! /usr/bin/env python3

import sys
from collections import defaultdict

name_map = sys.argv[1]
out_file = sys.argv[2]
in_fa = sys.argv[3]
out_fa = sys.argv[4] 

name_dict = defaultdict()

with open(name_map, "r") as names:
	for line in names:
		if not line.startswith("tRNAscan"):
			line = line.strip()
			name = line.split("\t")[1]
			name_dict[name] = line.split("\t")[0]

score_dict = defaultdict(list)

with open(out_file, "r") as scores:
	for line in scores:
		line = line.strip()
		if line.startswith("chr"):
			name = line.split("\t")[0].strip() + ".trna" + line.split("\t")[1]
			anticodon = line.split("\t")[5]
			score_dict[name].append(line.split("\t")[8])
			score_dict[name].append(anticodon)

with open(in_fa, "r") as fasta, open(out_fa, "w") as output:
	for line in fasta:
		line = line.strip()
		if line.startswith(">"):
			scan_id = name_dict[line.split(":")[0].replace(">","")]
			name = "Homo_sapiens_" + line.split(":")[0].replace(">","")
			if ("NNN" in name): #and not ("Und" or "Undet" in name):
				anticodon = score_dict[scan_id][1]
				name = "-".join(name.split("-")[0:2]) + "-" + anticodon + "-" + "-".join(name.split("-")[3:])
			else:
				anticodon = score_dict[scan_id][1]
			pos = line.split(":")[2:]
			start = int(pos[1].split("-")[0])
			stop = int(pos[1].split("-")[1].split("(")[0])
			length = str(stop - start)
			output.write(">" + name + " (tRNAscan-SE ID: " + scan_id + ") " + name.split("-")[1] + \
				" (" + anticodon + ") " + length + " bp Sc: " + score_dict[scan_id][0] + " " + ":".join(pos) + "\n")
		else:
			output.write(line + "\n")