#! /usr/bin/env python3

## Mature tRNA cluster alignments and secondary structure parsing

import subprocess, os, re
from Bio import AlignIO
from Bio.Alphabet import generic_rna
from collections import defaultdict, Counter
from itertools import groupby
from operator import itemgetter

stkname = ''

def aligntRNA(tRNAseqs):
# run cmalign to generate Stockholm file for tRNA sequences
	global stkname
	stkname = tRNAseqs.split(".fa")[0] + '_align.stk'
	cmfile ='/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/data/tRNAmatureseq.cm'
	cmcommand = 'cmalign -o ' + stkname + ' --nonbanded -g ' + cmfile + ' ' + tRNAseqs
	subprocess.call(cmcommand, shell = True)

def tRNAclassifier(out):

	struct_dict = structureParser()
	tRNA_struct = defaultdict(dict)

	# Loop thorugh every tRNA in alignment and create dictionary entry for pos - structure information (1-based to match to mismatchTable from mmQuant)
	stk = AlignIO.read(stkname, "stockholm", alphabet=generic_rna)
	for record in stk:
		tRNA = record.id
		seq = record.seq
		bases = ["A", "C", "G", "U"]

		for i, letter in enumerate(seq, 1):
			if letter.upper() in bases:
				tRNA_struct[tRNA][i] = struct_dict[i]
			else:
				tRNA_struct[tRNA][i] = 'gap'

	# Get canonical tRNA position numbering (cons_pos_list). Useful to retain cononical numbering of tRNA positions (i.e. anticodon at 34 - 36, m1A 58 etc...)
	# Return list of characters with pos or '-'. To be used in all plots with positional data such as heatmaps for stops or modifications.
	# cons_pos_dict is a dictionary of 1 based positions and canonical positions to create extra column in mods tables (mismatchTable and RTstopTable) mapping ungapped positions to canonical ones
	ss_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC SS_cons")])
	cons_pos = 0
	cons_pos_list = list()
	cons_pos_dict = defaultdict()
	openstem_count = 0
	closestem_count = 0
	for pos, char in enumerate(ss_cons):
		if not ss_cons[pos] == ".":
			if cons_pos < 46:
				if (not cons_pos == 17) and (not cons_pos == 20):
					cons_pos_list.append(str(cons_pos))
					cons_pos_dict[pos+1] = str(cons_pos)
					cons_pos += 1
				elif (cons_pos == 17) and not ('17' in cons_pos_list):
					cons_pos_dict[pos+1] = '17'
					cons_pos_list.append('17')
				elif (cons_pos == 17) and ('17' in cons_pos_list):
					cons_pos_dict[pos+1] = '17a'
					cons_pos_list.append('17a')
					cons_pos += 1

				elif cons_pos == 20: 
					if not '20' in cons_pos_list:
						cons_pos_dict[pos+1] = '20'
						cons_pos_list.append('20')
					elif not '20a' in cons_pos_list:
						cons_pos_dict[pos+1] = '20a'
						cons_pos_list.append('20a')
					elif not '20b' in cons_pos_list:
						cons_pos_dict[pos+1] = '20b'
						cons_pos_list.append('20b')
						cons_pos += 1

			elif cons_pos == 46:
				if (not closestem_count == openstem_count) or (closestem_count == 0 or openstem_count == 0):
					if ss_cons[pos] == "<":
						openstem_count += 1
						cons_pos_dict[pos+1] = 'e'
						cons_pos_list.append("e")
					elif ss_cons[pos] == ">":
						closestem_count += 1
						cons_pos_dict[pos+1] = 'e'
						cons_pos_list.append("e")
					elif ss_cons[pos] == "_":
						cons_pos_dict[pos+1] = 'e'
						cons_pos_list.append("e")

				elif (closestem_count == openstem_count) and (not closestem_count == 0 or not openstem_count == 0):
					cons_pos_dict[pos+1] = str(cons_pos)
					cons_pos_list.append(str(cons_pos))
					cons_pos += 1

			elif cons_pos > 46:
				cons_pos_dict[pos+1] = str(cons_pos)
				cons_pos_list.append(str(cons_pos))
				cons_pos += 1

		elif ss_cons[pos] == ".":
			cons_pos_dict[pos+1] = '-'
			cons_pos_list.append("-")

	cons_pos_list = "_".join(cons_pos_list)

	return(tRNA_struct, cons_pos_list, cons_pos_dict)

def tRNAclassifier_nogaps():

	struct_dict = structureParser()
	tRNA_struct = defaultdict(dict)

	# Loop thorugh every tRNA in alignment and create dictionary entry for pos - structure information (1-based to match to mismatchTable from mmQuant)
	stk = AlignIO.read(stkname, "stockholm", alphabet=generic_rna)
	for record in stk:
		tRNA = record.id
		seq = record.seq
		pos = 0
		bases = ["A", "C", "G", "U"]

		for i, letter in enumerate(seq):
			if letter.upper() in bases:
				tRNA_struct[tRNA][pos] = struct_dict[i+1]
				pos += 1

	return(tRNA_struct)

def getAnticodon():
	# return anticodon position from conserved alignment

	anticodon = list()
	rf_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC RF")])
	# use '*' in rf_cons from stk to delimit the anticodon positions
	for pos, char in enumerate(rf_cons):
	 	if char == "*":
	 		anticodon.append(pos)

	return(anticodon)

def getAnticodon_1base():
	# return anticodon position from conserved alignment in 1-based positions - specifically for modContext()

	anticodon = list()
	rf_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC RF")])
	# use '*' in rf_cons from stk to delimit the anticodon positions
	for pos, char in enumerate(rf_cons, 1):
	 	if char == "*":
	 		anticodon.append(pos)

	return(anticodon)

def clusterAnticodon(cons_anticodon, cluster):
	# return anticodon position without gaps for specific cluster

	bases = ["A", "C", "G", "U"]
	stk = AlignIO.read(stkname, "stockholm", alphabet=generic_rna)
	cluster_anticodon = list()
	for record in stk:
		if record.id == cluster:
			for pos in cons_anticodon:
				gapcount = 0
				for char in record.seq[:pos]:
					if char.upper() not in bases:
						gapcount += 1
				cluster_anticodon.append(pos - gapcount)

	return(cluster_anticodon)

def modContext(out):
# outputs file of defined mods of interest pos, identity and context sequence for each cluster

	tRNA_struct, cons_pos_list, cons_pos_dict = tRNAclassifier(out)
	anticodon = getAnticodon_1base()

	# Define positions of conserved mod sites in gapped alignment for each tRNA
	sites_dict = defaultdict(dict)
	for tRNA in tRNA_struct:
		# 58
		section = [pos for gene, data in tRNA_struct.items() for pos, value in data.items() if value == "T stem-loop" and gene == tRNA]
		section.sort()
		if section and len(section) >= 8:
			sites_dict[tRNA]['58'] = section[-8]
		# 26
		section	= [pos for gene, data in tRNA_struct.items() for pos, value in data.items() if value == 'bulge2' and gene == tRNA]
		section.sort()
		if section:
			sites_dict[tRNA]['26'] = section[0]
		# 9
		section = [pos for gene, data in tRNA_struct.items() for pos, value in data.items() if value == 'bulge1' and gene == tRNA]
		section.sort()
		if section:
			sites_dict[tRNA]['9'] = section[-1]
		# 37
		anti_37 = max(anticodon) + 1
		while tRNA_struct[tRNA][anti_37] == 'gap':
			anti_37 += 1
		sites_dict[tRNA]['37'] = anti_37
		# 32
		anti_32 = min(anticodon) -2
		while tRNA_struct[tRNA][anti_32] == 'gap':
			anti_32 -= 1
		sites_dict[tRNA]['32'] = anti_32
		#47
		section	= [pos for gene, data in tRNA_struct.items() for pos, value in data.items() if value == 'Variable loop' and gene == tRNA]
		section.sort()
		if section and len(section) >= 2:
			sites_dict[tRNA]['47'] = section[-2]
		#34
		anti_34 = min(anticodon)
		sites_dict[tRNA]['34'] = anti_34

	mod_sites = ['9', '26', '32', '34', '37', '47', '58']

	upstream_dict = defaultdict(lambda: defaultdict(list))

	stk = AlignIO.read(stkname, "stockholm", alphabet=generic_rna) 
	for record in stk:
		gene = record.id
		seq = record.seq
		for tRNA, data in sites_dict.items():
			for site in data.keys():
				if gene == tRNA:
					pos = sites_dict[tRNA][site]
					up = pos - 2 # pos is 1 based from struct, therefore -1 to make 0 based and -1 to get upstream nucl
					down = pos
					while seq[up].upper() not in ['A','C','G','U','T']:
						up -= 1
					while seq[down].upper() not in ['A','C','G','U','T']:
						down += 1
					upstream_dict[tRNA][pos].append(seq[pos-1]) # identity of base at modification position
					upstream_dict[tRNA][pos].append(seq[up]) # upstream base
					upstream_dict[tRNA][pos].append(seq[down]) # downstream base


	with open(out + "mods/modContext.txt", 'w') as outfile:
		outfile.write("cluster\tpos\tidentity\tupstream\tdownstream\n")
		for cluster, data in upstream_dict.items():
			for pos, base in data.items():
				outfile.write(cluster + "\t" + str(pos) + "\t" + base[0] + "\t" + base[1] + "\t" + base[2] + "\n")

	# return 7 most abundant (gapped) positions in upstream_dict as consensus positions for 7 modifications retrieved above
	# pass these positions out to plotting script to facet plots by positions of modifications 

	counter = Counter()
	for cluster in upstream_dict:
		counter.update(upstream_dict[cluster].keys())

	cons_mod_pos = list()
	for pos, count in counter.most_common(7):
		cons_mod_pos.append(pos)

	cons_mod_pos = str("_".join(str(e) for e in cons_mod_pos))
	mod_sites = str("_".join(str(e) for e in mod_sites))

	# cons_mod_pos are ungapped positions of modification sites of interest, mod_sites are canonical numberings of these positions
	# cons_pos_list is the full list of tRNA positions numbered according to canonical numbering scheme obtained from multiple seq alignments
	return(cons_mod_pos, mod_sites, cons_pos_list, cons_pos_dict)

def structureParser():
# read in stk file generated above and define structural regions for each tRNA input
	
	struct_dict = dict()
	# get conserved tRNA structure from alignment
	ss_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC SS_cons")])

	acc = defaultdict()

	term = defaultdict()
	term_type = "5'"

	bulges = defaultdict()
	bulge_list = []
	bulge_items = []
	bulge_count = 0

	stemloops = defaultdict()
	stemloops_count = 0
	stemloops_type = ['D stem-loop','Anticodon stem-loop','Variable loop','T stem-loop']
	open_count = 0
	close_count = 0

	for pos, char in enumerate(ss_cons):
		# terminal ends
		if char == ':':
			if pos < 10:
				term[pos+1] = term_type
			else:
				term_type = "3'"
				term[pos+1] = term_type

		# Acceptor stem
		if char == '(':
			acc[pos+1] = "Acceptor stem 5'"

		if char == ')':
			acc[pos+1] = "Acceptor stem 3'"

		# Internal stem loops
		if char == "<":
			open_count += 1
			stemloops[pos+1] = stemloops_type[stemloops_count]
		if char == ">":
			close_count +=1
			stemloops[pos+1] = stemloops_type[stemloops_count]
			if close_count == open_count: # when the stems on either side have equal base pairs...
				stemloops_count += 1
				open_count = 0
				close_count = 0

	# create full ranges for acceptor stem
	for i in ["Acceptor stem 5'","Acceptor stem 3'"]:
		pos_list = [k for k, v in acc.items() if v == i]
		start = min(pos_list)
		stop = max(pos_list)
		acc.update([(n, i) for n in range(start, stop +1)])

	# create full ranges for stem loops
	for i in stemloops_type:
		pos_list = [k for k, v in stemloops.items() if v == i]
		start = min(pos_list)
		stop = max(pos_list)
		stemloops.update([(n, i) for n in range(start,stop+1)])

	# combine all into one dict
	struct_dict = {**term, **acc}
	struct_dict.update(stemloops)

	# bulges classification - i.e. everyhting that isn't already classified as a strucutral element
	for pos, char in enumerate(ss_cons):
		if (pos+1) not in [x for x in struct_dict.keys()]:
			bulge_list.append(pos+1)

	# separate bulges into distinct lists based on consecutive positions
	for k, g in groupby(enumerate(bulge_list), lambda x:x[0] - x[1]):
		group = map(itemgetter(1), g)
		group = list(map(int,group))
		bulge_items.append(group)

	# add bulges to struct_dict with naming
	for bulge in bulge_items:
		bulge_count += 1
		# bulges 3 and 4 form part of the variable loop
		if bulge_count == 3 or bulge_count == 4:
			for pos in bulge:
				struct_dict[pos] = 'Variable loop'
		# everything else is a bulge
		elif bulge_count == 5:
			for pos in bulge:
				struct_dict[pos] = 'bulge3' # set bulge 5 to bulge 3 because of variable loop bulges above
		else:
			for pos in bulge:
				struct_dict[pos] = 'bulge' + str(bulge_count)

	return(struct_dict)

