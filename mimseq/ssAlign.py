#! /usr/bin/env python3

## Mature tRNA cluster alignments and secondary structure parsing

import subprocess, os, re
from Bio import AlignIO
from collections import defaultdict, Counter
from itertools import groupby
from operator import itemgetter

stkname = ''

def aligntRNA(tRNAseqs, out):
# run cmalign to generate Stockholm file for tRNA sequences
	global stkname
	stkname = tRNAseqs.split(".fa")[0] + '_align.stk'
	cmfile = os.path.dirname(os.path.realpath(__file__)) + '/data/tRNAmatureseq.cm'
	cmcommand = ['cmalign', '-o', stkname, '--nonbanded', '-g', cmfile, tRNAseqs]
	subprocess.check_call(cmcommand, stdout = open(out + 'cm.log', 'w'))

def extraCCA():
	# look for extra CCA's added spuriously that fall outside of canonical tRNA structure
	# Seems to be a problem in certain sequences in mouse - either an artifact from gtRNAdb or tRNAScan, or CCA is genomically encoded for these tRNAs?
	extra_cca = list()
	stk = AlignIO.read(stkname, "stockholm")
	for record in stk:
		if record.seq[-3:] == 'cca': #lowercase here indicates alignment issue to other clusters
			extra_cca.append(record.name)

	os.remove(stkname)

	return(extra_cca)

def tRNAclassifier(ungapped = False):

	struct_dict = structureParser()
	stk = AlignIO.read(stkname, "stockholm")

	# Get canonical tRNA position numbering (cons_pos_list). Useful to retain cononical numbering of tRNA positions (i.e. anticodon at 34 - 36, m1A 58 etc...)
	# Return list of characters with pos or '-'. To be used in all plots with positional data such as heatmaps for stops or modifications.
	# cons_pos_dict is a dictionary of 1 based positions and canonical positions to create extra column in mods tables (mismatchTable and RTstopTable) mapping ungapped positions to canonical ones
	ss_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC SS_cons")])
	if ungapped:
		ss_cons_orig = ss_cons
		ss_cons = ss_cons.replace(".", "")
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

	# Loop thorugh every tRNA in alignment and create dictionary entry for pos-structure and pos-canonpos information (1-based to match to mismatchTable from mmQuant)
	tRNA_struct = defaultdict(dict)
	tRNA_ungap2canon = defaultdict(dict)

	for record in stk:
		tRNA = record.id
		seq = record.seq
		ungapped_pos = 0
		bases = ["A", "C", "G", "U"]

		if not ungapped:
			for i, letter in enumerate(seq, 1):
				if letter.upper() in bases:
					tRNA_ungap2canon[tRNA][ungapped_pos] = cons_pos_dict[i]
					ungapped_pos += 1
					tRNA_struct[tRNA][i] = struct_dict[i]
				else:
					tRNA_struct[tRNA][i] = 'gap'
		elif ungapped:
			n = 0
			for i, letter in enumerate(seq,1):
				if letter.upper() in bases:
					n += 1
					try:
						tRNA_ungap2canon[tRNA][ungapped_pos] = cons_pos_dict[n]
					except KeyError:
						break
					tRNA_struct[tRNA][n] = struct_dict[i]
					ungapped_pos += 1
				elif letter == "-" and ss_cons_orig[i-1] != ".":
					n += 1
					tRNA_struct[tRNA][n] = 'gap'
					#tRNA_ungap2canon[tRNA][ungapped_pos] = cons_pos_dict[n]
					#ungapped_pos += 1

	return(tRNA_struct, tRNA_ungap2canon, cons_pos_list, cons_pos_dict)

def tRNAclassifier_nogaps():

	struct_dict = structureParser()
	tRNA_struct = defaultdict(dict)

	# Loop thorugh every tRNA in alignment and create dictionary entry for pos - structure information (1-based to match to mismatchTable from mmQuant)
	stk = AlignIO.read(stkname, "stockholm")
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
	stk = AlignIO.read(stkname, "stockholm")
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

	cons_pos_list, cons_pos_dict = tRNAclassifier()[2:4]
	#anticodon = getAnticodon_1base()

	# Define positions of conserved mod sites in gapped alignment for each tRNA
	sites_dict = defaultdict()
	mod_sites = ['9', '20', '26', '32', '34', '37', '58']
	
	for mod in mod_sites:
		sites_dict[mod] = list(cons_pos_dict.keys())[list(cons_pos_dict.values()).index(mod)]

	upstream_dict = defaultdict(lambda: defaultdict(list))

	stk = AlignIO.read(stkname, "stockholm") 
	for record in stk:
		gene = record.id
		seq = record.seq
		for pos in sites_dict.values():
			identity = seq[pos-1] # identity of base at modification position
			if identity in ['A','C','G','U','T']:
				up = pos - 2 # pos is 1 based from struct, therefore -1 to make 0 based and -1 to get upstream nucl
				down = pos
				while seq[up].upper() not in ['A','C','G','U','T']:
					up -= 1
				while seq[down].upper() not in ['A','C','G','U','T']:
					down += 1
				upstream_dict[gene][pos].append(identity) 
				upstream_dict[gene][pos].append(seq[up]) # upstream base
				upstream_dict[gene][pos].append(seq[down]) # downstream base

	try:
		os.mkdir(out + "mods")
	except FileExistsError:
		pass

	with open(out + "mods/modContext.txt", 'w') as outfile:
		outfile.write("cluster\tpos\tidentity\tupstream\tdownstream\n")
		for cluster, data in upstream_dict.items():
			for pos, base in data.items():
				outfile.write(cluster + "\t" + str(pos) + "\t" + base[0] + "\t" + base[1] + "\t" + base[2] + "\n")

	mod_sites = str("_".join(str(e) for e in mod_sites))

	# mod_sites are canonical numberings of these positions
	# cons_pos_list is the full list of tRNA positions numbered according to canonical numbering scheme obtained from multiple seq alignments
	# cons_pos_dict is the same as list but a dictionary with matching gapped tRNA positions to each canonically numbered position
	return(mod_sites, cons_pos_list, cons_pos_dict)

def structureParser():
# read in stk file generated above and define structural regions for each tRNA input
	
	struct_dict = dict()
	# get conserved tRNA structure from alignment
	ss_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC SS_cons")])

	acc = defaultdict()

	term = defaultdict()
	term_type = "5'"

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
