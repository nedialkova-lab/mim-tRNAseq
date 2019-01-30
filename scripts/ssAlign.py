#! /usr/bin/env python3

## Mature tRNA cluster alignments and secondary structure parsing

import subprocess, os, re
from Bio import AlignIO
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

stkname = ''

def aligntRNA(tRNAseqs):
# run cmalign to generate Stockholm file for tRNA sequences
	global stkname
	print(tRNAseqs)
	stkname = tRNAseqs.split(".fa")[0] + '_align.stk'
	cmfile ='/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/data/tRNAmatureseq.cm'
	cmcommand = 'cmalign -o ' + stkname + ' --nonbanded -g ' + cmfile + ' ' + tRNAseqs
	subprocess.call(cmcommand, shell = True)
	print(stkname)

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
				term[pos] = term_type
			else:
				term_type = "3'"
				term[pos] = term_type

		# Acceptor stem
		if char == '(':
			acc[pos] = "Acceptor stem 5'"

		if char == ')':
			acc[pos] = "Acceptor stem 3'"

		# Internal stem loops
		if char == "<":
			open_count += 1
			stemloops[pos] = stemloops_type[stemloops_count]
		if char == ">":
			close_count +=1
			stemloops[pos] = stemloops_type[stemloops_count]
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
		if pos not in [x for x in struct_dict.keys()]:
			bulge_list.append(pos)

	# separate bulges into distinct lists based on consecutive positions
	for k, g in groupby(enumerate(bulge_list), lambda x:x[0] - x[1]):
		group = map(itemgetter(1), g)
		group = list(map(int,group))
		bulge_items.append(group)

	# add to struct_dict with naming
	for bulge in bulge_items:
		bulge_count += 1
		for pos in bulge:
			struct_dict[pos] = 'bulge' + str(bulge_count)

	return(struct_dict)

#def tRNAclassifier():

