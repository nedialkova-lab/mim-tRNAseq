#! /usr/bin/env python3

###########################################################################
# Analysis of modifications/misincorporations and stops per tRNA position #
###########################################################################

import os, logging
import re
import pysam
from getCoverage import getBamList
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
from collections import defaultdict

log = logging.getLogger(__name__)

def countMods_mp(out_dir, cov_table, info, inputs):
# modification counting and table generation
	
	# generate DataFrame with the references (as index) and positions imported from cov_table
	modTable = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
	stopTable = defaultdict(lambda: defaultdict(int))
	geneCov = defaultdict(int)
	condition = info[inputs][0]
	bam_file = pysam.AlignmentFile(inputs, "rb")
	log.info('Analysing misincorporations and stops for {}...'.format(inputs))
	for read in bam_file.fetch(until_eof=True):
		query = read.query_name
		reference = read.reference_name
		geneCov[reference] += 1

		# get MD tags for mismatches, split into list of integers and characters
		md_tag = read.get_tag('MD')
		md_list = re.split('(.*?)([A-Za-z]|[\^][A-Za-z]+)', md_tag)
		md_list = list(filter(None, md_list))
		
		# get offset of read mapping position to reference start in order to adjust mismatch position 
		# (offset is simply start position of read alignment realtive to reference)
		offset = read.reference_start
		
		# offset + 1 (0 to 1 based) is start of alignment - i.e. any start > 1 indicates a stop to RT at this position
		stopTable[reference][offset+1] += 1
		
		# build mismatch dictionary ignoring inserts in reference ('^') due to alignment algorithm and adding offset only to first interval of MD tag
		ref_pos = 0
		read_pos = 0
		for index, interval in enumerate(md_list):
			if index == 0 and not interval.startswith('^'):
				if interval.isdigit(): # stretch of matches
					read_pos += int(interval)
					interval = int(interval) + offset
					ref_pos += interval
				elif interval.isalpha(): # is a mismatch
					identity = read.query_sequence[read_pos]
					ref_pos += 1
					read_pos += 1
					modTable[reference][ref_pos][identity] += 1 # log the identity of the misincorporated base
					if interval == interval.lower():
						modTable[reference][ref_pos]['known'] = True # log if the mismatch was a known modified position by checking for lowercase letter in MD tag (see --md-lowercase-snp in GSNAP parameters)
					else:
						modTable[reference][ref_pos]['known'] = False
			elif not interval.startswith('^'):
				if interval.isdigit(): # stretch of matches
					read_pos += int(interval)
					ref_pos += int(interval)
				elif interval.isalpha(): # is a mismatch
					identity = read.query_sequence[read_pos] # identity is misincorporated nucleotide
					read_pos += 1
					ref_pos += 1
					modTable[reference][ref_pos][identity] += 1
					if interval == interval.lower():
						modTable[reference][ref_pos]['known'] = True
					else:
						modTable[reference][ref_pos]['known'] = False
			elif interval.startswith('^'):
				insertion = len(interval) - 1
				ref_pos += insertion

	# build dictionaries for mismatches, known modification sites, and stops, normalizing to total coverage per nucleotide
	modTable_prop = {cluster: {pos: {
				group: count / cov_table[(cov_table.pos == pos) & (cov_table.condition == condition) & (cov_table.bam == inputs)].loc[cluster]['cov']
				for group, count in data.items() if group in ['A','C','G','T']
								}
			for pos, data in values.items()
							}
		for cluster, values in modTable.items()
					}

	knownTable = {cluster: {pos: test['known']
			for pos, test in values.items()
						}
		for cluster, values in modTable.items()
				}

	stopTable_prop = {cluster: {
			pos: count / geneCov[cluster]
			for pos, count in values.items()
						}
		for cluster, values in stopTable.items()
				}

	# reformat modTable and save to temp file
	reform = {(outerKey, innerKey): values for outerKey, innerDict in modTable_prop.items() for innerKey, values in innerDict.items()}
	modTable_prop_df = pd.DataFrame.from_dict(reform)
	modTable_prop_df['type'] = modTable_prop_df.index
	modTable_prop_melt = modTable_prop_df.melt(id_vars=['type'], var_name=['cluster','pos'], value_name='proportion')
	modTable_prop_melt['condition'] = condition
	modTable_prop_melt['bam'] = inputs
	modTable_prop_melt = modTable_prop_melt[['cluster','pos','type','proportion','condition', 'bam']]
	modTable_prop_melt.to_csv(inputs + "mismatchTable.csv", sep = "\t", index = False, na_rep = '0')

	# reformat knownTable and save to temp file
	knownTable_df = pd.DataFrame.from_dict(knownTable)
	knownTable_df['pos'] = knownTable_df.index
	knownTable_df_melt = knownTable_df.melt(id_vars='pos', var_name='cluster', value_name='known')
	knownTable_df_melt['condition'] = condition
	knownTable_df_melt['bam'] = inputs
	knownTable_df_melt = knownTable_df_melt[['cluster', 'pos', 'known', 'condition', 'bam']]
	knownTable_df_melt.to_csv(inputs + "knownModSites.csv", sep = "\t", index = False, na_rep = '0')

	# reformat stopTable and save to temp file
	stopTable_prop_df = pd.DataFrame.from_dict(stopTable_prop)
	stopTable_prop_df['pos'] = stopTable_prop_df.index
	stopTable_prop_melt = stopTable_prop_df.melt(id_vars='pos', var_name='cluster', value_name='proportion')
	stopTable_prop_melt['condition'] = condition
	stopTable_prop_melt['bam'] = inputs
	stopTable_prop_melt = stopTable_prop_melt[['cluster', 'pos', 'proportion', 'condition', 'bam']]
	stopTable_prop_melt.to_csv(inputs + "RTstopTable.csv", sep = "\t", index = False, na_rep = '0')

	log.info('Analysis complete for {}...'.format(inputs))

def generateModsTable(sampleGroups, out_dir, mod_lists, threads, cov_table):
# Wrapper function to call countMods_mp with multiprocessing

	log.info("\n+---------------------------------------------+\
	\n| Analysing misincorporations and stops to RT |\
	\n+---------------------------------------------+")

	# Get bam info using function in getCoverage
	baminfo, bamlist = getBamList(sampleGroups)

	if len(baminfo) > threads:
		multi = threads
	else:
		multi = len(baminfo)

	# initiate multiprocessing pool and run with bam names
	pool = Pool(multi)
	func = partial(countMods_mp, out_dir, cov_table, baminfo)
	pool.map(func, bamlist)
	pool.close()
	pool.join()

	modTable_total = pd.DataFrame()
	knownTable_total = pd.DataFrame()
	stopTable_total = pd.DataFrame()

	for bam in bamlist:
		# read in temp files and then delete
		modTable = pd.read_table(bam + "mismatchTable.csv", header = None)
		os.remove(bam + "mismatchTable.csv")
		knownTable = pd.read_table(bam + "knownModSites.csv", header = None)
		os.remove(bam + "knownModSites.csv")
		stopTable = pd.read_table(bam + "RTstopTable.csv", header = None)
		os.remove(bam + "RTstopTable.csv")

		# add individual temp files to main big table and save
		modTable_total = modTable_total.append(modTable)
		knownTable_total = knownTable_total.append(knownTable)
		stopTable_total = stopTable_total.append(stopTable)

	modTable_total.to_csv(out_dir + "mismatchTable.csv", sep = "\t", index = False, header = False)
	knownTable_total.to_csv(out_dir + "knownModsTable.csv", sep = "\t", index = False, header = False)
	stopTable_total.to_csv(out_dir + "RTstopTable.csv", sep = "\t", index = False, header = False)	





			