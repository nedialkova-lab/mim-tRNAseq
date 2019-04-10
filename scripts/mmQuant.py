#! /usr/bin/env python3

###########################################################################
# Analysis of modifications/misincorporations and stops per tRNA position #
#    also includes counting of CCA vs CC ends required for CCA analysis   #
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
import ssAlign

log = logging.getLogger(__name__)

def unknownMods(inputs, out_dir, knownTable, modTable, misinc_thresh):
# find unknown modifications with a total misincorporation threshold >= misinc_thresh

	log.info('Finding potential unannotated mods for {}'.format(inputs))
	new_mods =  defaultdict(list)
	for cluster, data in modTable.items():
		for pos, type in data.items():
			if (sum(modTable[cluster][pos].values()) >= misinc_thresh and pos-1 not in knownTable[cluster]):
				new_mods[cluster].append(pos-1) #modTable had 1 based values - convert back to 0 based for snp index

	with open(out_dir + "mods/predictedModstemp.csv", "w") as predMods:
		for cluster, data in new_mods.items():
			for pos in data:
				predMods.write(cluster + "\t" + \
					str(pos) + "\t" + \
					str(sum(modTable[cluster][pos+1].values())) + "\n")

	return(new_mods)

def countMods_mp(out_dir, cov_table, info, mismatch_dict, cca, filtered_list, tRNA_struct, remap, misinc_thresh, knownTable, inputs):
# modification counting and table generation, and CCA analysis
	
	modTable = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
	stopTable = defaultdict(lambda: defaultdict(int))
	geneCov = defaultdict(int)
	condition = info[inputs][0]

	if cca:
		aln_count = 0
		cca_dict = defaultdict(lambda: defaultdict(int))
		dinuc_dict = defaultdict(int)
		dinuc_prop = open(inputs + "_dinuc.csv", "w")
		CCAvsCC_counts = open(inputs + "_CCAcounts.csv", "w")

	bam_file = pysam.AlignmentFile(inputs, "rb")
	log.info('Analysing {}...'.format(inputs))
	for read in bam_file.fetch(until_eof=True):
		query = read.query_name
		reference = read.reference_name
		if reference in filtered_list:
			continue
		geneCov[reference] += 1

		#########################
		# Modification analysis #
		#########################

		# get MD tags for mismatches, split into list of integers and characters
		md_tag = read.get_tag('MD')
		md_list = re.split('(.*?)([A-Za-z]|[\^][A-Za-z]+)', md_tag)
		md_list = list(filter(None, md_list))

		# get cigar string, split as above
		cigar = read.cigarstring
		cigar_list = re.split('(.*?)([A-Za-z]|[\^][A-Za-z]+)', cigar)
		cigar_list = list(filter(None, cigar_list))

		# check cigar for softclipping at 5' of read (usually because of terminal transferase of RT) - if so remove soft-clipped bases from sequence
		if cigar_list[1].upper() == "S".upper():
			soft_clip = int(cigar_list[0])
			read_seq = read.query_sequence[soft_clip:]
		else:
			read_seq = read.query_sequence 

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
					identity = read_seq[read_pos]
					ref_pos += offset 
					# check for position in mismatch dictionary from clustering
					# only include these positions if they aren't registered mismatches between clusters, or if they are known modified sites (lowercase)
					if (ref_pos not in mismatch_dict[reference]) or (ref_pos in mismatch_dict[reference] and identity.islower()):
						modTable[reference][ref_pos+1][identity] += 1 # log the identity of the misincorporated base
						# if interval == interval.lower():
						# 	modTable[reference][ref_pos]['known'] = 1 # log if the mismatch was a known modified position by checking for lowercase letter in MD tag (see --md-lowercase-snp in GSNAP parameters)
						# else:
						# 	modTable[reference][ref_pos]['known'] = 0
					# move forward
					read_pos += 1
					ref_pos += 1

			elif not interval.startswith('^'):
				if interval.isdigit(): # stretch of matches
					read_pos += int(interval)
					ref_pos += int(interval)
				elif interval.isalpha(): # is a mismatch
					identity = read_seq[read_pos] # identity is misincorporated nucleotide
					# check for position in mismatch dictionary from clustering
					# only include these positions if they aren't registered mismatches between clusters, or if they are known modified sites (lowercase)
					if (ref_pos not in mismatch_dict[reference]) or (ref_pos in mismatch_dict[reference] and identity.islower()):
						modTable[reference][ref_pos+1][identity] += 1
						# if interval == interval.lower():
						# 	modTable[reference][ref_pos]['known'] = 1
						# else:
						# 	modTable[reference][ref_pos]['known'] = 0
					# move forward
					read_pos += 1
					ref_pos += 1
			elif interval.startswith('^'):
				insertion = len(interval) - 1
				ref_pos += insertion

		################
		# CCA analysis #
		################

		if cca:
			aln_count += 1
			dinuc = read.query_sequence[-2:]
			cca_dict[reference][dinuc] += 1
			dinuc_dict[dinuc] += 1

	if cca:
		# write dinuc proportions for current bam
		for dinuc, count in dinuc_dict.items():
			dinuc_prop.write(dinuc + "\t" + str(count/aln_count) + "\t" + inputs.split("/")[-1] + "\n")

		# write CCA outputs for current bam
		for cluster, data in cca_dict.items():
			for dinuc, count in data.items():
				if (dinuc.upper() == "CC") or (dinuc.upper() == "CA"):
					CCAvsCC_counts.write(cluster + "\t" + dinuc + "\t" + inputs + "\t" + condition + "\t" + str(count) + "\n")

		dinuc_prop.close()
		CCAvsCC_counts.close()

	## Edit misincorportation and stop data before writing

	# build dictionaries for mismatches, known modification sites, and stops, normalizing to total coverage per nucleotide

	modTable_prop = {cluster: {pos: {
				group: count / cov_table[(cov_table.pos == pos) & (cov_table.condition == condition) & (cov_table.bam == inputs)].loc[cluster]['cov']
				  for group, count in data.items() if group in ['A','C','G','T']
								}
			for pos, data in values.items()
							}
		for cluster, values in modTable.items()
					}

	# knownTable = {cluster: {pos: test['known']
	# 		for pos, test in values.items()
	# 					}
	# 	for cluster, values in modTable.items()
	# 			}

	stopTable_prop = {cluster: {
			pos: count / geneCov[cluster]
			for pos, count in values.items()
						}
		for cluster, values in stopTable.items()
				}

	# if remapping is enabled, find uknown mod sites
	if remap:
		new_mods = unknownMods(inputs, out_dir, knownTable, modTable_prop, misinc_thresh)
	else:
		new_mods = {}

	# reformat modTable, add gaps and structure, and save to temp file
	reform = {(outerKey, innerKey): values for outerKey, innerDict in modTable_prop.items() for innerKey, values in innerDict.items()}
	modTable_prop_df = pd.DataFrame.from_dict(reform)
	modTable_prop_df['type'] = modTable_prop_df.index
	modTable_prop_melt = modTable_prop_df.melt(id_vars=['type'], var_name=['cluster','pos'], value_name='proportion')
	modTable_prop_melt['condition'] = condition
	modTable_prop_melt['bam'] = inputs
	modTable_prop_melt.pos = pd.to_numeric(modTable_prop_melt.pos)
	
	grouped = modTable_prop_melt.groupby('cluster')
	for name, group in grouped:
		for pos in tRNA_struct.loc[name].index:
			if tRNA_struct.loc[name].iloc[pos-1].struct == 'gap':
				modTable_prop_melt.loc[(modTable_prop_melt.cluster == name) & (modTable_prop_melt.pos >= pos), 'pos'] += 1
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1]})
				modTable_prop_melt = modTable_prop_melt.append(new)
			if not any(modTable_prop_melt.loc[modTable_prop_melt.cluster == name].pos == pos):
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1]})
				modTable_prop_melt = modTable_prop_melt.append(new)


	modTable_prop_melt = modTable_prop_melt[['cluster','pos', 'type','proportion','condition', 'bam']]
	modTable_prop_melt = modTable_prop_melt.join(tRNA_struct, on=['cluster', 'pos'])

	modTable_prop_melt.to_csv(inputs + "mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')

	# reformat knownTable and save to temp file
	# knownTable_df = pd.DataFrame.from_dict(knownTable)
	# knownTable_df['pos'] = knownTable_df.index
	# knownTable_df_melt = knownTable_df.melt(id_vars='pos', var_name='cluster', value_name='known')
	# knownTable_df_melt['condition'] = condition
	# knownTable_df_melt['bam'] = inputs
	# knownTable_df_melt = knownTable_df_melt[['cluster', 'pos', 'known', 'condition', 'bam']]
	# knownTable_df_melt.to_csv(inputs + "knownModSites.csv", sep = "\t", index = False, na_rep = 'NA')

	# reformat stopTable, add gaps and structure, and save to temp file
	stopTable_prop_df = pd.DataFrame.from_dict(stopTable_prop)
	stopTable_prop_df['pos'] = stopTable_prop_df.index
	stopTable_prop_melt = stopTable_prop_df.melt(id_vars='pos', var_name='cluster', value_name='proportion')
	stopTable_prop_melt['condition'] = condition
	stopTable_prop_melt['bam'] = inputs
	stopTable_prop_melt.pos = pd.to_numeric(stopTable_prop_melt.pos)

	grouped = stopTable_prop_melt.groupby('cluster')
	for name, group in grouped:
		for pos in tRNA_struct.loc[name].index:
			if tRNA_struct.loc[name].iloc[pos-1].struct == 'gap':
				stopTable_prop_melt.loc[(stopTable_prop_melt.cluster == name) & (stopTable_prop_melt.pos >= pos), 'pos'] += 1
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1]}, index=[0])
				stopTable_prop_melt = stopTable_prop_melt.append(new)
			if not any(stopTable_prop_melt.loc[stopTable_prop_melt.cluster == name].pos == pos):
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1]}, index=[0])
				stopTable_prop_melt = stopTable_prop_melt.append(new)

	stopTable_prop_melt = stopTable_prop_melt[['cluster', 'pos', 'proportion', 'condition', 'bam']]
	stopTable_prop_melt = stopTable_prop_melt.join(tRNA_struct, on = ['cluster', 'pos'])

	stopTable_prop_melt.to_csv(inputs + "RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')

	log.info('Analysis complete for {}...'.format(inputs))

	return(new_mods)

def generateModsTable(sampleGroups, out_dir, threads, cov_table, mismatch_dict, filtered_list, cca, remap, misinc_thresh, knownTable):
# Wrapper function to call countMods_mp with multiprocessing

	if cca:
	
		log.info("\n+--------------------------------------------------------------------+\
		\n| Analysing misincorporations and stops to RT, and analysing 3' ends |\
		\n+--------------------------------------------------------------------+")

		try:
			os.mkdir(out_dir + "CCAanalysis")
			os.mkdir(out_dir + "mods")
		except FileExistsError:
			log.warning("Rewriting over old mods and CCA files...")

	else:

		log.info("\n+---------------------------------------------+\
		\n| Analysing misincorporations and stops to RT |\
		\n+---------------------------------------------+")

		try:
			os.mkdir(out_dir + "mods")
		except FileExistsError:
				log.warning("Rewriting over old mods files...")

	if remap:
		log.info("** Discovering unannotated modifications for realignment **")

	# Get bam info using function in getCoverage
	baminfo, bamlist = getBamList(sampleGroups)

	if len(baminfo) > threads:
		multi = threads
	else:
		multi = len(baminfo)

	# get tRNA struct info from ssAlign
	tRNA_struct = ssAlign.tRNAclassifier(out_dir)
	tRNA_struct_df = pd.DataFrame(tRNA_struct).unstack().rename_axis(('cluster', 'pos')).rename('struct')
	tRNA_struct_df = pd.DataFrame(tRNA_struct_df)

	# initiate multiprocessing pool and run with bam names
	pool = Pool(multi)
	func = partial(countMods_mp, out_dir, cov_table, baminfo, mismatch_dict, cca, filtered_list, tRNA_struct_df, remap, misinc_thresh, knownTable)
	new_mods = pool.map(func, bamlist)
	pool.close()
	pool.join()

	modTable_total = pd.DataFrame()
	#knownTable_total = pd.DataFrame()
	stopTable_total = pd.DataFrame()

	dinuc_table = pd.DataFrame()
	CCAvsCC_table = pd.DataFrame()

	for bam in bamlist:
		# read in temp files and then delete
		modTable = pd.read_table(bam + "mismatchTable.csv", header = 0)
		os.remove(bam + "mismatchTable.csv")
		#knownTable = pd.read_table(bam + "knownModSites.csv", header = 0)
		#os.remove(bam + "knownModSites.csv")
		stopTable = pd.read_table(bam + "RTstopTable.csv", header = 0)
		os.remove(bam + "RTstopTable.csv")

		# add individual temp files to main big table and save
		modTable_total = modTable_total.append(modTable)
		#knownTable_total = knownTable_total.append(knownTable)
		stopTable_total = stopTable_total.append(stopTable)

		if cca:
			# same for CCA analysis files
			dinuc = pd.read_table(bam + "_dinuc.csv", header = None)
			os.remove(bam + "_dinuc.csv")
			CCA = pd.read_table(bam + "_CCAcounts.csv", header = None)
			os.remove(bam + "_CCAcounts.csv")

			dinuc_table = dinuc_table.append(dinuc)
			CCAvsCC_table = CCAvsCC_table.append(CCA)

	modTable_total.to_csv(out_dir + "mods/mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')
	with open(out_dir + "mods/knownModsTable.csv", "w") as known:
		known.write("cluster\tpos\n")
		for cluster, data in knownTable.items():
			for pos in data:
				known.write(cluster + "\t" + str(pos) + "\n")
	#knownTable_total.to_csv(out_dir + "mods/knownModsTable.csv", sep = "\t", index = False, na_rep = 'NA')
	stopTable_total.to_csv(out_dir + "mods/RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')	

	if cca:
		dinuc_table.columns = ['dinuc', 'proportion', 'sample']
		dinuc_table.to_csv(out_dir + "CCAanalysis/AlignedDinucProportions.csv", sep = "\t", index = False)
		CCAvsCC_table.columns = ['gene', 'end', 'sample', 'condition', 'count']
		#CCAvsCC_table = CCAvsCC_table.pivot_table(index = 'gene', columns = 'sample', values = 'count', fill_value = 0)
		CCAvsCC_table.to_csv(out_dir + "CCAanalysis/CCAcounts.csv", sep = "\t", index = False)

	return(new_mods)

			