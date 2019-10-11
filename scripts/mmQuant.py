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

def unknownMods(inputs, out_dir, knownTable, modTable, misinc_thresh, cov_table, min_cov, tRNA_dict, cons_pos_dict):
# find unknown modifications with a total misincorporation threshold >= misinc_thresh

	log.info('Finding potential unannotated mods for {}'.format(inputs))
	new_mods =  defaultdict(list)
	new_Inosines = defaultdict(list)
	cons_anticodon = ssAlign.getAnticodon()
	
	for cluster, data in modTable.items():
		anticodon = ssAlign.clusterAnticodon(cons_anticodon, cluster)
		for pos, type in data.items():
			cov = cov_table[(cov_table.pos == pos) & (cov_table.bam == inputs)].loc[cluster]['cov']
			if (sum(modTable[cluster][pos].values()) >= misinc_thresh and cov >= min_cov and pos-1 not in knownTable[cluster]): # misinc above threshold, cov above threshold and not previously known
				# if one nucleotide dominates misinc. pattern (i.e. >= 0.9 of all misinc, likely a true SNP or misalignment)
				if (max(modTable[cluster][pos].values()) / sum(modTable[cluster][pos].values()) >= 0.90):
					# if mod seems to be an inosine (i.e. A with G misinc at 34) add to list and modification SNPs file (see tRNAtools.ModsParser())
					if (tRNA_dict[cluster]['sequence'][pos-1] == 'A' and list(modTable[cluster][pos].keys())[list(modTable[cluster][pos].values()).index(max(modTable[cluster][pos].values()))] == 'G' and pos-1 == min(anticodon)):
						new_Inosines[cluster].append(pos-1)
				else:
					new_mods[cluster].append(pos-1) #modTable had 1 based values - convert back to 0 based for snp index

	with open(out_dir + "mods/predictedModstemp.csv", "w") as predMods:
		for cluster, data in new_mods.items():
			for pos in data:
				predMods.write(cluster + "\t" + \
					str(pos) + "\t" + \
					str(sum(modTable[cluster][pos+1].values())) + "\n")

	return(new_mods, new_Inosines)

def countMods_mp(out_dir, cov_table, min_cov, info, mismatch_dict, cca, filtered_list, tRNA_struct, remap, misinc_thresh, knownTable, tRNA_dict, cons_pos_dict, inputs):
# modification counting and table generation, and CCA analysis
	
	modTable = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
	clusterMMTable = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
	stopTable = defaultdict(lambda: defaultdict(int))
	geneCov = defaultdict(int)
	condition = info[inputs][0]

	# initialise structures and outputs if CCA analysis in on
	if cca:
		aln_count = 0
		cca_dict = defaultdict(lambda: defaultdict(int))
		dinuc_dict = defaultdict(int)
		dinuc_prop = open(inputs + "_dinuc.csv", "w")
		CCAvsCC_counts = open(inputs + "_CCAcounts.csv", "w")

	# process mods by looping through alignments in bam file
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

		# check cigar for softclipping at 5' and 3' of read (usually because of terminal transferase of RT) - if so remove soft-clipped bases from sequence
		if cigar_list[1].upper() == "S".upper():
			soft_clip = int(cigar_list[0])
			read_seq = read.query_sequence[soft_clip:]
		else:
			read_seq = read.query_sequence 

		if cigar_list[-1].upper() == "S".upper():
			soft_clip = int(cigar_list[-2])
			read_seq = read_seq[:-soft_clip]

		# get offset of read mapping position to reference start in order to adjust mismatch position 
		# (offset is simply start position of read alignment realtive to reference)
		offset = read.reference_start
		
		# offset + 1 (0 to 1 based) is start of alignment - i.e. any start > 1 indicates a stop to RT at this position
		stopTable[reference][offset+1] += 1
		
		# build mismatch dictionary ignoring inserts in reference ('^') due to alignment algorithm and adding offset only to first interval of MD tag
		# build clusterMMTable to count reads with known mismatches between cluster members to later split reads by isodecoders
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
					elif ref_pos in mismatch_dict[reference]:
						clusterMMTable[reference][ref_pos][identity] += 1 # here 0 based numbering used for pos as these positions will refer back to sequences positions and not used for referring to canonical numbering etc...
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
					elif ref_pos in mismatch_dict[reference]:
						clusterMMTable[reference][ref_pos][identity] += 1
					# move forward
					read_pos += 1
					ref_pos += 1
			elif interval.startswith('^'):
				insertion = len(interval) - 1
				identity = "insertion"
				clusterMMTable[reference][ref_pos][identity] += 1
				ref_pos += insertion

		################
		# CCA analysis #
		################

		if cca:
			aln_count += 1
			dinuc = read.query_sequence[-2:]
			dinuc_dict[dinuc] += 1

			ref_length = bam_file.get_reference_length(reference)
			if ref_pos in [ref_length, ref_length - 1]:
				cca_dict[reference][dinuc] += 1
			elif ref_pos == ref_length - 2:
				dinuc = read.query_sequence[-1:]
				cca_dict[reference][dinuc] += 1
			elif ref_pos == ref_length - 3:
				dinuc = "Absent"
				cca_dict[reference][dinuc] += 1

	if cca:
		# write dinuc proportions for current bam
		for dinuc, count in dinuc_dict.items():
			dinuc_prop.write(dinuc + "\t" + str(count/aln_count) + "\t" + inputs.split("/")[-1] + "\n")

		# add missing ends to dict to prevent issues with plotting in R
		for end in ["CA", "CC", "C", "Absent"]:
			for cluster, data in cca_dict.items():
				if not end in data.keys():
					cca_dict[cluster][end] = 0

		# write CCA outputs for current bam
		for cluster, data in cca_dict.items():
			for dinuc, count in data.items():
				if (dinuc.upper() == "CC") or (dinuc.upper() == "CA") or (dinuc.upper() == "C") or (dinuc == "Absent"):
					CCAvsCC_counts.write(cluster + "\t" + dinuc + "\t" + inputs + "\t" + condition + "\t" + str(count) + "\n")

		dinuc_prop.close()
		CCAvsCC_counts.close()

	## Edit misincorportation and stop data before writing

	# build dictionaries for mismatches, isodecoder counts, and stops, normalizing to total coverage per nucleotide
	modTable_prop = {cluster: {pos: {
				group: count / cov_table[(cov_table.pos == pos) & (cov_table.condition == condition) & (cov_table.bam == inputs)].loc[cluster]['cov']
				  for group, count in data.items() if group in ['A','C','G','T']
								}
			for pos, data in values.items()
							}
		for cluster, values in modTable.items()
					}

	clusterMMTable_prop = {cluster: {pos: {
				group: count / cov_table[(cov_table.pos == pos+1) & (cov_table.condition == condition) & (cov_table.bam == inputs)].loc[cluster]['cov']
				  for group, count in data.items() if group in ['A','C','G','T', 'insertion']
								}
			for pos, data in values.items()
							}
		for cluster, values in clusterMMTable.items()
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
		new_mods, new_Inosines = unknownMods(inputs, out_dir, knownTable, modTable_prop, misinc_thresh, cov_table, min_cov, tRNA_dict, cons_pos_dict)
	else:
		new_mods = {}
		new_Inosines = {}

	# reformat modTable, add gaps and structure, and save to temp file
	reform = {(outerKey, innerKey): values for outerKey, innerDict in modTable_prop.items() for innerKey, values in innerDict.items()}
	modTable_prop_df = pd.DataFrame.from_dict(reform)
	modTable_prop_df['type'] = modTable_prop_df.index
	modTable_prop_melt = modTable_prop_df.melt(id_vars=['type'], var_name=['cluster','pos'], value_name='proportion')
	modTable_prop_melt['condition'] = condition
	modTable_prop_melt['bam'] = inputs
	modTable_prop_melt.pos = pd.to_numeric(modTable_prop_melt.pos)
	# add coverage per nucelotide from cov_table
	cov_merge = cov_table[['pos', 'cov', 'bam']]
	cov_merge['cluster'] = cov_merge.index
	#cov_merge.columns = ['canon_pos', 'cov', 'bam', 'cluster']
	cov_merge['pos'] = cov_merge['pos'].astype(int)
	modTable_prop_melt = pd.merge(modTable_prop_melt, cov_merge, on = ['cluster', 'pos', 'bam'], how = 'left')
	
	grouped = modTable_prop_melt.groupby('cluster')
	for name, group in grouped:
		for pos in tRNA_struct.loc[name].index:
			if tRNA_struct.loc[name].iloc[pos-1].struct == 'gap':
				modTable_prop_melt.loc[(modTable_prop_melt.cluster == name) & (modTable_prop_melt.pos >= pos), 'pos'] += 1
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1], 'cov':'NA'})
				modTable_prop_melt = modTable_prop_melt.append(new)
			if not any(modTable_prop_melt.loc[modTable_prop_melt.cluster == name].pos == pos):
				new = pd.DataFrame({'cluster':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1], 'cov':'NA'})
				modTable_prop_melt = modTable_prop_melt.append(new)


	modTable_prop_melt = modTable_prop_melt[['cluster','pos', 'type','proportion','condition', 'bam', 'cov']]
	modTable_prop_melt = modTable_prop_melt.join(tRNA_struct, on=['cluster', 'pos'])
	modTable_prop_melt = modTable_prop_melt.dropna(subset=['struct'])

	modTable_prop_melt.to_csv(inputs + "mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')

	# reformat clusterMMTable and save to temp
	reform = {(outerKey, innerKey): values for outerKey, innerDict in clusterMMTable_prop.items() for innerKey, values in innerDict.items()}
	clusterMMTable_prop_df = pd.DataFrame.from_dict(reform)
	clusterMMTable_prop_df['type'] = clusterMMTable_prop_df.index
	clusterMMTable_prop_melt = clusterMMTable_prop_df.melt(id_vars=['type'], var_name=['cluster','pos'], value_name='proportion')
	clusterMMTable_prop_melt['bam'] = inputs
	clusterMMTable_prop_melt.pos = pd.to_numeric(clusterMMTable_prop_melt.pos)
	clusterMMTable_prop_melt = clusterMMTable_prop_melt.dropna()
	clusterMMTable_prop_melt.index = clusterMMTable_prop_melt['cluster']
	clusterMMTable_prop_melt = clusterMMTable_prop_melt.drop(columns = ['cluster'])
	clusterMMTable_prop_melt.to_csv(inputs + "clusterMMTable.csv", sep = "\t", index = True)

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
	stopTable_prop_melt = stopTable_prop_melt.dropna(subset=['struct'])

	stopTable_prop_melt.to_csv(inputs + "RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')

	log.info('Analysis complete for {}...'.format(inputs))

	return(new_mods, new_Inosines)

def generateModsTable(sampleGroups, out_dir, threads, cov_table, min_cov, mismatch_dict, filtered_list, cca, remap, misinc_thresh, knownTable, tRNA_dict, Inosine_clusters):
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
	tRNA_struct, cons_pos_list, cons_pos_dict = ssAlign.tRNAclassifier(out_dir)
	tRNA_struct_df = pd.DataFrame(tRNA_struct).unstack().rename_axis(('cluster', 'pos')).rename('struct')
	tRNA_struct_df = pd.DataFrame(tRNA_struct_df)

	# initiate multiprocessing pool and run with bam names
	pool = Pool(multi)
	func = partial(countMods_mp, out_dir, cov_table, min_cov, baminfo, mismatch_dict, cca, filtered_list, tRNA_struct_df, remap, misinc_thresh, knownTable, tRNA_dict, cons_pos_dict)
	new_mods, new_Inosines = zip(*pool.map(func, bamlist))
	pool.close()
	pool.join()

	modTable_total = pd.DataFrame()
	clusterMMTable_total = pd.DataFrame()
	#knownTable_total = pd.DataFrame()
	stopTable_total = pd.DataFrame()

	dinuc_table = pd.DataFrame()
	CCAvsCC_table = pd.DataFrame()

	for bam in bamlist:
		# read in temp files and then delete
		modTable = pd.read_csv(bam + "mismatchTable.csv", header = 0, sep = "\t")
		modTable['canon_pos'] = modTable['pos'].map(cons_pos_dict)
		for cluster in Inosine_clusters:
			modTable.at[(modTable.canon_pos == '34') & (modTable['type'] == 'G') & (modTable.cluster == cluster), 'proportion'] = 1 - sum(modTable[(modTable.canon_pos == '34') & (modTable['type'] != 'G') & (modTable.cluster == cluster)]['proportion'].dropna())
		os.remove(bam + "mismatchTable.csv")

		clusterMMTable = pd.read_csv(bam + "clusterMMTable.csv", header = 0, sep = "\t", index_col = 0)
		os.remove(bam + "clusterMMTable.csv")

		stopTable = pd.read_csv(bam + "RTstopTable.csv", header = 0, sep = "\t")
		stopTable['canon_pos'] = stopTable['pos'].map(cons_pos_dict)
		os.remove(bam + "RTstopTable.csv")

		# add individual temp files to main big table and save
		modTable_total = modTable_total.append(modTable)
		clusterMMTable_total = clusterMMTable_total.append(clusterMMTable)
		stopTable_total = stopTable_total.append(stopTable)

		if cca:
			# same for CCA analysis files
			dinuc = pd.read_csv(bam + "_dinuc.csv", header = None, keep_default_na=False, sep = "\t")
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
		dinuc_table.to_csv(out_dir + "CCAanalysis/AlignedDinucProportions.csv", sep = "\t", index = False, na_rep = 'NA')
		CCAvsCC_table.columns = ['gene', 'end', 'sample', 'condition', 'count']
		#CCAvsCC_table = CCAvsCC_table.pivot_table(index = 'gene', columns = 'sample', values = 'count', fill_value = 0)
		CCAvsCC_table.to_csv(out_dir + "CCAanalysis/CCAcounts.csv", sep = "\t", index = False)

	return(new_mods, new_Inosines, clusterMMTable_total)

			