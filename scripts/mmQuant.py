#! /usr/bin/env python3

###########################################################################
# Analysis of modifications/misincorporations and stops per tRNA position #
#    also includes counting of CCA vs CC ends required for CCA analysis   #
###########################################################################

import os, logging
import re
import pysam
import tRNAtools
from getCoverage import getBamList
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
from collections import defaultdict
import ssAlign

log = logging.getLogger(__name__)

def unknownMods(inputs, out_dir, knownTable, cluster_dict, modTable, misinc_thresh, cov_table, min_cov, tRNA_dict, cons_pos_dict):
# find unknown modifications with a total misincorporation threshold >= misinc_thresh

	log.info('Finding potential unannotated mods for {}'.format(inputs))
	new_mods =  defaultdict(list)
	new_Inosines = defaultdict(list)
	cons_anticodon = ssAlign.getAnticodon()
	
	for isodecoder, data in modTable.items():
		cluster = [parent for parent, child in cluster_dict.items() if isodecoder in child][0]
		#print(cluster)
		anticodon = ssAlign.clusterAnticodon(cons_anticodon, cluster)
		for pos, type in data.items():
			cov = cov_table[isodecoder][pos]
			if (sum(modTable[isodecoder][pos].values()) >= misinc_thresh and cov >= min_cov and pos-1 not in knownTable[cluster]): # misinc above threshold, cov above threshold and not previously known
				# if one nucleotide dominates misinc. pattern (i.e. >= 0.9 of all misinc, likely a true SNP or misalignment)
				if (max(modTable[isodecoder][pos].values()) / sum(modTable[isodecoder][pos].values()) >= 0.90):
					# if mod seems to be an inosine (i.e. A with G misinc at 34) add to list and modification SNPs file (see tRNAtools.ModsParser())
					if (tRNA_dict[isodecoder]['sequence'][pos-1] == 'A' and list(modTable[isodecoder][pos].keys())[list(modTable[isodecoder][pos].values()).index(max(modTable[isodecoder][pos].values()))] == 'G' and pos-1 == min(anticodon)):
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

def bamMods_mp(out_dir, min_cov, info, mismatch_dict, cluster_dict, cca, filtered_list, tRNA_struct, remap, misinc_thresh, knownTable, tRNA_dict, cons_pos_dict, unique_isodecoderMMs, splitBool, isodecoder_sizes, inputs):
# modification counting and table generation, and CCA analysis
	
	modTable = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
	stopTable = defaultdict(lambda: defaultdict(int))
	counts = defaultdict(lambda: defaultdict(int))
	cov = defaultdict(lambda: defaultdict(int))
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

		# get offset of read mapping position to reference start in order to adjust mismatch position, and end of alignment for coverage calculation
		# (offset is simply start position of read alignment realtive to reference)
		offset = read.reference_start
		aln_end = read.reference_end
		
		# map mismatches to reference and count mods
		ref_pos = 0
		read_pos = 0
		temp = defaultdict()
		old_reference = reference
		temp, ref_pos, read_pos, reference = countMods(temp, ref_pos, read_pos, read_seq, offset, reference, md_list, unique_isodecoderMMs, mismatch_dict)

		# read counts, stops and coverage
		counts[inputs][reference] += 1
		geneCov[reference] += 1
		# offset + 1 (0 to 1 based) is start of alignment - i.e. any start > 1 indicates a stop to RT at this position
		stopTable[reference][offset+1] += 1
		for i in range(offset+1, aln_end+1):
			cov[reference][i] += 1

		for pos, identity in temp.items():
			modTable[reference][pos][identity] += 1

		################
		# CCA analysis #
		################

		if cca:
			aln_count += 1
			dinuc = read.query_sequence[-2:]
			dinuc_dict[dinuc] += 1

			mapped_ref = read.reference_name

			ref_length = bam_file.get_reference_length(mapped_ref)
			if ref_pos in [ref_length, ref_length - 1]:
				cca_dict[reference][dinuc] += 1
			elif ref_pos == ref_length - 2:
				dinuc = read.query_sequence[-1:]
				cca_dict[reference][dinuc] += 1
			elif ref_pos == ref_length - 3:
				dinuc = "Absent"
				cca_dict[reference][dinuc] += 1

	## Edit misincorportation and stop data before writing

	# build dictionaries for mismatches and stops, normalizing to total coverage per nucleotide
	modTable_prop = {isodecoder: {pos: {
				group: count / cov[isodecoder][pos]
				  for group, count in data.items() if group in ['A','C','G','T']
								}
			for pos, data in values.items()
							}
		for isodecoder, values in modTable.items()
					}

	stopTable_prop = {isodecoder: {
			pos: count / geneCov[isodecoder]
			for pos, count in values.items()
						}
		for isodecoder, values in stopTable.items()
				}

	# if remapping is enabled, find uknown mod sites
	if remap:
		new_mods, new_Inosines = unknownMods(inputs, out_dir, knownTable, cluster_dict, modTable_prop, misinc_thresh, cov, min_cov, tRNA_dict, cons_pos_dict)
	else:
		new_mods = {}
		new_Inosines = {}

		log.info('Building modification, stop and count data tables for {}'.format(inputs))

		# reformat modTable, add gaps and structure, and save to temp file
		reform = {(outerKey, innerKey): values for outerKey, innerDict in modTable_prop.items() for innerKey, values in innerDict.items()}
		modTable_prop_df = pd.DataFrame.from_dict(reform)
		modTable_prop_df['type'] = modTable_prop_df.index
		modTable_prop_melt = modTable_prop_df.melt(id_vars=['type'], var_name=['isodecoder','pos'], value_name='proportion')
		modTable_prop_melt['condition'] = condition
		modTable_prop_melt['bam'] = inputs
		modTable_prop_melt.pos = pd.to_numeric(modTable_prop_melt.pos)
		# add coverage per nucelotide from cov
		cov_table = pd.DataFrame.from_dict(cov)
		cov_table['pos'] = cov_table.index
		cov_table['pos'] = cov_table['pos'].astype(int)
		cov_table_melt = cov_table.melt(id_vars='pos', var_name='isodecoder', value_name='cov')
		cov_table_melt['bam'] = inputs
		modTable_prop_melt = pd.merge(modTable_prop_melt, cov_table_melt, on = ['isodecoder', 'pos', 'bam'], how = 'left')

		modTable_prop_melt = addNA(modTable_prop_melt, tRNA_struct, cluster_dict, "mods")
		modTable_prop_melt = modTable_prop_melt[['isodecoder','pos', 'type','proportion','condition', 'bam', 'cov']]
		#modTable_prop_melt = modTable_prop_melt.join(tRNA_struct, on=['cluster', 'pos'])
		#modTable_prop_melt = modTable_prop_melt.dropna(subset=['struct'])

		modTable_prop_melt.to_csv(inputs + "mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')

		# reformat stopTable, add gaps and structure, and save to temp file
		stopTable_prop_df = pd.DataFrame.from_dict(stopTable_prop)
		stopTable_prop_df['pos'] = stopTable_prop_df.index
		stopTable_prop_melt = stopTable_prop_df.melt(id_vars='pos', var_name='isodecoder', value_name='proportion')
		stopTable_prop_melt['condition'] = condition
		stopTable_prop_melt['bam'] = inputs
		stopTable_prop_melt.pos = pd.to_numeric(stopTable_prop_melt.pos)

		stopTable_prop_melt = addNA(stopTable_prop_melt, tRNA_struct, cluster_dict, "stops")
		stopTable_prop_melt = stopTable_prop_melt[['isodecoder', 'pos', 'proportion', 'condition', 'bam']]
		#stopTable_prop_melt = stopTable_prop_melt.join(tRNA_struct, on = ['cluster', 'pos'])
		#stopTable_prop_melt = stopTable_prop_melt.dropna(subset=['struct'])

		stopTable_prop_melt.to_csv(inputs + "RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')

		# build temp counts DataFrame
		counts_table = pd.DataFrame.from_dict(counts)
		counts_table['isodecoder'] = counts_table.index
		counts_table = counts_table[['isodecoder', inputs]]
		# add 0 count isodecoders to table
		temp_add = pd.DataFrame(columns = ['isodecoder', inputs])
		for isodecoder in isodecoder_sizes.keys():
			if not isodecoder in counts_table['isodecoder']:
				temp_add = temp_add.append({'isodecoder':isodecoder, inputs:0}, ignore_index = True)
		counts_table = counts_table.append(temp_add)
		counts_table.to_csv(inputs + "countTable.csv", sep = "\t", index = False, na_rep = "NA")

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
					for dinuc, count in data.items():
						if (dinuc.upper() == "CC") or (dinuc.upper() == "CA") or (dinuc.upper() == "C") or (dinuc == "Absent"):
							CCAvsCC_counts.write(cluster + "\t" + dinuc + "\t" + inputs + "\t" + condition + "\t" + str(count) + "\n")

			dinuc_prop.close()
			CCAvsCC_counts.close()

	log.info('Analysis complete for {}...'.format(inputs))

	return(new_mods, new_Inosines)

def countMods(temp, ref_pos, read_pos, read_seq, offset, reference, md_list, unique_isodecoderMMs, mismatch_dict):
	
	for index, interval in enumerate(md_list):
		if not index == 0:
			new_offset = 0
		else:
			new_offset = offset
		if not interval.startswith('^'): 
			if interval.isdigit(): # stretch of matches
				read_pos += int(interval)
				interval = int(interval) + new_offset
				ref_pos += interval
			elif interval.isalpha(): # is a mismatch
				identity = read_seq[read_pos]
				ref_pos += new_offset 
				# check for position in mismatch dictionary from clustering
				# only include these positions if they aren't registered mismatches between clusters, or if they are known modified sites (lowercase)
				if identity in unique_isodecoderMMs[reference][ref_pos]:
					reference = unique_isodecoderMMs[reference][ref_pos][identity]
				elif (ref_pos not in mismatch_dict[reference]) or (ref_pos in mismatch_dict[reference] and interval.islower()):
					temp[ref_pos+1] = identity
				# move forward
				read_pos += 1
				ref_pos += 1
		elif interval.startswith('^'):
			insertion = len(interval) - 1
			ref_pos += insertion

	return(temp, ref_pos, read_pos, reference)

def addNA(table, tRNA_struct, cluster_dict, data_type):
	
	grouped = table.groupby('isodecoder')
	for name, group in grouped:
		cluster = [parent for parent, child in cluster_dict.items() if name in child][0]
		for pos in tRNA_struct.loc[cluster].index:
			if data_type == 'mods':
				new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1], 'cov':'NA'})
			elif data_type == 'stops':
				new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1]}, index=[0])
			if tRNA_struct.loc[cluster].iloc[pos-1].struct == 'gap':
				table.loc[(table.isodecoder == name) & (table.pos >= pos), 'pos'] += 1
				table = table.append(new)
			if not any(table.loc[table.isodecoder == name].pos == pos):
				table = table.append(new)

	return(table)

def generateModsTable(sampleGroups, out_dir, threads, min_cov, mismatch_dict, cluster_dict, filtered_list, cca, remap, misinc_thresh, knownTable, tRNA_dict, Inosine_clusters, unique_isodecoderMMs, splitBool, isodecoder_sizes):
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
	func = partial(bamMods_mp, out_dir, min_cov, baminfo, mismatch_dict, cluster_dict, cca, filtered_list, tRNA_struct_df, remap, misinc_thresh, knownTable, tRNA_dict, cons_pos_dict, unique_isodecoderMMs, splitBool, isodecoder_sizes)
	new_mods, new_Inosines = zip(*pool.map(func, bamlist))
	pool.close()
	pool.join()

	if not remap:

		modTable_total = pd.DataFrame()
		countsTable_total = pd.DataFrame()
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

			stopTable = pd.read_csv(bam + "RTstopTable.csv", header = 0, sep = "\t")
			stopTable['canon_pos'] = stopTable['pos'].map(cons_pos_dict)
			os.remove(bam + "RTstopTable.csv")

			countsTable = pd.read_csv(bam + "countTable.csv", header = 0, sep = "\t")
			os.remove(bam + "countTable.csv")

			# add individual temp files to main big table and save
			modTable_total = modTable_total.append(modTable)
			stopTable_total = stopTable_total.append(stopTable)
			if countsTable_total.empty:
				countsTable_total = countsTable_total.append(countsTable)
			else:
				countsTable_total = pd.merge(countsTable_total, countsTable, on = "isodecoder", how = "left")

			if cca:
				# same for CCA analysis files
				dinuc = pd.read_csv(bam + "_dinuc.csv", header = None, keep_default_na=False, sep = "\t")
				os.remove(bam + "_dinuc.csv")
				CCA = pd.read_table(bam + "_CCAcounts.csv", header = None)
				os.remove(bam + "_CCAcounts.csv")

				dinuc_table = dinuc_table.append(dinuc)
				CCAvsCC_table = CCAvsCC_table.append(CCA)

		# output tables
		modTable_total.to_csv(out_dir + "mods/mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')
		with open(out_dir + "mods/knownModsTable.csv", "w") as known:
			known.write("cluster\tpos\n")
			for cluster, data in knownTable.items():
				for pos in data:
					known.write(cluster + "\t" + str(pos) + "\n")

		stopTable_total.to_csv(out_dir + "mods/RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')	
		countsTable_total.to_csv(out_dir + "Isodecoder_counts.txt", sep = "\t", index = False)

		if cca:
			dinuc_table.columns = ['dinuc', 'proportion', 'sample']
			dinuc_table.to_csv(out_dir + "CCAanalysis/AlignedDinucProportions.csv", sep = "\t", index = False, na_rep = 'NA')
			CCAvsCC_table.columns = ['gene', 'end', 'sample', 'condition', 'count']
			CCAvsCC_table.to_csv(out_dir + "CCAanalysis/CCAcounts.csv", sep = "\t", index = False)

		# Anticodon counts
		tRNAtools.countReadsAnticodon(out_dir + "Isodecoder_counts.txt", out_dir)

		log.info("** Read counts per isodecoder saved to " + out_dir + "counts/Isodecoder_counts.txt **")

	return(new_mods, new_Inosines)

			