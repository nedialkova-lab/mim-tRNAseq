#! /usr/bin/env python3

###########################################################################
# Analysis of modifications/misincorporations and stops per tRNA position #
#    also includes counting of CCA vs CC ends required for CCA analysis   #
###########################################################################

from __future__ import absolute_import
import os, logging
import re
import pysam
from .tRNAtools import countReads
from .getCoverage import getBamList, filterCoverage
from multiprocessing import Pool
import multiprocessing.pool
from functools import partial
import pandas as pd
import numpy as np
from collections import defaultdict
from .ssAlign import getAnticodon, clusterAnticodon, tRNAclassifier

log = logging.getLogger(__name__)

# custom classes to allow non-demonic processes allowing children processes to spawn more children sub-processes (i.e. multiprocessing within multiprocessing)
class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def unknownMods(inputs, out_dir, knownTable, cluster_dict, modTable, misinc_thresh, cov_table, min_cov, tRNA_dict, remap):
# find unknown modifications with a total misincorporation threshold >= misinc_thresh

	log.info('Finding potential unannotated mods for {}'.format(inputs))
	new_mods_cluster =  defaultdict(list)
	new_inosines_cluster = defaultdict(list)
	new_mods_isodecoder = defaultdict(list)
	new_inosines_isodecoder = defaultdict(list)
	cons_anticodon = getAnticodon()
	
	for isodecoder, data in modTable.items():
		if cluster_dict:
			cluster = [parent for parent, child in cluster_dict.items() if isodecoder in child][0]
		# if cluster_dict is empty, then clustering is disabled and in this case isodecoder in modTable is the "cluster"
		else:
			cluster = isodecoder
		anticodon = clusterAnticodon(cons_anticodon, cluster)
		for pos, type in data.items():
			cov = cov_table[isodecoder][pos]
			if (sum(modTable[isodecoder][pos].values()) >= misinc_thresh and cov >= min_cov and pos-1 not in knownTable[cluster]): # misinc above threshold, cov above threshold and not previously known
				# if one nucleotide dominates misinc. pattern (i.e. >= 0.9 of all misinc, likely a true SNP or misalignment)
				if (max(modTable[isodecoder][pos].values()) / sum(modTable[isodecoder][pos].values()) >= 0.90):
					# if mod seems to be an inosine (i.e. A with G misinc at 34) add to list and modification SNPs file (see tRNAtools.ModsParser())
					if (tRNA_dict[isodecoder]['sequence'][pos-1] == 'A' and list(modTable[isodecoder][pos].keys())[list(modTable[isodecoder][pos].values()).index(max(modTable[isodecoder][pos].values()))] == 'G' and pos-1 == min(anticodon)):
						new_inosines_cluster[cluster].append(pos-1)
						new_inosines_isodecoder[isodecoder].append(pos-1)
				else:
					new_mods_cluster[cluster].append(pos-1) #modTable had 1 based values - convert back to 0 based for snp index
					new_mods_isodecoder[isodecoder].append(pos-1)

	with open(inputs + "_predictedModstemp.csv", "a") as predMods:
		#predMods.write("isodecoder\tpos\tidentity\tbam\n")
		for isodecoder, data in new_mods_isodecoder.items():
			for pos in data:
				predMods.write(isodecoder + "\t" + \
					str(pos) + "\t" + \
					str(tRNA_dict[isodecoder]['sequence'][int(pos)]) + "\t" + \
					inputs.split("/")[-1] + "\n")
		if not remap: # only write inosines on second round (or if remap is disabled) to prevent writing entries twice
			for isodecoder, data in new_inosines_isodecoder.items():
				for pos in data:
					predMods.write(isodecoder + "\t" + \
						str(pos) + "\t" + \
						"I" + "\t" + \
						inputs.split("/")[-1] + "\n")

	return(new_mods_cluster, new_inosines_cluster)

def bamMods_mp(out_dir, min_cov, info, mismatch_dict, insert_dict, cluster_dict, cca, tRNA_struct, remap, misinc_thresh, knownTable, tRNA_dict, unique_isodecoderMMs, splitBool, isodecoder_sizes, threads, inputs):
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
			del cigar_list[0:2]
		else:
			read_seq = read.query_sequence 

		if cigar_list[-1].upper() == "S".upper():
			soft_clip = int(cigar_list[-2])
			read_seq = read_seq[:-soft_clip]
			del cigar_list[-3:-1]

		# remove insertions in read from the sequence (MD tags do not account for deletions and effect misinc. identity matching)
		read_ins_pos = 0
		for index, i in enumerate(cigar_list):
			if i.isdigit():
				read_ins_pos += int(i)
			elif i.isalpha() and i.upper() == 'I':
				ins_size = int(cigar_list[index-1])
				read_seq = read_seq[0:read_ins_pos-ins_size] + read_seq[read_ins_pos:]
				read_ins_pos -= ins_size

		# get offset of read mapping position to reference start in order to adjust mismatch position, and end of alignment for coverage calculation
		# (offset is simply start position of read alignment realtive to reference)
		offset = read.reference_start
		aln_end = read.reference_end
		
		# map mismatches to reference and count mods
		ref_pos = 0
		read_pos = 0
		temp = defaultdict()
		temp, ref_pos, read_pos, reference = countMods(temp, ref_pos, read_pos, read_seq, offset, reference, md_list, unique_isodecoderMMs, mismatch_dict, insert_dict, remap)

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
	# readthroughTable is similar to stops but normalised by coverage at each base
	# This reflects the proportion of reads at a given site that stop at this site, as opposed to the proportion of all reads for the reference that stop here

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

	readthroughTable = {isodecoder: {
			pos: 1 - (count / cov[isodecoder][pos]) # 1 - stops gives readthrough
			for pos, count in values.items()
						}
		for isodecoder, values in stopTable.items()
				}
				
	# find unknown mod sites
	new_mods, new_Inosines = unknownMods(inputs, out_dir, knownTable, cluster_dict, modTable_prop, misinc_thresh, cov, min_cov, tRNA_dict, remap)
	
	# format and output mods and stops if remap is disabled (i.e. also occurs after round 2 of mapping)
	if not remap:
		new_mods = {}
		new_Inosines = {}

		log.info('Building modification, stop and count data tables for {}'.format(inputs))

		# reformat modTable and add gaps
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
		cov_table_melt.dropna(inplace = True)
		cov_table_melt['bam'] = inputs
		cov_table_melt = cov_table_melt[['isodecoder', 'pos', 'bam', 'cov']]
		cov_table_melt.to_csv(out_dir + inputs.split("/")[-1] + "_coverage.txt", sep = "\t", index = False)

		# split and parallelize addNA
		names, dfs = splitTable(modTable_prop_melt)
		pool = Pool(threads)
		func = partial(addNA, tRNA_struct, cluster_dict, "mods")
		modTable_prop_melt = pd.concat(pool.starmap(func, zip(names, dfs)))
		pool.close()
		pool.join()
		#modTable_prop_melt = addNA(modTable_prop_melt, tRNA_struct, cluster_dict, "mods")
		modTable_prop_melt = pd.merge(modTable_prop_melt, cov_table_melt, on = ['isodecoder', 'pos', 'bam'], how = 'left')
		modTable_prop_melt = modTable_prop_melt[['isodecoder','pos', 'type','proportion','condition', 'bam', 'cov']]

		# reformat stopTable and add gaps 
		stopTable_prop_df = pd.DataFrame.from_dict(stopTable_prop)
		stopTable_prop_df['pos'] = stopTable_prop_df.index
		stopTable_prop_melt = stopTable_prop_df.melt(id_vars='pos', var_name='isodecoder', value_name='proportion')
		stopTable_prop_melt['condition'] = condition
		stopTable_prop_melt['bam'] = inputs
		stopTable_prop_melt.pos = pd.to_numeric(stopTable_prop_melt.pos)

		# split and parallelize addNA
		stopTable_prop_melt.dropna(inplace = True)
		names, dfs = splitTable(stopTable_prop_melt)
		pool = Pool(threads)
		func = partial(addNA, tRNA_struct, cluster_dict, "stops")
		stopTable_prop_melt = pd.concat(pool.starmap(func, zip(names, dfs)))
		pool.close()
		pool.join()
		#stopTable_prop_melt = addNA(stopTable_prop_melt, tRNA_struct, cluster_dict, "stops")
		stopTable_prop_melt = stopTable_prop_melt[['isodecoder', 'pos', 'proportion', 'condition', 'bam']]

		# reformat readthroughTable
		readthroughTable = pd.DataFrame.from_dict(readthroughTable)
		readthroughTable['pos'] = readthroughTable.index
		readthroughTable_melt = readthroughTable.melt(id_vars='pos', var_name='isodecoder', value_name='proportion')
		readthroughTable_melt['condition'] = condition
		readthroughTable_melt['bam'] = inputs
		readthroughTable_melt.pos = pd.to_numeric(readthroughTable_melt.pos)

		# split and parallelize addNA
		readthroughTable_melt.dropna(inplace = True)
		names, dfs = splitTable(readthroughTable_melt)
		pool = Pool(threads)
		func = partial(addNA, tRNA_struct, cluster_dict, "stops") #"readthrough")
		readthroughTable_melt = pd.concat(pool.starmap(func, zip(names, dfs)))
		pool.close()
		pool.join()

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

		# save tables to temp files per sample
		counts_table.to_csv(inputs + "countTable.csv", sep = "\t", index = False, na_rep = "0")
		modTable_prop_melt.to_csv(inputs + "mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')
		stopTable_prop_melt.to_csv(inputs + "RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')
		readthroughTable_melt.to_csv(inputs + "readthroughTable.csv", sep = "\t", index = False, na_rep = 'NA')

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

	log.info('Analysis complete for {}...'.format(inputs))

	return(new_mods, new_Inosines)

def splitTable(table):
# splits various tables so that addNA can be run in parallel on groupby objects

	split_table = table.groupby('isodecoder')
	dfs = [group for name, group in split_table]
	names = [name for name, group in split_table]

	return(names, dfs)

def countMods(temp, ref_pos, read_pos, read_seq, offset, reference, md_list, unique_isodecoderMMs, mismatch_dict, insert_dict, remap):
# Loop though mismatches in read and return new reference for splitting reads and/or count mods
	
	old_reference = reference
	insertions = list()
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
				# update reference for isodecoder splitting if cluster_id not 1 and remap is disabed or this is round 2 of alignment (avoid errors in adding new mods for clusters)
				if (unique_isodecoderMMs) and (identity.upper() in unique_isodecoderMMs[old_reference][ref_pos]) and (not remap):
					reference = unique_isodecoderMMs[old_reference][ref_pos][identity]
				# only include these positions if they aren't registered mismatches between clusters, or if they are known modified sites (lowercase)
				elif (ref_pos not in mismatch_dict[old_reference]): #or (ref_pos in mismatch_dict[old_reference] and interval.islower()):
					temp[ref_pos+1] = identity
				# move forward
				read_pos += 1
				ref_pos += 1
		elif interval.startswith('^'):
			identity = 'insertion'
			#insertions.extend([x + ref_pos for x in range(len(interval) - 1)]) # register all insertions (including consecutive insertions) to be checked later against cluster parent
			insertions.append(ref_pos)
			if (unique_isodecoderMMs) and (identity in unique_isodecoderMMs[old_reference][ref_pos]) and (not remap):
				reference = unique_isodecoderMMs[old_reference][ref_pos][identity]
			insertion = len(interval) - 1
			ref_pos += insertion

	# handle insertions between cluster parent and members present in read (different from insertions in read only)
	# i.e. once new ref is found above and this member has an insertion relative to parent, subtract 1 from all misinc positions after the insertion
	# corrects for difference in length between member and parent. 
	# Note that 1 bp up and down from the insertion are checked to account for differences in insertion due to short read aligner, and usearch clustering
	if (not reference == old_reference) and (insertions):
		for i in insertions:
			if reference in (insert_dict[old_reference][i] or insert_dict[old_reference][i+1] or insert_dict[old_reference][i-1]):
				temp = {(k - 1 if k > i else k):v for k,v in temp.items()}

	return(temp, ref_pos, read_pos, reference)

def addNA(tRNA_struct, cluster_dict, data_type, name, table):
# fill mods and stops tables with 'NA' for gapped alignment

	#cluster = [parent for parent, child in cluster_dict.items() if name in child][0]
	for pos in tRNA_struct.loc[name].index:
		if data_type == 'mods':
			new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':table.condition.iloc[1], 'bam':table.bam.iloc[1]})
		elif data_type == 'stops':
			new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'proportion':'NA', 'condition':table.condition.iloc[0], 'bam':table.bam.iloc[0]}, index=[0])
		if tRNA_struct.loc[name].iloc[pos-1].struct == 'gap':
			if not pos == max(tRNA_struct.loc[name].index):
				table.loc[(table.isodecoder == name) & (table.pos >= pos), 'pos'] += 1
			#if not data_type == "readthrough":
			table = table.append(new)
		if not any(table.loc[table.isodecoder == name].pos == pos): #and not (data_type == "readthrough"):
			table = table.append(new)

	return(table)

# def addNAOLD(tRNA_struct, cluster_dict, data_type_dict, table):
# # fill mods and stops tables with 'NA' for gapped alignment

# 	grouped = table.groupby('isodecoder')
# 	for name, group in grouped:
# 		#cluster = [parent for parent, child in cluster_dict.items() if name in child][0]
# 		for pos in tRNA_struct.loc[name].index:
# 			if data_type == 'mods':
# 				new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'type':pd.Categorical(['A','C','G','T']), 'proportion':'NA', 'condition':group.condition.iloc[1], 'bam':group.bam.iloc[1], 'cov':'NA'})
# 			elif data_type == 'stops':
# 				new = pd.DataFrame({'isodecoder':name, 'pos':pos, 'proportion':'NA', 'condition':group.condition.iloc[0], 'bam':group.bam.iloc[0]}, index=[0])
# 			if tRNA_struct.loc[name].iloc[pos-1].struct == 'gap':
# 				if not pos == max(tRNA_struct.loc[name].index):
# 					table.loc[(table.isodecoder == name) & (table.pos >= pos), 'pos'] += 1
# 				table = table.append(new)
# 			if not any(table.loc[table.isodecoder == name].pos == pos):
# 				table = table.append(new)

# 	return(table)

def generateModsTable(sampleGroups, out_dir, threads, min_cov, mismatch_dict, insert_dict, cluster_dict, cca, remap, misinc_thresh, knownTable, tRNA_dict, Inosine_clusters, unique_isodecoderMMs, splitBool, isodecoder_sizes, clustering):
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
	tRNA_struct, tRNA_ungap2canon, cons_pos_list, cons_pos_dict = tRNAclassifier(out_dir)
	tRNA_struct_df = pd.DataFrame(tRNA_struct).unstack().rename_axis(('cluster', 'pos')).rename('struct')
	tRNA_struct_df = pd.DataFrame(tRNA_struct_df)

	# initiate custom non-daemonic multiprocessing pool and run with bam names
	pool = MyPool(multi)
	# to avoid assigning too many threads, divide available threads by number of processes
	threadsForMP = int(threads/multi)
	func = partial(bamMods_mp, out_dir, min_cov, baminfo, mismatch_dict, insert_dict, cluster_dict, cca, tRNA_struct_df, remap, misinc_thresh, knownTable, tRNA_dict, unique_isodecoderMMs, splitBool, isodecoder_sizes, threadsForMP)
	new_mods, new_Inosines = zip(*pool.map(func, bamlist))
	pool.close()
	pool.join()

	filtered = list()
	
	if not remap:

		modTable_total = pd.DataFrame()
		countsTable_total = pd.DataFrame()
		stopTable_total = pd.DataFrame()
		readthroughTable_total = pd.DataFrame()
		newMods_total = pd.DataFrame()

		dinuc_table = pd.DataFrame()
		CCAvsCC_table = pd.DataFrame()

		for bam in bamlist:
			# read in temp files and then delete
			modTable = pd.read_csv(bam + "mismatchTable.csv", header = 0, sep = "\t")
			modTable['canon_pos'] = modTable['pos'].map(cons_pos_dict)
			# edit misinc. propoportions of inosines to reflect true level of Gs
			for cluster in Inosine_clusters:
				for isodecoder in cluster_dict[cluster]:
					if any(modTable.isodecoder.str.contains(isodecoder)):
						modTable.at[(modTable.canon_pos == '34') & (modTable['type'] == 'G') & (modTable.isodecoder == isodecoder), 'proportion'] = 1 - sum(modTable[(modTable.canon_pos == '34') & (modTable['type'] != 'G') & (modTable.isodecoder == isodecoder)]['proportion'].dropna())
			os.remove(bam + "mismatchTable.csv")

			stopTable = pd.read_csv(bam + "RTstopTable.csv", header = 0, sep = "\t")
			stopTable['canon_pos'] = stopTable['pos'].map(cons_pos_dict)
			os.remove(bam + "RTstopTable.csv")

			readthroughTable = pd.read_csv(bam + "readthroughTable.csv", header = 0, sep = "\t")
			readthroughTable['canon_pos'] = readthroughTable['pos'].map(cons_pos_dict)
			os.remove(bam + "readthroughTable.csv")

			countsTable = pd.read_csv(bam + "countTable.csv", header = 0, sep = "\t")
			os.remove(bam + "countTable.csv")

			newModsTable = pd.read_csv(bam + "_predictedModstemp.csv", header = None, names = ['isodecoder', 'pos', 'identity', 'bam'], sep = "\t")
			os.remove(bam + "_predictedModstemp.csv")

			# add individual temp files to main concatenated table
			modTable_total = modTable_total.append(modTable)
			stopTable_total = stopTable_total.append(stopTable)
			readthroughTable_total = readthroughTable_total.append(readthroughTable)
			if countsTable_total.empty:
				countsTable_total = countsTable_total.append(countsTable)
			else:
				countsTable_total = pd.merge(countsTable_total, countsTable, on = "isodecoder", how = "left")
			newMods_total = newMods_total.append(newModsTable)

			if cca:
				# same for CCA analysis files
				dinuc = pd.read_csv(bam + "_dinuc.csv", header = None, keep_default_na=False, sep = "\t")
				os.remove(bam + "_dinuc.csv")
				CCA = pd.read_table(bam + "_CCAcounts.csv", header = None)
				os.remove(bam + "_CCAcounts.csv")

				dinuc_table = dinuc_table.append(dinuc)
				CCAvsCC_table = CCAvsCC_table.append(CCA)

		# plot coverage and get isodecoders to filter from mods and stops
		countsTable_total.index = countsTable_total.isodecoder
		countsTable_total.drop(columns = ['isodecoder'], inplace = True)
		filtered = filterCoverage(countsTable_total, min_cov)

		# filter and output tables
		modTable_total = modTable_total[~modTable_total.isodecoder.isin(filtered)]
		#modTable_total.drop(modTable_total[modTable_total['isodecoder'] in filtered].index, inplace = True)
		modTable_total.to_csv(out_dir + "mods/mismatchTable.csv", sep = "\t", index = False, na_rep = 'NA')
		with open(out_dir + "mods/knownModsTable.csv", "w") as known:
			known.write("cluster\tpos\n")
			for cluster, data in knownTable.items():
				for pos in data:
					known.write(cluster + "\t" + str(pos) + "\n")

		stopTable_total = stopTable_total[~stopTable_total.isodecoder.isin(filtered)]
		stopTable_total.to_csv(out_dir + "mods/RTstopTable.csv", sep = "\t", index = False, na_rep = 'NA')

		readthroughTable_total = readthroughTable_total[~readthroughTable_total.isodecoder.isin(filtered)]
		readthroughTable_total.to_csv(out_dir + "mods/readthroughTable.csv", sep = "\t", index = False, na_rep = 'NA')	
		
		# add column to counts to indicate complete isodecoder split or not
		countsTable_total['Single_isodecoder'] = "NA"
		for cluster in countsTable_total.index:
			if cluster in splitBool:
				countsTable_total.at[cluster, 'Single_isodecoder'] = "False"
			else:
				countsTable_total.at[cluster, 'Single_isodecoder'] = "True"
		countsTable_total.to_csv(out_dir + "Isodecoder_counts.txt", sep = "\t", index = True, na_rep = "0")

		# map canon_pos for each isodecoder ungapped pos to newMods
		tRNA_ungap2canon_table = pd.DataFrame.from_dict(tRNA_ungap2canon, orient = "index")
		tRNA_ungap2canon_table = tRNA_ungap2canon_table.reset_index()
		tRNA_ungap2canon_table = tRNA_ungap2canon_table.melt(var_name='pos', value_name='canon_pos', id_vars='index')
		tRNA_ungap2canon_table.columns = ['isodecoder', 'pos', 'canon_pos']
		tRNA_ungap2canon_table['pos'] = tRNA_ungap2canon_table['pos'].astype(int)
		newMods_total = pd.merge(newMods_total, tRNA_ungap2canon_table, on = ['isodecoder', 'pos'], how = "left")
		newMods_total = newMods_total[~newMods_total.isodecoder.isin(filtered)]
		# make pivot table from mods and add A, C, G, T misinc. proportions for new mods
		pivot = modTable_total.pivot_table(index = ['isodecoder', 'bam', 'canon_pos'], columns = 'type', values = 'proportion')
		pivot = pivot.reset_index()
		pivot['bam'].replace(out_dir, "", regex = True, inplace = True)
		newMods_total = pd.merge(newMods_total, pivot, on = ['isodecoder', 'canon_pos', 'bam'], how = "left")
		newMods_total.drop(columns = ['pos'], inplace = True)
		newMods_total.to_csv(out_dir + 'mods/predictedMods.csv', sep = "\t", index = False, na_rep = "NA")

		if cca:
			dinuc_table.columns = ['dinuc', 'proportion', 'sample']
			dinuc_table.to_csv(out_dir + "CCAanalysis/AlignedDinucProportions.csv", sep = "\t", index = False, na_rep = 'NA')
			CCAvsCC_table.columns = ['gene', 'end', 'sample', 'condition', 'count']
			CCAvsCC_table = CCAvsCC_table[~CCAvsCC_table.gene.isin(filtered)]
			CCAvsCC_table.to_csv(out_dir + "CCAanalysis/CCAcounts.csv", sep = "\t", index = False)

		# Anticodon and/or isodecoder counts counts
		countReads(out_dir + "Isodecoder_counts.txt", out_dir, isodecoder_sizes, clustering, tRNA_dict)

		log.info("** Read counts per isodecoder saved to " + out_dir + "counts/Isodecoder_counts.txt **")

	return(new_mods, new_Inosines, filtered)
