 #! /usr/bin/env python3

from __future__ import absolute_import
import pandas as pd
import numpy as np
import logging
from .ssAlign import aligntRNA
from collections import defaultdict
import pickle

#########################################################
# Calculate read counts per isodecoder for each cluster #
#########################################################

log = logging.getLogger(__name__)

def dd():
	return(defaultdict(dict))

def splitIsodecoder(tRNA_dict, cluster_dict, mismatch_dict, insert_dict, del_dict, cluster_perPos_mismatchMembers, out_dir, experiment_name):

	log.info("\n+------------------------------------------------------------------------------+\
		\n| Characterizing cluster mismatches for read splitting by unique tRNA sequence |\
		\n+------------------------------------------------------------------------------+")

	isodecoder_sizes = defaultdict(int)
	unique_isodecoderMMs = defaultdict(dd)

	with open("cluster_dict.pkl", "wb") as cluster_dict_out, open("mismatch_dict.pkl", "wb") as mismatch_dict_out, open("insert_dict.pkl", "wb") as insert_dict_out, open("cluster_perPos_mismatchMembers.pkl", "wb") as cluster_perPos_mismatchMembers_out, open("tRNA_dict.pkl", "wb") as tRNA_dict_out, open("del_dict.pkl","wb") as del_dict_out:
		pickle.dump(cluster_dict, cluster_dict_out)
		pickle.dump(mismatch_dict, mismatch_dict_out)
		pickle.dump(insert_dict, insert_dict_out)
		pickle.dump(cluster_perPos_mismatchMembers, cluster_perPos_mismatchMembers_out)
		pickle.dump(tRNA_dict, tRNA_dict_out)
		pickle.dump(del_dict, del_dict_out)

	log.info("** Assessing mismatches between cluster members and parent... **")

	for cluster, mismatches in mismatch_dict.items():
		mismatches = sorted(mismatches, reverse = True)
		curr_isodecoders = 1 # automatically start at 1 which is the cluster parent
		detected_seqs = defaultdict(list)
		detected_seqs = [tRNA_dict[cluster]['sequence']] # list of detected sequences - add parent automatically
		detected_clusters = [cluster] # list of detected/accounted for isodecoders - also automatically add parent here
		cluster_members = {tRNA:data['sequence'] for tRNA, data in tRNA_dict.items() if tRNA in cluster_dict[cluster]}
		isodecoder_num = len(set([sequences.upper() for sequences in cluster_members.values()]))
		# for each mismatch position, find sequences with unique mismatch position and type comapred to other cluster members and record mismatch type
		# use this mismatch to find misinc proportion and split reads
		for pos in mismatches:
			if curr_isodecoders < isodecoder_num: # do this until all isodecoders for cluster have been found, then stop
				type_count = defaultdict(list)
				mismatch_members = {tRNA:sequence for tRNA, sequence in cluster_members.items() if tRNA in cluster_perPos_mismatchMembers[cluster][pos]}
				for tRNA, sequence in mismatch_members.items():
					# do not process same sequence or cluster twice
					if (not sequence.upper() in detected_seqs) and (not tRNA in detected_clusters):
						# catch IndexError exception when cluster parent is longer than member and mismatch position lies outiside cluster member length
						# in these cases just ignore this specific mismatch_member tRNA as the mismatch in question does not apply to this member
						try:
							# find number of inserts before mismatch in question to ensure that the correct identity in the member is sliced by subtracting from the mismatch pos in the parent
							ins_num = len(set([ins for ins in insert_dict[cluster] if tRNA in insert_dict[cluster][ins] and ins < pos]))
							type_count[sequence[pos-ins_num]].append(tRNA)
							detected_seqs.append(sequence.upper())
						except IndexError:
							continue
				
				# Only process tRNAs which are unique in a specific mismatch at a position, thereby using this mismacth to uniquely assign reads to only this isodecoder
				for identity, tRNAs in type_count.items():
					if len(tRNAs) == 1:
						detected_clusters.append(tRNAs[0])
						isodecoder_items = [tRNA for tRNA, sequence in cluster_members.items() if sequence.upper() == tRNA_dict[tRNAs[0]]['sequence'].upper()]
						isodecoder_sizes[tRNAs[0]] = len(isodecoder_items)
						# update cluster_dict by keeping only tRNAs not in current isodecoder group
						cluster_dict[cluster] = [tRNA for tRNA in cluster_dict[cluster] if not tRNA in isodecoder_items]
						curr_isodecoders += 1
						unique_isodecoderMMs[cluster][pos][identity.upper()] = tRNAs[0]

					# Otherwise remove the sequence from detected_seqs so that it can be processed again for another mismatch position at which it might be unique
					elif len(tRNAs) > 1:
						for tRNA in tRNAs:
							sequence = tRNA_dict[tRNA]['sequence'].upper()
							detected_seqs.remove(sequence)

		# Handle insertions in cluster parents similarly to mismatches as above
		for insertion, members in insert_dict[cluster].items():
			if curr_isodecoders < isodecoder_num:
				type_count = defaultdict(list)
				insert_members = {tRNA:sequence for tRNA, sequence in cluster_members.items() if tRNA in members}
				# check if insertion is specific to one isodecoder by finding the length of the set of isodecoder numbers from insert_members 
				# (isodecoder numbers given as second last value in tRNA name - e.g. tRNA-Phe-GAA-X-1 where X = isodecoder number)
				isodecoders_withIns = len(set([num for tRNA in list(insert_members.keys()) for num in tRNA.split("-")[-2]]))
				if isodecoders_withIns == 1:
					for tRNA, sequence in insert_members.items():
						if (not sequence.upper() in detected_seqs) and (not tRNA in detected_clusters):
							type_count['insertion'].append(tRNA)
							detected_seqs.append(sequence.upper())

				for identity, tRNAs in type_count.items():
					if len(tRNAs) == 1:
						detected_clusters.append(tRNAs[0])
						isodecoder_items = [tRNA for tRNA, sequence in cluster_members.items() if sequence.upper() == tRNA_dict[tRNAs[0]]['sequence'].upper()]
						isodecoder_size = len(isodecoder_items)
						isodecoder_sizes[tRNAs[0]] = isodecoder_size
						# update cluster_dict by keeping only tRNAs not in current isodecoder group
						cluster_dict[cluster] = [tRNA for tRNA in cluster_dict[cluster] if not tRNA in isodecoder_items]
						curr_isodecoders += 1
						unique_isodecoderMMs[cluster][insertion][identity] = tRNAs[0]

					# Otherwise remove the sequence from detected_seqs so that it can be processed again for a mismatch position at which it is unique
					elif len(tRNAs) > 1:
						for tRNA in tRNAs:
							sequence = tRNA_dict[tRNA]['sequence'].upper()
							detected_seqs.remove(sequence)

	# for all clusters in cluster_dict, isodecoder size is number of members remaining after updating in above code
	# these include clusters with only one isodecoder - i.e. not in mismatch_dict, and clusters that could not be separated into isodecoders because no unique mismatch distinguishes them
	splitBool = list()
	for cluster, members in cluster_dict.items():
		cluster_size = len(members)
		isodecoder_sizes[cluster] = cluster_size
		remaining_isodecoders = set([data['sequence'].upper() for member, data in tRNA_dict.items() if member in members])
		if len(remaining_isodecoders) > 1:
			splitBool.append("-".join(cluster.split("-")[:-1]))

	# save isodecoder info for DESeq2
	total_detected_isodecoders = 0
	with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
		isodecoderInfo.write("Isodecoder\tsize\n")
		for isodecoder, size in isodecoder_sizes.items():
			isodecoder = "-".join(isodecoder.split("-")[:-1]) if not "chr" in isodecoder else isodecoder
			isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")
			if not isodecoder in splitBool:
				total_detected_isodecoders += 1

	# write isodecoder fasta for alignment and context analysis
	with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
		for seq in isodecoder_sizes.keys():
			shortname = "-".join(seq.split("-")[:-1]) if not "chr" in seq else seq
			tempSeqs.write(">" + shortname + "\n" + tRNA_dict[seq]['sequence'] + "\n")
	aligntRNA(tempSeqs.name, out_dir)

	total_isodecoders = len(set([data["sequence"].upper() for tRNA,data in tRNA_dict.items()]))

	log.info("Total unique tRNA sequenes in input: {}".format(total_isodecoders))
	log.info("Total deconvoluted unique tRNA sequences: {}".format(total_detected_isodecoders))

	return(unique_isodecoderMMs, splitBool, isodecoder_sizes)

def writeIsodecoderTranscripts(out_dir, experiment_name, cluster_dict, tRNA_dict):
	# write isodecoderTransripts.fa when cluster_id == 1 to avoid issues with shortened isodecoder names and output files

	with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
		for seq in cluster_dict.keys():
			shortname = "-".join(seq.split("-")[:-1]) if not "chr" in seq else seq
			tempSeqs.write(">" + shortname + "\n" + tRNA_dict[seq]['sequence'] + "\n")
	aligntRNA(tempSeqs.name, out_dir)

def getIsodecoderSizes(out_dir, experiment_name, tRNAdict):
	# get isodecoder sizes for tRNA sequences - useful for when clustering is disabled and above function is not applicable

	isodecoder_sizes = defaultdict(int)
	already_added = set()
	for tRNA in tRNAdict:
		if tRNA not in already_added:
			sameSeq = [tRNAs for tRNAs, data in tRNAdict.items() if data['sequence'] == tRNAdict[tRNA]['sequence']]
			already_added.update(sameSeq)
			isodecoder_sizes[tRNA] = len(sameSeq)

	with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
		isodecoderInfo.write("Isodecoder\tsize\n")
		for isodecoder, size in isodecoder_sizes.items():
			isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	return(isodecoder_sizes)
