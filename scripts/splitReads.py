 #! /usr/bin/env python3

import pandas as pd
import numpy as np
import logging
from collections import defaultdict

#########################################################
# Calculate read counts per isodecoder for each cluster #
#########################################################

log = logging.getLogger(__name__)

def splitReadsIsodecoder(isodecoder_counts, clusterMMprops, tRNA_dict, cluster_dict, mismatch_dict, insert_dict, cluster_perPos_mismatchMembers, out_dir, experiment_name):

	log.info("\n+--------------------------------------+\
		\n| Splitting read counts by isodecoders |\
		\n+--------------------------------------+")

	# read in counts from featureCounts
	counts = pd.read_csv(out_dir + "counts.txt", header = 0, sep = "\t", comment='#', quotechar="'")
	counts.index = counts['Geneid']
	counts = counts.drop(columns=['Geneid', 'Chr','Start','End','Strand','Length'])

	isodecoder_sizes = defaultdict(int)
	total_isodecoders = 0

	for cluster, mismatches in mismatch_dict.items():
		mismatches = sorted(mismatches, reverse = True)
		isodecoder_num = isodecoder_counts[cluster] 
		total_isodecoders += isodecoder_num
		curr_isodecoders = 1 # automatically start at 1 which is the cluster parent
		detected_seqs = defaultdict(list)
		detected_seqs = [tRNA_dict[cluster]['sequence']] # list of detected sequences - add parent automatically
		detected_clusters = [cluster] # list of detected/accounted for isodecoders - also automatically add parent here
		cluster_members = {tRNA:data['sequence'] for tRNA, data in tRNA_dict.items() if tRNA in cluster_dict[cluster]}
		parent_size = len([tRNA for tRNA, sequence in cluster_members.items() if sequence == tRNA_dict[cluster]['sequence']]) # number of sequences identical to parent - i.e. size of parent isodecoder group
		
		# for each mismatch position, find sequences with unique mismatch position and type comapred to other cluster members and record mismatch type
		# use this msmatch to find misinc proportion and split reads
		for pos in mismatches:
			if curr_isodecoders < isodecoder_num: # do this until all isodecoders for cluster have been found, then stop
				type_count = defaultdict(list)
				mismatch_members = {tRNA:sequence for tRNA, sequence in cluster_members.items() if tRNA in cluster_perPos_mismatchMembers[cluster][pos]}
				for tRNA, sequence in mismatch_members.items():
					# do not process same sequence or cluster twice
					if (not sequence.upper() in detected_seqs) and (not tRNA in detected_clusters):
						type_count[sequence[pos]].append(tRNA)
						detected_seqs.append(sequence.upper())
				
				# Only process tRNAs which are unique in a specific mismatch at a position, thereby using this mismacth to uniquely assign reads to only this isodecoder
				for identity, tRNAs in type_count.items():
					if len(tRNAs) == 1:
						new_cluster_counts = list()
						detected_clusters.append(tRNAs[0])
						isodecoder_items = [tRNA for tRNA, sequence in cluster_members.items() if sequence == tRNA_dict[tRNAs[0]]['sequence']]
						isodecoder_size = len(isodecoder_items)
						isodecoder_sizes[tRNAs[0]] = isodecoder_size
						# update cluster_dict by keeping only tRNAs not in current isodecoder group
						cluster_dict[cluster] = [tRNA for tRNA in cluster_dict[cluster] if not tRNA in isodecoder_items]
						curr_isodecoders += 1
						for bam in clusterMMprops['bam'].unique().tolist():
							try:
								fraction = clusterMMprops[(clusterMMprops.bam == bam) & (clusterMMprops.pos == pos) & (clusterMMprops.type == identity)].loc[cluster]['proportion']
							except KeyError:
								fraction = 0 
							
							new_read_count = round(int(counts[bam].loc[cluster]) * fraction)
							new_cluster_counts.append(new_read_count)
							remaining_read_count = int(counts[bam].loc[cluster]) - new_read_count
							counts.at[cluster, bam] = remaining_read_count
						
						counts.loc[tRNAs[0]] = new_cluster_counts
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
				for tRNA, sequence in insert_members.items():
					if (not sequence.upper() in detected_seqs) and (not tRNA in detected_clusters):
						type_count['insertion'].append(tRNA)
						detected_seqs.append(sequence.upper())

				for identity, tRNAs in type_count.items():
					if len(tRNAs) == 1:
						new_cluster_counts = list()
						detected_clusters.append(tRNAs[0])
						isodecoder_items = [tRNA for tRNA, sequence in cluster_members.items() if sequence == tRNA_dict[tRNAs[0]]['sequence']]
						isodecoder_size = len(isodecoder_items)
						isodecoder_sizes[tRNAs[0]] = isodecoder_size
						# update cluster_dict by keeping only tRNAs not in current isodecoder group
						cluster_dict[cluster] = [tRNA for tRNA in cluster_dict[cluster] if not tRNA in isodecoder_items]
						curr_isodecoders += 1
						for bam in clusterMMprops['bam'].unique().tolist():
							try:
								fraction = clusterMMprops[(clusterMMprops.bam == bam) & (clusterMMprops.pos == insertion) & (clusterMMprops.type == identity)].loc[cluster]['proportion']
							except KeyError:
								fraction = 0 
							
							new_read_count = round(int(counts[bam].loc[cluster]) * fraction)
							new_cluster_counts.append(new_read_count)
							remaining_read_count = int(counts[bam].loc[cluster]) - new_read_count
							counts.at[cluster, bam] = remaining_read_count
						
						counts.loc[tRNAs[0]] = new_cluster_counts
					# Otherwise remove the sequence from detected_seqs so that it can be processed again for a mismatch position at which it is unique
					elif len(tRNAs) > 1:
						for tRNA in tRNAs:
							sequence = tRNA_dict[tRNA]['sequence'].upper()
							detected_seqs.remove(sequence)

	
	counts['Single_isodecoder'] = "NA"
	# for all clusters in cluster_dict, isdecoder size is number of members remaining after updating in above code
	# these include clusters with only one isodecoder - i.e. not in mismatch_dict, and clusters that could not be separated into isodecoders because no unique mismatch distinguishes them
	for cluster, members in cluster_dict.items():
		cluster_size = len(members)
		isodecoder_sizes[cluster] = cluster_size
		remaining_isodecoders = set([data['sequence'].upper() for member, data in tRNA_dict.items() if member in members])
		if len(remaining_isodecoders) > 1:
			counts.at[cluster, 'Single_isodecoder'] = "False"
		# if cluster not in mismatch_dict:
		# 	cluster_members = [gene for gene, data in tRNA_dict.items() if gene in cluster_dict[cluster]]
		# 	parent_size = len(cluster_members)
		# 	isodecoder_sizes[cluster] = parent_size
		# else:
		# 	#for member in cluster_dict[cluster]:
		# 	#	cluster_uniqeseqs.add(tRNA_dict[member]['sequence'].upper())
		# 	cluster_size = len(cluster_dict)
		# 	isodecoder_sizes[cluster] = cluster_size

	# save isodecoder info for DESeq2
	with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
		isodecoderInfo.write("Isodecoder\tsize\n")
		for isodecoder, size in isodecoder_sizes.items():
			isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	counts.to_csv(out_dir + "Isodecoder_counts.txt", sep = "\t", na_rep = "NA")
	log.info("Read counts per isodecoder saved to " + out_dir + "counts/Isodecoder_counts.txt")


