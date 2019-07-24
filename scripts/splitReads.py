#! /usr/bin/env python3

import pandas as pd
import logging
from collections import defaultdict

#########################################################
# Calculate read counts per isodecoder for each cluster #
#########################################################

log = logging.getLogger(__name__)

def splitReadsIsodecoder(isodecoder_counts, clusterMMprops, tRNA_dict, cluster_dict, mismatch_dict, cluster_perPos_mismatchMembers, out_dir):

	log.info("\n+--------------------------------------+\
		\n| Splitting read counts by isodecoders |\
		\n+--------------------------------------+")

	# read in counts from featureCounts
	counts = pd.read_csv(out_dir + "counts.txt", header = 0, sep = "\t", comment='#', quotechar="'")
	counts.index = counts['Geneid']
	counts = counts.drop(columns=['Geneid', 'Chr','Start','End','Strand','Length'])

	isodecoder_sizes = defaultdict(int)

	for cluster, mismatches in mismatch_dict.items():
		mismatches = sorted(mismatches, reverse = True)
		isodecoder_num = isodecoder_counts[cluster]
		curr_isodecoders = 1 # automatically start at 1 which is the cluster parent
		detected_seqs = defaultdict(list)
		detected_seqs = [tRNA_dict[cluster]['sequence']] # list of detected sequences - add parent automatically
		detected_clusters = [cluster] # list of detected/accounted for isodecoders - also automatically add parent here
		cluster_members = {tRNA:data['sequence'] for tRNA, data in tRNA_dict.items() if tRNA in cluster_dict[cluster]}
		for pos in mismatches:
			if curr_isodecoders < isodecoder_num:
				type_count = defaultdict(list)
				mismatch_members = {tRNA:sequence for tRNA, sequence in cluster_members.items() if tRNA in cluster_perPos_mismatchMembers[cluster][pos]}
				#mismatch_members = {tRNA:sequence for tRNA, sequence in cluster_members.items() if not sequence[pos].upper() == tRNA_dict[cluster]['sequence'][pos].upper()}
				for tRNA, sequence in mismatch_members.items():
					if (not sequence.upper() in detected_seqs) and (not tRNA in detected_clusters):
						type_count[sequence[pos]].append(tRNA)
						detected_seqs.append(sequence.upper())
						#isodecoder_sizes[sequence] = 
				for identity, tRNAs in type_count.items():
					if len(tRNAs) == 1:
						new_cluster_counts = list()
						detected_clusters.append(tRNAs[0])
						curr_isodecoders += 1
						for bam in clusterMMprops['bam'].unique().tolist():
							try:
								fraction = clusterMMprops[(clusterMMprops.bam == bam) & (clusterMMprops.pos == pos) & (clusterMMprops.type == identity)].loc[cluster]['proportion']
							except KeyError:
								fraction = 0 
							
							new_read_count = int(counts[bam].loc[cluster]) * fraction
							new_cluster_counts.append(new_read_count)
							remaining_read_count = int(counts[bam].loc[cluster]) - new_read_count
							counts.at[cluster, bam] = remaining_read_count
						
						counts.loc[tRNAs[0]] = new_cluster_counts

					elif len(tRNAs) > 1:
						for tRNA in tRNAs:
							sequence = tRNA_dict[tRNA]['sequence'].upper()
							detected_seqs.remove(sequence)

	print(counts.shape)
	counts.to_csv(out_dir + "Isodecoder_counts.txt", sep = "\t")
	log.info("Read counts per isodecoder saved to " + out_dir + "counts/Isodecoder_counts.txt")


