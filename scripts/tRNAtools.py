#!/usr/bin/env python3

##################################################################################
# Utilities for tRNA modification parsing, transcript building, and SNP indexing #
##################################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import re, copy, sys, os, shutil, subprocess, logging
from pathlib import Path
import urllib.request
from collections import defaultdict
import ssAlign

log = logging.getLogger(__name__)

# This function is used to create nested defaultdicts like that needed for tRNA_dict
# Note this can be done with lambda functions (e.g. tRNA_dict = defaultdict(lambda: defaultdict()))
# But lambda functions cannot be pickled, and pickling is required for parallelization with multiprocessing. tRNA_dict is passed to such a multiprocessing pool later
def dd():
	return(defaultdict())

def tRNAparser (gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, posttrans_mod_off):
# tRNA sequence files parser and dictionary building

	# Generate modification reference table
	modifications = modificationParser(modifications_table)
	temp_name = gtRNAdb.split("/")[-1]
                
	log.info("\n+" + ("-" * (len(temp_name)+24)) + "+\
		\n| Starting analysis for {} |\
		\n+".format(temp_name) + ("-" * (len(temp_name)+24)) + "+")
	      
	log.info("Processing tRNA sequences...")

	# Build dictionary of sequences from gtRNAdb fasta
	tRNA_dict = defaultdict(dd)
	temp_dict = SeqIO.to_dict(SeqIO.parse(gtRNAdb,"fasta"))

	# Initialise intron dictionary
	Intron_dict = initIntronDict(tRNAscan_out)

	species = set()
	for seq in temp_dict:
		# Get species of input tRNA seqs to subset full Modomics table
		species.add(' '.join(seq.split('_')[0:2]))
		# only add to dictionary if not nmt or undetermined sequence
		if not (re.search('Und', seq) or re.search('nmt', seq)):
			tRNAseq = intronRemover(Intron_dict, temp_dict, seq, posttrans_mod_off)
			# if re.search('nmt', seq):
			# 	seq = seq.split("-")[0] + "_" + "-".join(seq.split("-")[1:])
			# 	loc_type = "mitochondrial"
			# else:
			loc_type = "cytosolic"
			tRNA_dict[seq]['sequence'] = tRNAseq
			tRNA_dict[seq]['species'] = ' '.join(seq.split('_')[0:2])
			tRNA_dict[seq]['type'] = loc_type

	# add mitochondrial tRNAs if given
	if mitotRNAs:
		temp_dict = SeqIO.to_dict(SeqIO.parse(mitotRNAs,"fasta"))
		mito_count = defaultdict(int)
		# read each mito tRNA, edit sequence header to match nuclear genes as above and add to tRNA_dict
		for seq in temp_dict:
			seq_parts = seq.split("|")
			species = ' '.join(seq_parts[1].split("_")[0:])
			anticodon = seq_parts[4]
			amino = re.search("[a-zA-z]+", seq_parts[3]).group(0)
			mito_count[anticodon] += 1
			new_seq = seq_parts[1] + "_mito_tRNA-" + amino + "-" + seq_parts[4] + "-1-" + str(mito_count[anticodon])
			tRNAseq = str(temp_dict[seq].seq) + "CCA"
			tRNA_dict[new_seq]['sequence'] = tRNAseq
			tRNA_dict[new_seq]['type'] = 'mitochondrial'
			tRNA_dict[new_seq]['species'] = ' '.join(seq.split('_')[0:2])

		num_cytosilic = len([k for k in tRNA_dict.keys() if tRNA_dict[k]['type'] == "cytosolic"])
		num_mito = len([k for k in tRNA_dict.keys() if tRNA_dict[k]['type'] == "mitochondrial"])

		log.info("{} cytosolic and {} mitochondrial tRNA sequences imported".format(num_cytosilic, num_mito))

		# if 'nmt' in seq:
		# 	tRNA_dict[seq]['type'] = 'mitochondrial'
		# 	new_seq = str(seq.split('-')[0]) + '_' + '-'.join(seq.split('-')[1:])
		# 	tRNA_dict[new_seq] = tRNA_dict[seq]
		# 	del tRNA_dict[seq]
		# else:
		# 	tRNA_dict[seq]['type'] = 'cytosolic'


	# Read in and parse modomics file to contain similar headers to tRNA_dict
	# Save in new dict

	log.info("Processing modomics database...")
	modomics_file = getModomics()
	modomics_dict = {}
	for line in modomics_file.splitlines():
		line = line.strip()
		sameIDcount = 0

		if line.startswith('>'):
			line = line.replace(' | ','|')
			# Check if species matches those from gtRNAdb input, otherwise skip entry
			mod_species = line.split('|')[3]
			if not mod_species in species:
				continue
			else:
				# Replace modomics antidon with normal ACGT codon, keep "." for pattern matching
				anticodon = str(line.split('|')[2])
				new_anticodon = getUnmodSeq(anticodon, modifications)
				if "N" in new_anticodon:
					continue
				#new_anticodon = new_anticodon.replace("N", ".")

				# Check amino acid name in modomics - set to iMet if equal to Ini to match gtRNAdb
				amino = str(line.split('|')[1])
				if amino == 'Ini' :
					amino = 'iMet'

				curr_id = str(line.split('|')[3].split(' ')[0]) + '_' + str(line.split('|')[3].split(' ')[1]) + '_' + str(line.split('|')[0].split(' ')[1]) + '-' + amino + '-' + new_anticodon
				#Unique names for duplicates
				if curr_id in modomics_dict:
					sameIDcount += 1
					curr_id = curr_id + '-' + str(sameIDcount)

				tRNA_type = str(line.split('|')[4])
				modomics_dict[curr_id] = {'sequence':'','type':tRNA_type, 'anticodon':new_anticodon}
		
		else:
			if not mod_species in species:
				continue
			else:
				if "N" in new_anticodon:
					continue

				sequence = line.strip().replace('U','T').replace('-','')
				modomics_dict[curr_id]['sequence'] = sequence
				unmod_sequence = getUnmodSeq(sequence, modifications)

				# Return list of modified nucl indices and add to modomics_dict
				# add unmodified seq to modomics_dict by lookup to modifications
				nonMod = ['A','C','G','T','-']
				modPos = [i for i, x in enumerate(modomics_dict[curr_id]['sequence']) if x not in nonMod]
				inosinePos = [i for i, x in enumerate(modomics_dict[curr_id]['sequence']) if x == 'I']
				modomics_dict[curr_id]['modified'] = modPos
				modomics_dict[curr_id]['unmod_sequence'] = unmod_sequence
				modomics_dict[curr_id]["InosinePos"] = inosinePos

				# If 'N' in anticodon then duplicate entry 4 times for each possibility
				# anticodon = curr_id.split('-')[2]
				# if 'N' in anticodon:
				# 	for rep in ['A','C','G','T']:
				# 		duplicate_item = str(curr_id.split('-')[0]) + '-' + str(curr_id.split('-')[1]) + '-' + str(anticodon.replace('N', rep))
				# 		duplicate_unmod_seq = modomics_dict[curr_id]['unmod_sequence'].replace('N',rep)
				# 		duplicate_anticodon = modomics_dict[curr_id]['anticodon'].replace('N', rep)
				# 		modomics_dict[duplicate_item] = copy.deepcopy(modomics_dict[curr_id])
				# 		modomics_dict[duplicate_item]['unmod_sequence'] = duplicate_unmod_seq
				# 		modomics_dict[duplicate_item]['anticodon'] = duplicate_anticodon
				# 	del modomics_dict[curr_id]
				# else:
				# 	duplicate_item = curr_id

	log.info('Number of Modomics entries for species of interest: {}'.format(len(modomics_dict)))

	return(tRNA_dict,modomics_dict, species)

def getModomics():
	# Get full Modomics modified tRNA data from web

	with urllib.request.urlopen('http://modomics.genesilico.pl/sequences/list/?type_field=tRNA&subtype=all&species=all&display_ascii=Display+as+ASCII&nomenclature=abbrev') as response:
		modomics = response.read().decode()

	return modomics


def modsToSNPIndex(gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, experiment_name, out_dir, snp_tolerance = False, cluster = False, cluster_id = 0.95, posttrans_mod_off = False):
# Builds SNP index needed for GSNAP based on modificaiton data for each tRNA and clusters tRNAs

	nomatch_count = 0
	match_count = 0
	total_snps = 0
	total_inosines = 0
	snp_records = list()
	seq_records = defaultdict()
	anticodon_list = list()
	tRNAbed = open(out_dir + experiment_name + "_maturetRNA.bed","w")
	# generate modomics_dict and tRNA_dict
	tRNA_dict, modomics_dict, species = tRNAparser(gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, posttrans_mod_off)
	temp_dir = out_dir + "/tmp/"

	try:
		os.mkdir(temp_dir)
	except FileExistsError:
		log.warning("Temp folder present - previous run interrupted? Overwriting old temp files...\n")

	###################################################
	## Main code for matching and SNP index building ##
	###################################################

	# remove spurious CCA seen in some sequences (i.e. genomically encoded or artifact from gtRNAdb or tRNAScan-SE prediction)
	# to do this, remove previously added CCA (intronRemover()) that falls out of canonical structure using ssAlign module

	with open(out_dir + 'temptRNAseqs.fa', 'w') as tempSeqs:
		for seq in tRNA_dict:
			tempSeqs.write(">" + seq + "\n" + tRNA_dict[seq]['sequence'] + "\n")

	ssAlign.aligntRNA(tempSeqs.name, out_dir)
	extra_cca = ssAlign.extraCCA()

	for record in extra_cca:
		tRNA_dict[record]['sequence'] = tRNA_dict[record]['sequence'][:-3]

	os.remove(tempSeqs.name)

	# match each sequence in tRNA_dict to value in modomics_dict using BLAST

	log.info("\n+------------------------+ \
		\n| Beginning SNP indexing |\
		\n+------------------------+")	

	for seq in tRNA_dict:
		# Initialise list of modified sites for each tRNA
		tRNA_dict[seq]['modified'] = []
		tRNA_dict[seq]['InosinePos'] = []
		temp_tRNAFasta = open(temp_dir + seq + ".fa","w")
		temp_tRNAFasta.write(">" + seq + "\n" + tRNA_dict[seq]['sequence'] + "\n")
		temp_tRNAFasta.close()
		tRNA_dict[seq]['anticodon'] = anticodon = re.search('.*tRNA-.*?-(.*?)-', seq).group(1)
		if not anticodon in anticodon_list:
			anticodon_list.append(anticodon)
		# find initial possible matches to modomics where anticodons match and types are the same (here regex is used to match anticodons with "." in modomics to all possible matching sequences from input tRNAs)
		match = {k:v for k,v in modomics_dict.items() if re.match("^" + v['anticodon'] + "+$", anticodon) and tRNA_dict[seq]['type'] == v['type']}
		if len(match) >= 1:
			temp_matchFasta = open(temp_dir + "modomicsMatch.fasta","w")
			for i in match:	
				temp_matchFasta.write(">" + i + "\n" + match[i]['unmod_sequence'] + "\n")
			temp_matchFasta.close()

			#blast
			blastn_cline = NcbiblastnCommandline(query = temp_tRNAFasta.name, subject = temp_matchFasta.name, task = 'blastn-short', out = temp_dir + "blast_temp.xml", outfmt = 5)
			blastn_cline()

			#parse XML result and store hit with highest bitscore	
			blast_record = NCBIXML.read(open(temp_dir + "blast_temp.xml","r"))
			maxbit = 0
			tophit = ''
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					if (hsp.bits > maxbit) and (hsp.align_length / alignment.length >= 0.9):
						maxbit = hsp.bits
						tophit = alignment.title.split(' ')[0]
			
			# return list of all modified positions for the match as long as there is only 1, add to tRNA_dict
			if tophit:
				match_count += 1
				tRNA_dict[seq]['modified'] = match[tophit]['modified']
				tRNA_dict[seq]['InosinePos'] = match[tophit]['InosinePos']
			elif len(tophit) == 0:
				nomatch_count += 1
		if len(match) == 0:
			nomatch_count += 1

		# Build seqrecord list for writing
		seq_records[str(seq)] = SeqRecord(Seq(tRNA_dict[seq]['sequence'].upper(), Alphabet.generic_dna), id = str(seq))

		tRNAbed.write(seq + "\t0\t" + str(len(tRNA_dict[seq]['sequence'])) + "\t" + seq + "\t1000\t+\n" )

	tRNAbed.close()

	log.info("{} total tRNA gene sequences (undetermined and nmt sequences excluded)".format(len(tRNA_dict)))
	log.info("{} sequences with a match to Modomics dataset".format(match_count))

	with open(str(out_dir + experiment_name + '_tRNATranscripts.fa'), "w") as temptRNATranscripts:
		SeqIO.write(seq_records.values(), temptRNATranscripts, "fasta")

	# if clustering is not activated then write full gff and report on total SNPs written 
	if not cluster:
		coverage_bed = tRNAbed.name
		mod_lists = dict()
		Inosine_lists = dict()
		with open(out_dir + experiment_name + "_tRNA.gff","w") as tRNAgff, open(out_dir + experiment_name + "isoacceptorInfo.txt","w") as isoacceptorInfo:	
			isoacceptor_dict = defaultdict(int)
			isoacceptorInfo.write("Isoacceptor\tsize\n")
			for seq in tRNA_dict:
				mod_lists[seq] = tRNA_dict[seq]['modified']
				Inosine_lists[seq] = tRNA_dict[seq]['InosinePos']
				tRNAgff.write(seq + "\ttRNAseq\texon\t1\t" + str(len(tRNA_dict[seq]['sequence'])) + "\t.\t+\t0\tgene_id '" + seq + "'\n")
				isoacceptor_group = '-'.join(seq.split("-")[:-2])
				isoacceptor_dict[isoacceptor_group] += 1
			for key, value in isoacceptor_dict.items():
				isoacceptorInfo.write(key + "\t" + str(value) + "\n")
		# generate Stockholm alignment file for all tRNA transcripts and parse additional mods file
		ssAlign.aligntRNA(temptRNATranscripts.name, out_dir)
		additionalMods, additionalInosines = additionalModsParser(species, out_dir)
		# add additional SNPs from extra file to list of modified positions, and ensure non-redundancy with set()
		# index SNPs
		# Format for SNP index (space separated):
		# >snpID chromosomeName:position(1-based) RefMod
		# e.g. >rs111 Homo_sapiens_nmt_tRNA-Leu-TAA-1-1_exp0:29 GN
		for seq in mod_lists:
			isodecoder = "-".join(seq.split("-")[1:3])
			additionalMods_sub = {k:v for k, v in additionalMods.items() if k == isodecoder and tRNA_dict[seq]['species'] in v['species']}
			if additionalMods_sub:
				tRNA_dict[seq]['modified'] = list(set(tRNA_dict[seq]['modified'] + additionalMods_sub[isodecoder]['mods']))
				mod_lists[seq] = list(set(mod_lists[seq] + additionalMods_sub[isodecoder]['mods']))

			total_snps += len(mod_lists[seq])

			# Build snp_records as before but with cluster names and non-redundant sets of modifications
			# Position is 1-based for iit_store i.e. pos + 1
			for (index, pos) in enumerate(mod_lists[seq]):
				snp_records.append(">" + seq + "_snp" + str(index) + " " + seq + ":" + str(pos + 1) + " " + tRNA_dict[seq]['sequence'][pos].upper() + "N")

		for seq in Inosine_lists:
			isodecoder = "-".join(seq.split("-")[1:3])
			additionalInosines_sub = {k:v for k, v in additionalInosines.items() if k == isodecoder and tRNA_dict[seq]['species'] in v['species']}
			if additionalInosines_sub:
				tRNA_dict[seq]['InosinePos'] = list(set(tRNA_dict[seq]['InosinePos'] + additionalInosines_sub[isodecoder]['InosinePos']))
				Inosine_lists[seq] = list(set(Inosine_lists[seq] + additionalInosines_sub[isodecoder]['InosinePos']))

			total_inosines += len(Inosine_lists[seq])

		# edit ref seqs A to G at inosine positions
		for seq in Inosine_lists:
			for pos in Inosine_lists[seq]:
				seq_records[seq].seq = seq_records[seq].seq[0:pos] + "G" + seq_records[seq].seq[pos+1:]

		with open(str(out_dir + experiment_name + '_tRNATranscripts.fa'), "w") as temptRNATranscripts:
			SeqIO.write(seq_records.values(), temptRNATranscripts, "fasta")

		if total_snps == 0:
			snp_tolerance = False

		log.info("{:,} modifications written to SNP index".format(total_snps))
		log.info("{:,} A to G replacements in reference sequences for inosine modifications".format(total_inosines))
		# empty mismatch dict to avoid error when returning it from this function
		mismatch_dict = defaultdict(list)

	##########################
	# Cluster tRNA sequences #
	##########################

	elif cluster:

		log.info("**** Clustering tRNA sequences ****")
		log.info("Clustering tRNA sequences by {:.0%} similarity...".format(cluster_id))
		# dictionary of final centroid sequences
		final_centroids = defaultdict()
		# get dictionary of sequences for each anticodon and write to fastas
		for anticodon in anticodon_list:
			seq_set = {k:{'sequence':v['sequence'],'modified':v['modified']} for k,v in tRNA_dict.items() if v['anticodon'] == anticodon}
			with open(temp_dir + anticodon + "_allseqs.fa","w") as anticodon_seqs:
				for sequence in seq_set:
					anticodon_seqs.write(">" + sequence + "\n" + seq_set[sequence]['sequence'] + "\n")
			# run usearch on each anticodon sequence fatsa to cluster
			#cluster_cmd = "usearch -cluster_fast " + temp_dir + anticodon + "_allseqs.fa -id " + str(cluster_id) + " -sizeout -centroids " + temp_dir + anticodon + "_centroids.fa -uc " + temp_dir + anticodon + "_clusters.uc &> /dev/null" 
			cluster_cmd = "usearch -cluster_fast " + temp_dir + anticodon + "_allseqs.fa -sort length -id " + str(cluster_id) + " -centroids " + temp_dir + anticodon + "_centroids.fa -uc " + temp_dir + anticodon + "_clusters.uc &> /dev/null"
			subprocess.call(cluster_cmd, shell = True)
			# sort clusters by size (i.e. number of members in cluster)
			#sort_cmd = "usearch -sortbysize " + temp_dir + anticodon + "_centroids.fa -fastaout " + temp_dir + anticodon + "_centroids_sort.fa &> /dev/null"
			#subprocess.call(sort_cmd, shell = True)
			# recluster based on sorted by size clusters
			#final_cluster_cmd = "usearch -cluster_smallmem " + temp_dir + anticodon + "_centroids_sort.fa -sortedby size -id " + str(cluster_id) + " -centroids " + temp_dir + anticodon + "_centroidsFinal.fa &> /dev/null" 
			#subprocess.call(final_cluster_cmd, shell = True)
		# combine centroids files into one file
		combine_cmd = "cat " + temp_dir + "*_centroids.fa > " + temp_dir + "all_centroids.fa"
		subprocess.call(combine_cmd, shell = True)
		centroids = SeqIO.parse(temp_dir + "all_centroids.fa", "fasta")
		for centroid in centroids:
			centroid.id = centroid.id.split(";")[0]
			final_centroids[centroid.id] = SeqRecord(Seq(str(centroid.seq).upper(), Alphabet.generic_dna), id = centroid.id) 

		# read cluster files, get nonredudant set of mod positions of all members of a cluster, create snp_records for writing SNP index
		cluster_pathlist = Path(temp_dir).glob("**/*_clusters.uc")
		mod_lists = dict() # stores non-redundant sets of mismatches and mod positions for clusters
		Inosine_lists = dict() # stores positions of inosines for clusters
		snp_records = list()
		cluster_dict = dict() # info about clusters
		mismatch_dict = defaultdict(list) # dictionary of mismatches only (not mod positions - required for misincorporation analysis in mmQuant)
		cluster_num = 0
		total_snps = 0
		total_inosines = 0
		clusterbed = open(out_dir + experiment_name + "_clusters.bed","w")
		coverage_bed = clusterbed.name
		clustergff = open(out_dir + experiment_name + "_tRNA.gff","w")
		for path in cluster_pathlist:
			with open(path,"r") as cluster_file:
				for line in cluster_file:
					line = line.strip()

					# Handle cluster centroids and initialise modified positions list
					if line.split("\t")[0] == "S":
						cluster_num += 1
						cluster_name = line.split("\t")[8].split(";")[0]
						mod_lists[cluster_name] = tRNA_dict[cluster_name]["modified"]
						Inosine_lists[cluster_name] = tRNA_dict[cluster_name]['InosinePos']
						clusterbed.write(cluster_name + "\t0\t" + str(len(tRNA_dict[cluster_name]['sequence'])) + "\t" + cluster_name + "\t1000\t+\n" )
						clustergff.write(cluster_name + "\ttRNAseq\texon\t1\t" + str(len(tRNA_dict[cluster_name]['sequence'])) + "\t.\t+\t0\tgene_id '" + cluster_name + "'\n")
						cluster_dict[cluster_name] = cluster_num
				
					# Handle members of clusters
					elif line.split("\t")[0] == "H":
						member_name = line.split("\t")[8].split(";")[0]
						cluster_name = line.split("\t")[9].split(";")[0]
						compr_aln = line.split("\t")[7]
						# if member of cluster is 100% identical (i.e. "=" in 8th column of cluster file)
						if compr_aln == "=":
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + tRNA_dict[member_name]["modified"]))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + tRNA_dict[member_name]["InosinePos"]))
							cluster_dict[member_name] = cluster_dict[cluster_name]
						
						# if there are insertions or deletions in the centroid, edit member or centroid sequences to ignore these positions
						# and edit modified positions list in order to make non-redundant positions list, similar to next else statement
						elif re.search("[ID]", compr_aln):
							cluster_seq = tRNA_dict[cluster_name]["sequence"]
							member_seq = tRNA_dict[member_name]["sequence"]
							pos = 0
							adjust_pos_del = 0
							adjust_pos_ins = 0
							insertion_pos = list()
							deletion_pos = list()
							aln_list = re.split('(.*?[A-Z])', compr_aln)
							aln_list = list(filter(None, aln_list))

							for phrase in aln_list:
								
								if ("M" in phrase) and (phrase.split("M")[0] != ""):
									pos += int(phrase.split("M")[0])
								elif ("M" in phrase) and (phrase.split("M")[0] == ""):
									pos += 1

								if ("I" in phrase) and (phrase.split("I")[0] != ""):
									insert_len = int(phrase.split("I")[0])

									for i in range(insert_len):
										insertion_pos.append(pos+i)
										#pos += 1
								elif ("I" in phrase) and (phrase.split("I")[0] == ""):
									insertion_pos.append(pos)
									#pos += 1

								if ("D" in phrase) and (phrase.split("D")[0] != ""):
									delete_len = int(phrase.split("D")[0])

									for i in range(delete_len):
										deletion_pos.append(pos)
										pos += 1
								elif ("D" in phrase) and (phrase.split("D")[0] == ""):
									deletion_pos.append(pos)
									pos += 1

							for delete in deletion_pos:
								# if delete in tRNA_dict[member_name]["modified"]:
								# 	tRNA_dict[member_name]["modified"].remove(delete)
								new_delete = delete + adjust_pos_del
								member_seq = member_seq[ :new_delete] + member_seq[new_delete+1: ]
								adjust_pos_del -= 1

							for index, insert in enumerate(insertion_pos):
								adjust_pos_len = 0
								# if insert in mod_lists[cluster_name]:
								# 	mod_lists[cluster_name].remove(insert)
								for delete in deletion_pos:
									if delete < insert:
										adjust_pos_len -= 1
								if index != 0:
									if insert == insertion_pos[index-1] + 1:
										adjust_pos_ins -= 1
								new_insert = insert + adjust_pos_len + adjust_pos_ins
								cluster_seq = cluster_seq[ :new_insert] + cluster_seq[new_insert+1: ]

							mismatches = [i for i in range(len(member_seq)) if member_seq[i].upper() != cluster_seq[i].upper()]
							mismatch_dict[cluster_name] = list(set(mismatch_dict[cluster_name] + mismatches))
							member_mods = list(set(tRNA_dict[member_name]["modified"] + mismatches))
							member_Inosines = tRNA_dict[member_name]["InosinePos"]
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + member_mods))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + member_Inosines))
							cluster_dict[member_name] = cluster_dict[cluster_name]

						# handle members that are not exact sequence matches but have no indels either
						# find mismatches and build non-redundant set
						else:
							cluster_seq = tRNA_dict[cluster_name]["sequence"]
							member_seq = tRNA_dict[member_name]["sequence"]
							mismatches = [i for i in range(len(member_seq)) if member_seq[i].upper() != cluster_seq[i].upper()]
							mismatch_dict[cluster_name] = list(set(mismatch_dict[cluster_name] + mismatches))
							member_mods = list(set(tRNA_dict[member_name]["modified"] + mismatches))
							member_Inosines = tRNA_dict[member_name]["InosinePos"]
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + member_mods))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + member_Inosines))
							cluster_dict[member_name] = cluster_dict[cluster_name]
		
		clusterbed.close()

		# Write cluster information to tsv
		with open(out_dir + experiment_name + "clusterInfo.txt","w") as clusterInfo, open(out_dir + experiment_name + "isoacceptorInfo.txt","w") as isoacceptorInfo:
			isoacceptor_dict = defaultdict(int)
			clusterInfo.write("tRNA\tcluster_num\tcluster_size\n")
			isoacceptorInfo.write("Isoacceptor\tsize\n")
			for key, value in cluster_dict.items():
				clusterInfo.write("{}\t{}\t{}\n".format(key, value, sum(clusters == value for clusters in cluster_dict.values())))
				isoacceptor_group = '-'.join(key.split("-")[:-2])
				isoacceptor_dict[isoacceptor_group] += 1
			for key, value in isoacceptor_dict.items():
				isoacceptorInfo.write(key + "\t" + str(value) + "\n")

		with open(str(out_dir + experiment_name + '_clusterTranscripts.fa'), "w") as clusterTranscripts:
			SeqIO.write(final_centroids.values(), clusterTranscripts, "fasta")

		# generate Stockholm alignment file for cluster transcripts and process additional mods file
		ssAlign.aligntRNA(clusterTranscripts.name, out_dir)
		additionalMods, additionalInosines = additionalModsParser(species, out_dir)

		log.info("{} clusters created from {} tRNA sequences".format(cluster_num,len(tRNA_dict)))

		# update mod_lists with additional mods and write SNP index
		for cluster in mod_lists:
			isodecoder = "-".join(cluster.split("-")[1:3])
			additionalMods_sub = {k:v for k, v in additionalMods.items() if k == isodecoder and tRNA_dict[cluster]['species'] in v['species']}
			if additionalMods_sub:
				tRNA_dict[cluster]['modified'] = list(set(tRNA_dict[cluster]['modified'] + additionalMods_sub[isodecoder]['mods']))
				mod_lists[cluster] = list(set(mod_lists[cluster] + additionalMods_sub[isodecoder]['mods']))

			total_snps += len(mod_lists[cluster])

			for (index, pos) in enumerate(mod_lists[cluster]):
				# Build snp_records as before but with cluster names and non-redundant sets of modifications
				# Position is 1-based for iit_store i.e. pos + 1
				snp_records.append(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "N")

		for cluster in Inosine_lists:
			isodecoder = "-".join(cluster.split("-")[1:3])
			additionalInosines_sub = {k:v for k, v in additionalInosines.items() if k == isodecoder and tRNA_dict[cluster]['species'] in v['species']}
			if additionalInosines_sub:
				tRNA_dict[cluster]['InosinePos'] = list(set(tRNA_dict[cluster]['InosinePos'] + additionalInosines_sub[isodecoder]['InosinePos']))
				Inosine_lists[cluster] = list(set(Inosine_lists[cluster] + additionalInosines_sub[isodecoder]['InosinePos']))

			total_inosines += len(Inosine_lists[cluster])

		# edit ref seqs A to G at inosine positions
		for cluster in Inosine_lists:
			for pos in Inosine_lists[cluster]:
				final_centroids[cluster].seq = final_centroids[cluster].seq[0:pos] + "G" + final_centroids[cluster].seq[pos+1:]

		# rewrite edited cluster transcripts
		with open(str(out_dir + experiment_name + '_clusterTranscripts.fa'), "w") as clusterTranscripts:
			SeqIO.write(final_centroids.values(), clusterTranscripts, "fasta")

		if total_snps == 0:
			snp_tolerance = False		

		log.info("{:,} modifications written to SNP index".format(total_snps))
		log.info("{:,} A to G replacements in reference sequences for inosine modifications".format(total_inosines))		
	
	# write outputs for indexing 
	with open(out_dir + experiment_name + "_modificationSNPs.txt", "w") as snp_file:
		for item in snp_records:
			snp_file.write('{}\n'.format(item))
	
	shutil.rmtree(temp_dir)
	
	# Return coverage_bed (either tRNAbed or clusterbed depending on --cluster) for coverage calculation method
	return(coverage_bed, snp_tolerance, mismatch_dict, mod_lists, Inosine_lists, tRNA_dict)

def newModsParser(out_dir, experiment_name, new_mods_list, new_Inosines, mod_lists, Inosine_lists, tRNA_dict, cluster):
# Parses new mods (from remap) into mod_lists, rewrites SNP index

	log.info("\n+------------------+ \
		\n| Parsing new mods |\
		\n+------------------+")	

	new_snps = 0
	new_inosines = 0

	# add new predicted inosines to inosine list and tRNA_dict
	for l in new_Inosines:
		for cluster, inosines in l.items():
			tRNA_dict[cluster]['InosinePos'] = list(set(tRNA_dict[cluster]['InosinePos'] + l[cluster]))
			old_inosines = len(Inosine_lists[cluster])
			Inosine_lists[cluster] = list(set(Inosine_lists[cluster] + l[cluster]))
			new_inosines += len(Inosine_lists[cluster]) - old_inosines

	log.info("{} new predicted position 34 inosines".format(new_inosines))

	# add new predicted mods to mods_list and tRNA_dict
	for l in new_mods_list:
		for cluster, mods in l.items():
			tRNA_dict[cluster]['modified'] = list(set(tRNA_dict[cluster]['modified'] + l[cluster]))
			old_mods = len(mod_lists[cluster])
			mod_lists[cluster] = list(set(mod_lists[cluster] + l[cluster]))
			new_snps += len(mod_lists[cluster]) - old_mods

	log.info("{} new predicted modifications".format(new_snps))

	# write file of new predicted mods
	predictedMods = open(out_dir + "mods/predictedMods.csv", "w")
	predictedMods.write("cluster\tpos\tidentity\tmisinc\n")
	with open(out_dir + "mods/predictedModstemp.csv", "r") as predictedTemp:
		for line in predictedTemp:
			cluster, pos, misinc = line.split("\t")
			identity = str(tRNA_dict[cluster]['sequence'][int(pos)])
			predictedMods.write(cluster + "\t" + str(pos) + "\t" + identity + "\t" + str(misinc) + "\n")
	os.remove(out_dir + "mods/predictedModstemp.csv")

	# rewrite SNP index
	total_snps = 0
	with open(out_dir + experiment_name + "_modificationSNPs.txt", "w") as snp_file:
		for cluster in mod_lists:
			for (index, pos) in enumerate(mod_lists[cluster]):
				snp_file.write(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "N\n")
			total_snps += len(mod_lists[cluster])
		for cluster in Inosine_lists:
		 	for (index, pos) in enumerate(Inosine_lists[cluster]):
		 		snp_file.write(">" + cluster + "_snp" + str(index) + "  " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "G\n")
		 	total_snps += len(Inosine_lists[cluster])

	# read in reference transcripts for inosine editing - use cluster to determine correct file
	if cluster:
		tRNA_ref = out_dir + experiment_name + '_clusterTranscripts.fa'
	else:
		tRNA_ref = out_dir + experiment_name + '_tRNATranscripts.fa'
	
	tRNA_seqs = SeqIO.to_dict(SeqIO.parse(tRNA_ref, 'fasta'))	

	# edit A to G for updated inosines list
	# for cluster in Inosine_lists:
	# 	for pos in Inosine_lists[cluster]:
	# 		tRNA_seqs[cluster].seq = tRNA_seqs[cluster].seq[0:pos] + "G" + tRNA_seqs[cluster].seq[pos+1:]

	# rewrite tRNA transcript reference
	with open(tRNA_ref, "w") as transcript_fasta:
		SeqIO.write(tRNA_seqs.values(), transcript_fasta, "fasta")
		
	log.info("{:,} modifications written to SNP index".format(total_snps))	

def additionalModsParser(input_species, out_dir):
# Reads in manual addition of modifcations in /data/additionalMods.txt

	mods ='/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/data/additionalMods.txt'
	mods = open(mods, 'r')
	additionalMods = defaultdict(lambda: defaultdict(list))

	# build dict of additional mods
	for line in mods:
		line = line.strip()
		species, tRNA, mods = line.split("\t")
		if species in input_species:
			additionalMods[tRNA]['mods'] = mods.split(";")
			additionalMods[tRNA]['species'] = species

	# initialise dictionaries of structure (with and without gapped numbering) and anticodon positions to define canonical mod sites
	tRNA_struct = ssAlign.tRNAclassifier(out_dir)
	tRNA_struct_nogap = ssAlign.tRNAclassifier_nogaps()
	cons_anticodon = ssAlign.getAnticodon()

	# dictionary storing additional mods per tRNA cluster with corrected positions
	additionalMods_parse = defaultdict(lambda: defaultdict(list))
	additionalInosines = defaultdict(lambda: defaultdict(list))

	# for each additional modification in dictionary, define mod site based on conserved location relative to structural features
	for isodecoder, data in additionalMods.items():
		struct = [value for key, value in tRNA_struct_nogap.items() if isodecoder in key and 'nmt' not in key]
		if struct: # test if struct is not empty, i.e. if isodecoder from additional mods exists in tRNA dictionary
			struct = struct[0] # assume same structure for all members of an isodecoder/cluster and therefore just use first one
		else:
			continue 
		cluster = [key for key, value in tRNA_struct_nogap.items() if isodecoder in key and 'nmt' not in key]
		if cluster:
			cluster = cluster[0]
		else:
			continue
		anticodon = ssAlign.clusterAnticodon(cons_anticodon, cluster)
		for mod in data['mods']:
			if mod == "m1A58":
				section = [pos for pos, value in struct.items() if value == "T stem-loop"]
				mod_site = section[-8]
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == "m2,2G26":
				section	= [pos for pos, value in struct.items() if value == 'bulge2']
				mod_site = section[0]
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == "m1G9":
				section = [pos for pos, value in struct.items() if value == 'bulge1']
				mod_site = section[-1]
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == "m1G37":
				mod_site = max(anticodon) + 1
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == 'm3C32':
				mod_site = min(anticodon) - 2
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == 'm3C47':
				section = [pos for pos, value in struct.items() if value == 'Variable loop']
				mod_site = section[-2]
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
			if mod == 'I34':
				mod_site = min(anticodon)
				additionalMods_parse[isodecoder]['mods'].append(mod_site)
				additionalMods_parse[isodecoder]['species'].append(data['species'])
				additionalInosines[isodecoder]['InosinePos'].append(mod_site)
				additionalInosines[isodecoder]['species'].append(data['species'])

			#if mod == 'm3C20':
			
	return(additionalMods_parse, additionalInosines)


def generateGSNAPIndices(experiment_name, out_dir, map_round, snp_tolerance = False, cluster = False):
# Builds genome and snp index files required by GSNAP

	if map_round == 1:
		log.info("\n+--------------------------+ \
		 \n| Generating GSNAP indices |\
		 \n+--------------------------+")

	else:
		log.info("\n+----------------------------+ \
		\n| Regenerating GSNAP indices |\
		\n+----------------------------+")

	genome_index_path = out_dir + experiment_name + "_tRNAgenome"
	genome_index_name = genome_index_path.split("/")[-1]
	
	try:
		os.mkdir(genome_index_path)
	except FileExistsError:
		log.warning("Genome index folder found! Rebuilding index anyway...")
	
	if cluster:
		genome_file = out_dir + experiment_name + "_clusterTranscripts.fa"
	else:
		genome_file = out_dir + experiment_name + "_tRNATranscripts.fa"

	index_cmd = "gmap_build -q 1 -D " + out_dir + " -d " + experiment_name + "_tRNAgenome " + genome_file + \
				" &> " + out_dir + "genomeindex.log"
	subprocess.call(index_cmd, shell = True) 
	log.info("Genome indices done...")

	snp_index_path = out_dir + experiment_name + "snp_index"

	if snp_tolerance:

		try:
			os.mkdir(snp_index_path)
		except FileExistsError:
			log.warning("SNP index folder found! Rebuilding index anyway...")

		snp_file = out_dir + experiment_name + "_modificationSNPs.txt"
		snp_index_name = snp_file.split("/")[-1]. split(".txt")[0]
		index_cmd = "cat " + snp_file + " | iit_store -o " + snp_index_path + "/" + snp_index_name + " &> " + out_dir + "snpindex.log"
		subprocess.call(index_cmd, shell = True)
		index_cmd = "snpindex -q 1 -D " + genome_index_path + " -d " + experiment_name + "_tRNAgenome -V " + snp_index_path + \
					" -v " + experiment_name + "_modificationSNPs " + snp_index_path + "/" + experiment_name + \
					"_modificationSNPs.iit &>> " + out_dir + "snpindex.log"
		subprocess.call(index_cmd, shell = True)
		log.info("SNP indices done...")
		return(genome_index_path, genome_index_name, snp_index_path, snp_index_name)

	elif not snp_tolerance:
		log.warning("SNP-tolerant alignment turned off or no modifications found for input tRNAs: SNP indices not built...\n")
		snp_index_name = ""
		return(genome_index_path, genome_index_name, snp_index_path, snp_index_name)

def generateSNPIndex(experiment_name, out_dir, snp_tolerance = False):
# Build SNP index only

	log.info("\n+------------------------------+ \
		 \n| Regenerating GSNAP SNP index |\
		 \n+------------------------------+")

	snp_index_path = out_dir + experiment_name + "snp_index"
	genome_index_path = out_dir + experiment_name + "_tRNAgenome"

	if snp_tolerance:

		try:
			os.mkdir(snp_index_path)
		except FileExistsError:
			log.warning("Rewriting SNP index...")

		snp_file = out_dir + experiment_name + "_modificationSNPs.txt"
		snp_index_name = snp_file.split("/")[-1]. split(".txt")[0]
		index_cmd = "cat " + snp_file + " | iit_store -o " + snp_index_path + "/" + snp_index_name + " &> " + out_dir + "snpindex.log"
		subprocess.call(index_cmd, shell = True)
		index_cmd = "snpindex -q 1 -D " + genome_index_path + " -d " + experiment_name + "_tRNAgenome -V " + snp_index_path + \
					" -v " + experiment_name + "_modificationSNPs " + snp_index_path + "/" + experiment_name + \
					"_modificationSNPs.iit &>> " + out_dir + "snpindex.log"
		subprocess.call(index_cmd, shell = True)
		log.info("SNP indices done...")

	elif not snp_tolerance:
		log.warning("SNP-tolerant alignment turned off or no modifications found for input tRNAs: SNP indices not built...\n")
		snp_index_name = ""

def modificationParser(modifications_table):
	# Read in modifications and build dictionary

		mods = open(modifications_table, 'r', encoding='utf-8')
		modifications = {}
		for line in mods:
			if not line.startswith("#"):
				name, abbr, ref, mod = line.split('\t')
				# replace unknown modifications with reference of N
				if not ref or ref.isspace():
					ref = 'N'
				if mod and not mod.isspace():
					modifications[mod.strip()] = {'name':name.strip(), 'abbr':abbr.strip(), 'ref':ref.strip()}
		return(modifications)

def getUnmodSeq (seq, modification_table):
# Change modified bases into standard ACGT in input sequence

	new_seq = []
	for char in seq:
		# for insertions ('_') make reference N - this is not described in the modifications table
		if char == '_':
			char = 'N'
		else:
			char = modification_table[char]['ref']
			# Change queuosine to G (reference is preQ0base in modification file)
			if char == 'preQ0base':
				char = 'G'

		new_seq.append(char)

	new_seq = ''.join(new_seq)
	new_seq = new_seq.replace('U','T')
	return(new_seq)

def initIntronDict(tRNAscan_out):
# Build dictionary of intron locations

	Intron_dict = {}
	tRNAscan = open(tRNAscan_out, 'r') 
	intron_count = 0
	for line in tRNAscan:
		if line.startswith("chr"):
			tRNA_ID = line.split()[0] + ".trna" + line.split()[1]
			tRNA_start = int(line.split()[2])
			intron_start = int(line.split()[6])
			intron_stop = int(line.split()[7])
			# if inton boundaries are not 0, i.e. there is an intron then add to dict
			if (intron_start > 0) & (intron_stop > 0):
				if tRNA_start > intron_start: # tRNA is on reverse strand
					intron_count += 1
					intron_start = tRNA_start - intron_start
					intron_stop = tRNA_start - intron_stop + 1 # needed for python 0 indexing and correct slicing of intron
				else: # tRNA is on forward strand
					intron_count += 1
					intron_start -= tRNA_start
					intron_stop -= tRNA_start
					intron_stop += 1 # python 0 indexing

				Intron_dict[tRNA_ID] = {}
				Intron_dict[tRNA_ID]['intron_start'] = intron_start
				Intron_dict[tRNA_ID]['intron_stop'] = intron_stop

	log.info("{} introns registered...".format(intron_count))
	return(Intron_dict)


def intronRemover (Intron_dict, seqIO_dict, seqIO_record, posttrans_mod_off):
# Use Intron_dict to find and remove introns plus add CCA and 5' G for His (if eukaryotic)

	# Find a match, slice intron and add G and CCA
	ID = re.search("tRNAscan-SE ID: (.*?)\).|\((chr.*?)-",seqIO_dict[seqIO_record].description).groups()
	ID = list(filter(None, ID))[0]
	if ID in Intron_dict:
		seq = str(seqIO_dict[seqIO_record].seq[:Intron_dict[ID]['intron_start']] + seqIO_dict[seqIO_record].seq[Intron_dict[ID]['intron_stop']:])
	else:
		seq = str(seqIO_dict[seqIO_record].seq)
	if 'His' in seqIO_record and posttrans_mod_off == False:
		seq = 'G' + seq + 'CCA'
	elif posttrans_mod_off == False:
		seq = seq + 'CCA'

	return(seq)

def tidyFiles (out_dir, cca):
	
	os.mkdir(out_dir + "annotation/")
	os.mkdir(out_dir + "align/")
	os.mkdir(out_dir + "indices/")
	os.mkdir(out_dir + "cov/")
	os.mkdir(out_dir + "counts/")

	files = os.listdir(out_dir)

	for file in files:
		full_file = out_dir + file
		if (file.endswith("bed") or file.endswith("stk") or file.endswith("gff") or file.endswith("fa") or "cm.log" in file or "clusterInfo" in file or "isoacceptorInfo" in file or "modificationSNPs" in file):
			shutil.move(full_file, out_dir + "annotation")
		if (file.endswith("tRNAgenome") or file.endswith("index") or "index.log" in file):
			shutil.move(full_file, out_dir + "indices")
		if (file.endswith("bam") or "align.log" in file or file == "mapping_stats.txt" or "alignstats.pdf" in file):
			shutil.move(full_file, out_dir + "align")
		if ("cov" in file):
			shutil.move(full_file, out_dir + "cov")
		if ("counts".upper() in file.upper()):
			shutil.move(full_file, out_dir + "counts")
