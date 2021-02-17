#!/usr/bin/env python3

##################################################################################
# Utilities for tRNA modification parsing, transcript building, and SNP indexing #
##################################################################################

from __future__ import absolute_import
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import re, copy, sys, os, shutil, subprocess, logging, glob
from pathlib import Path
from collections import defaultdict
import pandas as pd
import requests
from requests.models import HTTPError
from .ssAlign import aligntRNA, extraCCA, tRNAclassifier, tRNAclassifier_nogaps, getAnticodon, clusterAnticodon

log = logging.getLogger(__name__)

# This function is used to create nested defaultdicts like that needed for tRNA_dict
# Note this can be done with lambda functions (e.g. tRNA_dict = defaultdict(lambda: defaultdict()))
# But lambda functions cannot be pickled, and pickling is required for parallelization with multiprocessing. tRNA_dict is passed to such a multiprocessing pool later
def dd():
	return(defaultdict())

# same as above but specifically for list type defaultdicts (see insert_dict below)
def dd_list():
	return(defaultdict(list))

def tRNAparser (gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, posttrans_mod_off, double_cca, pretrnas, local_mod):
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
			if not pretrnas:
				tRNAseq = intronRemover(Intron_dict, temp_dict, seq, posttrans_mod_off, double_cca)
			else:
				tRNAseq = str(temp_dict[seq].seq)
			tRNA_dict[seq]['sequence'] = tRNAseq
			tRNA_dict[seq]['species'] = ' '.join(seq.split('_')[0:2])
			tRNA_dict[seq]['type'] = "cytosolic"

	# add mitochondrial tRNAs if given
	if mitotRNAs:
		temp_dict = SeqIO.to_dict(SeqIO.parse(mitotRNAs,"fasta"))
		mito_count = defaultdict(int)
		# read each mito tRNA, edit sequence header to match nuclear genes as above and add to tRNA_dict
		for seq in temp_dict:
			seq_parts = seq.split("|")
			anticodon = seq_parts[4]
			amino = re.search("[a-zA-z]+", seq_parts[3]).group(0)
			mito_count[anticodon] += 1
			new_seq = seq_parts[1] + "_mito_tRNA-" + amino + "-" + seq_parts[4] + "-" + str(mito_count[anticodon]) + "-1"
			tRNAseq = str(temp_dict[seq].seq) + "CCA" if not double_cca else str(temp_dict[seq].seq) + "CCACCA"
			tRNA_dict[new_seq]['sequence'] = tRNAseq
			tRNA_dict[new_seq]['type'] = 'mitochondrial'
			tRNA_dict[new_seq]['species'] = ' '.join(seq.split('_')[0:2])

		num_cytosilic = len([k for k in tRNA_dict.keys() if tRNA_dict[k]['type'] == "cytosolic"])
		num_mito = len([k for k in tRNA_dict.keys() if tRNA_dict[k]['type'] == "mitochondrial"])

		log.info("{} cytosolic and {} mitochondrial tRNA sequences imported".format(num_cytosilic, num_mito))

	# Read in and parse modomics file to contain similar headers to tRNA_dict
	# Save in new dict

	log.info("Processing modomics database...")
	modomics_file, fetch = getModomics(local_mod)
	modomics_dict, perSpecies_count = processModomics(modomics_file, fetch, species, modifications)

	for s in species:
		log.info('Number of Modomics entries for {}: {}'.format(s, perSpecies_count[s]))

	return(tRNA_dict,modomics_dict, species)

def processModomics(modomics_file, fetch, species, modifications):

	modomics_dict = dict()
	perSpecies_count = defaultdict(int)

	# build modomics_dict from JSON data through API
	if fetch:
		log.info("Parsing Modomics JSON data...")
		for data in modomics_file.values():
			sameIDcount = 0

			mod_species = data['organism']
			if not mod_species in species:
				continue
			else:
				perSpecies_count[mod_species] += 1
				anticodon = data['anticodon']
				new_anticodon = getUnmodSeq(anticodon, modifications)
				if "N" in new_anticodon:
					continue
				
				# Check amino acid name in modomics - set to iMet if equal to Ini to match gtRNAdb
				amino = data['subtype']
				if amino == 'Ini':
					amino = 'iMet'

				curr_id = str(mod_species.replace(' ','_') + '_tRNA-' + amino + '-' + new_anticodon)
				#Unique names for duplicates
				while curr_id in modomics_dict:
					sameIDcount += 1
					curr_id = "-".join(curr_id.split('-')[0:3]) + '-' + str(sameIDcount)

				tRNA_type = data['organellum']
				tRNA_type = "cytosolic" if re.search("cytosol",tRNA_type) else tRNA_type
				
				sequence = data['seq'].replace('U','T').replace('-','')
				unmod_sequence = getUnmodSeq(sequence, modifications)
				# Return list of modified nucl and inosines indices and add to modomics_dict
				# add unmodified seq to modomics_dict by lookup to modifications
				Mods = ['"', 'K', 'R', "'", 'O', 'Y', 'W', '⊆', 'X', '*', '[']
				modPos = [i for i, x in enumerate(sequence) if x in Mods]
				inosinePos = [i for i, x in enumerate(sequence) if x == 'I']
				
				modomics_dict[curr_id] = {'sequence':sequence,'type':tRNA_type, 'anticodon':new_anticodon, 'modified':modPos, 'unmod_sequence':unmod_sequence, 'InosinePos':inosinePos}

	# build modomics_dict if data was not fetched from API - i.e. using local txt version
	elif not fetch:
		log.info("Parsing local Modomics data...")
		for line in modomics_file:
			line = line.strip()
			sameIDcount = 0

			if line.startswith('>'):
				line = line.replace(' | ','|')
				# Check if species matches those from gtRNAdb input, otherwise skip entry
				mod_species = line.split('|')[3]
				if not mod_species in species:
					continue
				else:
					# Add to perSpecies_count to tally total modomics entries per species of interest
					perSpecies_count[mod_species] += 1
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
					while curr_id in modomics_dict:
						sameIDcount += 1
						curr_id = "-".join(curr_id.split('-')[0:3]) + '-' + str(sameIDcount)

					tRNA_type = str(line.split('|')[4])
					tRNA_type = "cytosolic" if re.search("cytosol",tRNA_type) else tRNA_type
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
					Mods = ['"', 'K', 'R', "'", 'O', 'Y', 'W', '⊆', 'X', '*', '[']
					modPos = [i for i, x in enumerate(modomics_dict[curr_id]['sequence']) if x in Mods]
					inosinePos = [i for i, x in enumerate(modomics_dict[curr_id]['sequence']) if x == 'I']
					modomics_dict[curr_id]['modified'] = modPos
					modomics_dict[curr_id]['unmod_sequence'] = unmod_sequence
					modomics_dict[curr_id]["InosinePos"] = inosinePos
		
	return(modomics_dict, perSpecies_count)

def getModomics(local_mod):
	# Get full Modomics modified tRNA data from web
	fetch = False
	if not local_mod:
		try:
			response = requests.get("http://www.genesilico.pl/modomics/api/sequences?&RNAtype=tRNA")
			response.raise_for_status()
			modomics = response.json()
			fetch = True
			log.info("Modomics retrieved...")
		except HTTPError as http_err:
			log.error("Unable to connect to Modomics database! HTTP error: {}. Check status of Modomics webpage. Using local Modomics files...".format(http_err))
			modomics_path = os.path.dirname(os.path.realpath(__file__)) + '/data/modomics'
			modomics = open(modomics_path, "r+", encoding = "utf-8")
		except Exception as err:
			log.error("Error in connecting to Modomics: {}. Using local Modomics files...".format(err))
			modomics_path = os.path.dirname(os.path.realpath(__file__)) + '/data/modomics'
			modomics = open(modomics_path, "r+", encoding = "utf-8")
	else:
		log.warning("Retrieval of Modomics database disabled. Using local files instead...")
		modomics_path = os.path.dirname(os.path.realpath(__file__)) + '/data/modomics'
		modomics = open(modomics_path, "r+", encoding = "utf-8")

	return modomics, fetch

def modsToSNPIndex(gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, experiment_name, out_dir, double_cca, snp_tolerance = False, cluster = False, cluster_id = 0.95, posttrans_mod_off = False, pretrnas = False, local_mod = False):
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
	tRNA_dict, modomics_dict, species = tRNAparser(gtRNAdb, tRNAscan_out, mitotRNAs, modifications_table, posttrans_mod_off, double_cca, pretrnas, local_mod)
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

	aligntRNA(tempSeqs.name, out_dir)
	extra_cca = extraCCA()

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
		tRNA_dict[seq]['anticodon'] = anticodon = re.search('.*tR(NA|X)-.*?-(.*?)-', seq).group(2)
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
					if (hsp.bits > maxbit) and (hsp.align_length / alignment.length == 1) and (hsp.identities / alignment.length == 1):
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
		seq_records[str(seq)] = SeqRecord(Seq(tRNA_dict[seq]['sequence'].upper()), id = str(seq))

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
		aligntRNA(temptRNATranscripts.name, out_dir)
		additionalMods, additionalInosines = additionalModsParser(species, out_dir)
		# add additional SNPs from extra file to list of modified positions, and ensure non-redundancy with set()
		# index SNPs
		# Format for SNP index (space separated):
		# >snpID chromosomeName:position(1-based) RefMod
		# e.g. >rs111 Homo_sapiens_nmt_tRNA-Leu-TAA-1-1_exp0:29 GN
		for seq in mod_lists:
			additionalMods_sub = {k:v for k, v in additionalMods.items() if k == seq and tRNA_dict[seq]['species'] in v['species']}
			if additionalMods_sub:
				tRNA_dict[seq]['modified'] = list(set(tRNA_dict[seq]['modified'] + additionalMods_sub[seq]['mods']))
				mod_lists[seq] = list(set(mod_lists[seq] + additionalMods_sub[seq]['mods']))

			total_snps += len(mod_lists[seq])

			# Build snp_records as before but with cluster names and non-redundant sets of modifications
			# Position is 1-based for iit_store i.e. pos + 1
			for (index, pos) in enumerate(mod_lists[seq]):
				#if "Gln-TTG" in seq and pos + 1 == 34:
				#	snp_records.append(">" + seq + "_snp" + str(index) + " " + seq + ":" + str(pos + 1) + " " + tRNA_dict[seq]['sequence'][pos].upper() + "C")
				#else:
				snp_records.append(">" + seq + "_snp" + str(index) + " " + seq + ":" + str(pos + 1) + " " + tRNA_dict[seq]['sequence'][pos].upper() + "N")

		for seq in Inosine_lists:
			additionalInosines_sub = {k:v for k, v in additionalInosines.items() if k == seq and tRNA_dict[seq]['species'] in v['species']}
			if additionalInosines_sub:
				tRNA_dict[seq]['InosinePos'] = list(set(tRNA_dict[seq]['InosinePos'] + additionalInosines_sub[seq]['InosinePos']))
				Inosine_lists[seq] = list(set(Inosine_lists[seq] + additionalInosines_sub[seq]['InosinePos']))

			total_inosines += len(Inosine_lists[seq])
			total_snps += len(Inosine_lists[seq])

			for (index, pos) in enumerate(Inosine_lists[seq]):
				# Add inosines to snp index (in addition to changing ref seqence - see below)
				# Ensure only "A" SNPs are tolerated. That is, a G in the reference allows inosines while an A in the snp index allows unmodified reads
				snp_records.append(">" + seq + "_snp" + str(index) + " " + seq + ":" + str(pos + 1) + " " + "GA")


		Inosine_clusters = [cluster for cluster, inosines in Inosine_lists.items() if len(inosines) > 0]

		# edit ref seqs A to G at inosine positions only if snp_tolerance is enabled
		if snp_tolerance:
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
		insert_dict = defaultdict(dd_list)
		del_dict = defaultdict(dd_list)
		mod_lists = {tRNA:list() for tRNA in tRNA_dict.keys()}
		Inosine_lists = {tRNA:data['InosinePos'] for tRNA, data in tRNA_dict.items()}
		cluster_perPos_mismatchMembers = defaultdict(dd_list)
		cluster_dict = defaultdict(list)

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
			cluster_cmd = ["usearch", "-cluster_fast", temp_dir + anticodon + "_allseqs.fa", "-id", str(cluster_id), "-sizeout" ,"-centroids", temp_dir + anticodon + "_centroids.fa", "-uc", temp_dir + anticodon + "_clusters.uc"]
			#cluster_cmd = ["usearch", "-cluster_smallmem", temp_dir + anticodon + "_allseqs.fa", "-id", str(cluster_id), "--sortedby", "other" ,"-sizeout" ,"-centroids", temp_dir + anticodon + "_centroids.fa", "-uc", temp_dir + anticodon + "_clusters.uc"]			#cluster_cmd = "usearch -cluster_fast " + temp_dir + anticodon + "_allseqs.fa -sort length -id " + str(cluster_id) + " -centroids " + temp_dir + anticodon + "_centroids.fa -uc " + temp_dir + anticodon + "_clusters.uc &> /dev/null"
			subprocess.check_call(cluster_cmd, stdout = subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			# sort clusters by size (i.e. number of members in cluster)
			sort_cmd = ["usearch", "-sortbysize", temp_dir + anticodon + "_centroids.fa", "-fastaout", temp_dir + anticodon + "_centroids_sort.fa"]
			subprocess.check_call(sort_cmd, stdout = subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			# recluster based on sorted by size clusters
			final_cluster_cmd = ["usearch", "-cluster_smallmem", temp_dir + anticodon + "_centroids_sort.fa", "-sortedby", "size", "-id", str(cluster_id), "-centroids", temp_dir + anticodon + "_centroidsFinal.fa"]
			subprocess.check_call(final_cluster_cmd, stdout = subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		# combine centroids files into one file
		for filename in glob.glob(temp_dir + "*_centroidsFinal.fa"):
			with open(filename, "r") as fileh:
				with open(temp_dir + "all_centroids.fa", "a") as outh:
					outh.write(fileh.read())
		centroids = SeqIO.parse(temp_dir + "all_centroids.fa", "fasta")
		for centroid in centroids:
			centroid.id = centroid.id.split(";")[0]
			final_centroids[centroid.id] = SeqRecord(Seq(str(centroid.seq).upper()), id = centroid.id) 

		# read cluster files, get nonredudant set of mod positions of all members of a cluster, create snp_records for writing SNP index
		cluster_pathlist = Path(temp_dir).glob("**/*_clusters.uc")
		mod_lists = dict() # stores non-redundant sets of mismatches and mod positions for clusters
		Inosine_lists = dict() # stores positions of inosines for clusters
		snp_records = list()
		cluster_dict = defaultdict(list) # info about clusters and members
		mismatch_dict = defaultdict(list) # dictionary of mismatches only (not mod positions - required for exclusion from misincorporation analysis in mmQuant)
		insert_dict = defaultdict(dd_list) # dictionary of insertions in cluster parents - these are not recorded as mismatches but are needed in order to split clusters into isodecoders
		del_dict = defaultdict(dd_list) # dictionary of deletions in cluster parents - these are not recorded as mismatches but are needed in order to split clusters into isodecoders
		cluster_perPos_mismatchMembers = defaultdict(dd_list) # for each cluster, list of members that mismatch at each position corresonding to mismatch dict - used for splitting clusters into isodecoders (splitClusters.py)
		cluster_num = 0
		total_snps = 0
		total_inosines = 0
		clusterbed = open(out_dir + experiment_name + "_clusters.bed","w")
		coverage_bed = clusterbed.name
		clustergff = open(out_dir + experiment_name + "_tRNA.gff","w")
		for path in cluster_pathlist:
			with open(str(path),"r") as cluster_file:
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
						cluster_dict[cluster_name].append(cluster_name)
				
					# Handle members of clusters
					elif line.split("\t")[0] == "H":
						member_name = line.split("\t")[8].split(";")[0]
						cluster_name = line.split("\t")[9].split(";")[0]
						compr_aln = line.split("\t")[7]
						# if member of cluster is 100% identical (i.e. "=" in 8th column of cluster file)
						if compr_aln == "=":
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + tRNA_dict[member_name]["modified"]))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + tRNA_dict[member_name]["InosinePos"]))
							cluster_dict[cluster_name].append(member_name)
						
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
								adjust_pos_len = 0
								for insert in insertion_pos:
									if insert < delete:
										adjust_pos_len -= 1
								new_delete = delete + adjust_pos_del + adjust_pos_len
								member_seq = member_seq[ :new_delete] + member_seq[new_delete+1: ]
								del_dict[cluster_name][new_delete].append(member_name)
								adjust_pos_del -= 1

							for index, insert in enumerate(insertion_pos):
								adjust_pos_len = 0
								for delete in deletion_pos:
									if delete < insert:
										adjust_pos_len -= 1
								#if index != 0:
								#	if insert == insertion_pos[index-1] + 1:
								#		adjust_pos_ins -= 1
								new_insert = insert + adjust_pos_len + adjust_pos_ins
								member_seq = member_seq[ :new_insert] + cluster_seq[insert] + member_seq[new_insert: ]
								#insert_dict[cluster_name][new_insert-1].append(member_name)
								insert_dict[cluster_name][new_insert].append(member_name)
								#cluster_seq = cluster_seq[ :new_insert] + cluster_seq[new_insert+1: ]

							mismatches = [i for i in range(len(member_seq)) if member_seq[i].upper() != cluster_seq[i].upper()]
							mismatch_dict[cluster_name] = list(set(mismatch_dict[cluster_name] + mismatches))
							for mismatch in mismatches:
								cluster_perPos_mismatchMembers[cluster_name][mismatch].append(member_name)
							member_mods = list(set(tRNA_dict[member_name]["modified"] + mismatches))
							member_Inosines = tRNA_dict[member_name]["InosinePos"]
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + member_mods))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + member_Inosines))
							cluster_dict[cluster_name].append(member_name)

						# handle members that are not exact sequence matches but have no indels either
						# find mismatches and build non-redundant set
						else:
							cluster_seq = tRNA_dict[cluster_name]["sequence"]
							member_seq = tRNA_dict[member_name]["sequence"]
							mismatches = [i for i in range(len(member_seq)) if member_seq[i].upper() != cluster_seq[i].upper()]
							mismatch_dict[cluster_name] = list(set(mismatch_dict[cluster_name] + mismatches))
							for mismatch in mismatches:
								cluster_perPos_mismatchMembers[cluster_name][mismatch].append(member_name)
							member_mods = list(set(tRNA_dict[member_name]["modified"] + mismatches))
							member_Inosines = tRNA_dict[member_name]["InosinePos"]
							mod_lists[cluster_name] = list(set(mod_lists[cluster_name] + member_mods))
							Inosine_lists[cluster_name] = list(set(Inosine_lists[cluster_name] + member_Inosines))
							cluster_dict[cluster_name].append(member_name)

		clusterbed.close()

		# Write cluster information to tsv and get number of unique sequences per cluster (i.e. isodecoders) for read count splitting
		with open(out_dir + experiment_name + "clusterInfo.txt","w") as clusterInfo, open(out_dir + experiment_name + "isoacceptorInfo.txt","w") as isoacceptorInfo:
			isoacceptor_dict = defaultdict(int)
			clusterInfo.write("tRNA\tcluster_num\tcluster_size\tparent\n")
			isoacceptorInfo.write("Isoacceptor\tsize\n")
			for key, value in cluster_dict.items():
				cluster_num = list(cluster_dict.keys()).index(key)
				for member in value:
					clusterInfo.write("{}\t{}\t{}\t{}\n".format(member, cluster_num, len(value), key))
					isoacceptor_group = '-'.join(member.split("-")[:-2])
					isoacceptor_dict[isoacceptor_group] += 1
			for key, value in isoacceptor_dict.items():
				isoacceptorInfo.write(key + "\t" + str(value) + "\n")
		with open(str(out_dir + experiment_name + '_clusterTranscripts.fa'), "w") as clusterTranscripts:
			SeqIO.write(final_centroids.values(), clusterTranscripts, "fasta")

		# generate Stockholm alignment file for cluster transcripts and process additional mods file
		aligntRNA(clusterTranscripts.name, out_dir)
		additionalMods, additionalInosines = additionalModsParser(species, out_dir)

		log.info("{} clusters created from {} tRNA sequences".format(len(cluster_dict),len(tRNA_dict)))

		# update mod_lists with additional mods and write SNP index
		for cluster in mod_lists:
			additionalMods_sub = {k:v for k, v in additionalMods.items() if k == cluster and tRNA_dict[cluster]['species'] in v['species']}
			if additionalMods_sub:
				tRNA_dict[cluster]['modified'] = list(set(tRNA_dict[cluster]['modified'] + additionalMods_sub[cluster]['mods']))
				mod_lists[cluster] = list(set(mod_lists[cluster] + additionalMods_sub[cluster]['mods']))

			total_snps += len(mod_lists[cluster])

			for (index, pos) in enumerate(mod_lists[cluster]):
				# Build snp_records as before but with cluster names and non-redundant sets of modifications
				# Position is 1-based for iit_store i.e. pos + 1
				#if "Gln-TTG" in cluster and pos + 1 == 34:
				#	snp_records.append(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "C")
				#else:
				snp_records.append(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "N")

		for cluster in Inosine_lists:
			additionalInosines_sub = {k:v for k, v in additionalInosines.items() if k == cluster and tRNA_dict[cluster]['species'] in v['species']}
			if additionalInosines_sub:
				tRNA_dict[cluster]['InosinePos'] = list(set(tRNA_dict[cluster]['InosinePos'] + additionalInosines_sub[cluster]['InosinePos']))
				Inosine_lists[cluster] = list(set(Inosine_lists[cluster] + additionalInosines_sub[cluster]['InosinePos']))

			total_inosines += len(Inosine_lists[cluster])
			total_snps += len(Inosine_lists[cluster])

			for (index, pos) in enumerate(Inosine_lists[cluster]):
				# Add inosines to snp index (in addition to changing ref seqence - see below)
				# Ensure only "A" SNPs are tolerated. That is, a G in the reference allows inosines while an A in the snp index allows unmodified reads
				snp_records.append(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + "GA")


		# edit ref seqs A to G at inosine positions if snp_tolerance is enabled
		if snp_tolerance:
			for cluster in Inosine_lists:
				for pos in Inosine_lists[cluster]:
					final_centroids[cluster].seq = final_centroids[cluster].seq[0:pos] + "G" + final_centroids[cluster].seq[pos+1:]

		Inosine_clusters = [cluster for cluster, inosines in Inosine_lists.items() if len(inosines) > 0]

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
	return(coverage_bed, snp_tolerance, mismatch_dict, insert_dict, del_dict, mod_lists, Inosine_lists, Inosine_clusters, tRNA_dict, cluster_dict, cluster_perPos_mismatchMembers)

def newModsParser(out_dir, experiment_name, new_mods_list, new_Inosines, mod_lists, Inosine_lists, tRNA_dict, clustering, remap, snp_tolerance):
# Parses new mods (from remap) into mod_lists, rewrites SNP index

	log.info("\n+------------------+ \
		\n| Parsing new mods |\
		\n+------------------+")	

	new_snps_num = 0
	new_inosines = 0
	newInosine_lists = defaultdict(list)

	# keep original inosine list (i.e. before addding new ones) to modify misincorporation proportion in mismatchTable
	Inosine_clusters = [cluster for cluster, inosines in Inosine_lists.items() if len(inosines) > 0]

	# add new predicted inosines to inosine list and tRNA_dict
	for l in new_Inosines:
		for cluster in l.keys():
			tRNA_dict[cluster]['InosinePos'] = list(set(tRNA_dict[cluster]['InosinePos'] + l[cluster]))
			old_inosines = len(Inosine_lists[cluster])
			Inosine_lists[cluster] = list(set(Inosine_lists[cluster] + l[cluster])) # total inosines per cluster
			newInosine_lists[cluster] = list(set(newInosine_lists[cluster] + l[cluster])) # only new inosines to be written to SNP index
			new_inosines += len(Inosine_lists[cluster]) - old_inosines

	log.info("{} new predicted position 34 inosines".format(new_inosines))

	# add new predicted mods to mods_list and tRNA_dict
	for l in new_mods_list:
		for cluster in l.keys():
			tRNA_dict[cluster]['modified'] = list(set(tRNA_dict[cluster]['modified'] + l[cluster]))
			old_mods = len(mod_lists[cluster])
			mod_lists[cluster] = list(set(mod_lists[cluster] + l[cluster]))
			new_snps_num += len(mod_lists[cluster]) - old_mods

	# update snp_tolerance if previously False - i.e. no Modomics data but new mods discovered after alignment
	if snp_tolerance == False and not new_snps_num == 0:
		snp_tolerance = True

	log.info("{} new predicted modifications".format(new_snps_num))

	if remap:
		# rewrite SNP index
		total_snps = 0
		with open(out_dir + experiment_name + "_modificationSNPs.txt", "w") as snp_file:
			for cluster in mod_lists:
				for (index, pos) in enumerate(mod_lists[cluster]):
					snp_file.write(">" + cluster + "_snp" + str(index) + " " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "N\n")
				total_snps += len(mod_lists[cluster])
			for cluster in newInosine_lists:
				for (index, pos) in enumerate(Inosine_lists[cluster]):
					snp_file.write(">" + cluster + "_snp" + str(index) + "  " + cluster + ":" + str(pos + 1) + " " + tRNA_dict[cluster]['sequence'][pos].upper() + "G\n")
				total_snps += len(Inosine_lists[cluster])

		# read in reference transcripts for inosine editing - use cluster to determine correct file
		if clustering:
			tRNA_ref = out_dir + experiment_name + '_clusterTranscripts.fa'
		else:
			tRNA_ref = out_dir + experiment_name + '_tRNATranscripts.fa'
		
		tRNA_seqs = SeqIO.to_dict(SeqIO.parse(tRNA_ref, 'fasta'))	

		# rewrite tRNA transcript reference
		with open(tRNA_ref, "w") as transcript_fasta:
			SeqIO.write(tRNA_seqs.values(), transcript_fasta, "fasta")
				
		log.info("{:,} modifications written to SNP index".format(total_snps))	

	return(Inosine_clusters, snp_tolerance, tRNA_dict, mod_lists)

def additionalModsParser(input_species, out_dir):
# Reads in manual addition of modifcations in /data/additionalMods.txt

	mods = os.path.dirname(os.path.realpath(__file__)) + '/data/additionalMods.txt'
	mods = open(mods, 'r')
	additionalMods = defaultdict(lambda: defaultdict(list))

	# build dict of additional mods
	for line in mods:
		line = line.strip()
		species, tRNA, mods = line.split("\t")
		if species in input_species:
			if "mito" in tRNA:
				tRNA = tRNA.replace("mito", "")
				additionalMods[tRNA]['type'] = "mitochondrial"
			else:
				additionalMods[tRNA]['type'] = "cytosolic"
			additionalMods[tRNA]['mods'] = mods.split(";")
			additionalMods[tRNA]['species'] = species

	# initialise dictionaries of structure (with and without gapped numbering) and anticodon positions to define canonical mod sites
	tRNA_struct = tRNAclassifier()[0]
	cons_pos_dict = tRNAclassifier()[3]
	tRNA_struct_nogap = tRNAclassifier_nogaps()
	cons_anticodon = getAnticodon()

	# dictionary storing additional mods per tRNA cluster with corrected positions
	additionalMods_parse = defaultdict(lambda: defaultdict(list))
	additionalInosines = defaultdict(lambda: defaultdict(list))

	# for each additional modification in dictionary, define mod site based on conserved location relative to structural features
	for isodecoder, data in additionalMods.items():
		if data['type'] == "mitochondrial":
			clusters = [key for key, value in tRNA_struct_nogap.items() if isodecoder in key and 'nmt' not in key and "mito" in key]
		else:
			clusters = [key for key, value in tRNA_struct_nogap.items() if isodecoder in key and 'nmt' not in key and "mito" not in key]
		for cluster in clusters:
			no_gap_struct = [value for key, value in tRNA_struct_nogap.items() if key == cluster and 'nmt' not in key]
			if not no_gap_struct: # test if struct is empty, i.e. if isodecoder from additional mods does not exist in tRNA dictionary
				continue		
			anticodon = clusterAnticodon(cons_anticodon, cluster)

			for mod in data['mods']:
				if not 'I' in mod:
					cons_pos = re.search('.*?[A|C|G|U]([0-9].*)', mod).group(1)
					mod_site = getModSite(cluster, cons_pos, cons_pos_dict, tRNA_struct, tRNA_struct_nogap)
					if not mod_site == 'NA':
						additionalMods_parse[cluster]['mods'].append(mod_site)
						additionalMods_parse[cluster]['species'].append(data['species'])

				if mod == 'I34':
					mod_site = min(anticodon)
					additionalInosines[cluster]['InosinePos'].append(mod_site)
					additionalInosines[cluster]['species'].append(data['species'])

	return(additionalMods_parse, additionalInosines)

def getModSite(cluster, cons_pos, cons_pos_dict, tRNA_struct, tRNA_struct_nogap):
# Gets mod site for individual clusters based on conserved position

	gapped_pos = list(cons_pos_dict)[list(cons_pos_dict.values()).index(cons_pos)]
	struct_element = tRNA_struct[cluster][gapped_pos]
	if not struct_element == 'gap':
		all_struct_element = {pos: element for pos, element in tRNA_struct[cluster].items() if element == struct_element}
		all_struct_element_list = [pos for pos in all_struct_element.keys()]
		all_struct_element_list = sorted(all_struct_element_list)
		index_struct_element = all_struct_element_list.index(gapped_pos)
		all_struct_element_nogap = {pos: element for pos, element in tRNA_struct_nogap[cluster].items() if element == struct_element}
		all_struct_element_list_nogap = [pos for pos in all_struct_element_nogap.keys()]
		all_struct_element_list_nogap = sorted(all_struct_element_list_nogap)
		mod_site = all_struct_element_list_nogap[index_struct_element]
	else:
		mod_site = 'NA'
	
	return(mod_site)

def generateGSNAPIndices(species, name, out_dir, map_round, snp_tolerance = False, cluster = False):
# Builds genome and snp index files required by GSNAP

	if map_round == 1:
		log.info("\n+--------------------------+ \
		 \n| Generating GSNAP indices |\
		 \n+--------------------------+")

	else:
		log.info("\n+----------------------------+ \
		\n| Regenerating GSNAP indices |\
		\n+----------------------------+")

	genome_index_path = out_dir + species + "_tRNAgenome"
	genome_index_name = genome_index_path.split("/")[-1]
	
	try:
		os.mkdir(genome_index_path)
	except FileExistsError:
		log.warning("Genome index folder found! Rebuilding index anyway...")
	
	if cluster:
		genome_file = out_dir + name + "_clusterTranscripts.fa"
	else:
		genome_file = out_dir + name + "_tRNATranscripts.fa"

	index_cmd = ["gmap_build", "-q", "1", "-D", out_dir, "-d", species + "_tRNAgenome", genome_file]
	subprocess.check_call(index_cmd, stderr = open(out_dir + "genomeindex.log", "w"), stdout = subprocess.DEVNULL) 
	log.info("Genome indices done...")

	snp_index_path = out_dir + species + "snp_index"

	if snp_tolerance:

		try:
			os.mkdir(snp_index_path)
		except FileExistsError:
			log.warning("SNP index folder found! Rebuilding index anyway...")

		snp_file = out_dir + name + "_modificationSNPs.txt"
		snp_index_name = snp_file.split("/")[-1]. split(".txt")[0]
		ps = subprocess.Popen(('cat', snp_file), stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)
		index_cmd = ["iit_store", "-o", snp_index_path + "/" + snp_index_name]
		subprocess.check_call(index_cmd, stdin = ps.stdout, stdout = open(out_dir + "snpindex.log", "w"), stderr = subprocess.DEVNULL)
		index_cmd = ["snpindex", "-q", "1", "-D", genome_index_path, "-d", species + "_tRNAgenome", "-V", snp_index_path, \
					"-v", name + "_modificationSNPs", snp_index_path + "/" + name + \
					"_modificationSNPs.iit"]
		subprocess.check_call(index_cmd, stderr = open(out_dir + "snpindex.log", "w"), stdout = subprocess.DEVNULL)
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
		ps = subprocess.Popen(('cat', snp_file), stdout = subprocess.PIPE)
		index_cmd = ["iit_store", "-o", snp_index_path + "/" + snp_index_name]
		subprocess.check_call(index_cmd, stdin = ps.stdout, stdout = open(out_dir + "snpindex.log", "w"))
		index_cmd = ["snpindex", "-q", "1", "-D", genome_index_path, "-d", experiment_name + "_tRNAgenome", "-V", snp_index_path, \
					"-v", experiment_name + "_modificationSNPs", snp_index_path + "/" + experiment_name + \
					"_modificationSNPs.iit"]
		subprocess.check_call(index_cmd, stderr = open(out_dir + "snpindex.log", "w"))
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

def getUnmodSeq(seq, modification_table):
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


def intronRemover (Intron_dict, seqIO_dict, seqIO_record, posttrans_mod_off, double_cca):
# Use Intron_dict to find and remove introns plus add CCA and 5' G for His (if eukaryotic)

	# Find a match, slice intron and add G and CCA
	ID = re.search("tRNAscan-SE ID: (.*?)\).|\((chr.*?)-",seqIO_dict[seqIO_record].description).groups()
	ID = list(filter(None, ID))[0]
	if ID in Intron_dict:
		seq = str(seqIO_dict[seqIO_record].seq[:Intron_dict[ID]['intron_start']] + seqIO_dict[seqIO_record].seq[Intron_dict[ID]['intron_stop']:])
	else:
		seq = str(seqIO_dict[seqIO_record].seq)
	if posttrans_mod_off == False:
		if double_cca:
			seq = seq + 'CCACCA'
		else:
			seq = seq + 'CCA'
		if 'His' in seqIO_record:
			seq = 'G' + seq

	return(seq)

def countReads(input_counts, out_dir, isodecoder_sizes, clustering, tRNA_dict, clusterInfo):

	# Counts per anticodon
	count_dict_anticodon = defaultdict(lambda: defaultdict(int))
	count_dict_isodecoder = defaultdict(lambda: defaultdict(int))

	with open(input_counts, "r") as counts_file:
		for line in counts_file:
			line = line.strip()
			if not line.startswith("#"):
				if line.startswith("isodecoder"):
					sample_list = [samples for samples in line.split("\t")[1:-3]]
				else:
					isodecoder = line.split("\t")[0]
					anticodon = "-".join(isodecoder.split("-")[:-1]) if not "chr" in isodecoder else "-".join(isodecoder.split("-")[:-2])
					col = 1
					for sample in sample_list:
						count_dict_anticodon[anticodon][sample] += float(line.split("\t")[col])
						if not clustering:
							count_dict_isodecoder[isodecoder][sample] = float(line.split("\t")[col])
						col += 1

	count_anticodon_pd = pd.DataFrame.from_dict(count_dict_anticodon, orient='index')
	count_anticodon_pd.index.name = 'Anticodon'
	count_anticodon_pd.to_csv(out_dir + 'Anticodon_counts.txt', sep = '\t')

	if not clustering:
		new_count_isodecoder = defaultdict(lambda: defaultdict(int))
		for isodecoder in isodecoder_sizes:
			sameSeq = [tRNAs for tRNAs in tRNA_dict.keys() if tRNA_dict[tRNAs]['sequence'] == tRNA_dict[isodecoder]['sequence']]
			for i in sameSeq:
				i = "-".join(i.split("-")[:-1])
				for lib in count_dict_isodecoder[isodecoder].keys():
					new_count_isodecoder[isodecoder][lib] += count_dict_isodecoder[i][lib]

		count_isodecoder_pd = pd.DataFrame.from_dict(new_count_isodecoder, orient='index')
		count_isodecoder_pd.index.name = 'isodecoder'
		count_isodecoder_pd['Single_isodecoder'] = "True"

		isodecoder_sizes_short = defaultdict()
		for iso, size in isodecoder_sizes.items():
			if not "chr" in iso:
				short = "-".join(iso.split("-")[:-1])
			else:
				short = iso
			isodecoder_sizes_short[short] = size
		for cluster in count_isodecoder_pd.index:
			count_isodecoder_pd.at[cluster, 'size'] = isodecoder_sizes_short[cluster]
		count_isodecoder_pd = count_isodecoder_pd.join(clusterInfo)
		count_isodecoder_pd.to_csv(out_dir + 'Isodecoder_counts.txt', sep = '\t')

	log.info("** Read counts per anticodon saved to " + out_dir + "counts/Anticodon_counts.txt **")

def tidyFiles (out_dir, cca):
	
	os.mkdir(out_dir + "annotation/")
	os.mkdir(out_dir + "align/")
	os.mkdir(out_dir + "indices/")
	os.mkdir(out_dir + "cov/")
	os.mkdir(out_dir + "counts/")

	files = os.listdir(out_dir)

	for file in files:
		full_file = out_dir + file
		if (file.endswith("bed") or file.endswith("stk") or file.endswith("gff") or file.endswith("fa") or "cm.log" in file or "clusterInfo" in file or "isoacceptorInfo" in file or "isodecoderInfo" in file or "modificationSNPs" in file):
			shutil.move(full_file, out_dir + "annotation")
		if (file.endswith("tRNAgenome") or file.endswith("index") or "index.log" in file):
			shutil.move(full_file, out_dir + "indices")
		if (file.endswith("bam") or file.endswith("bam.bai") or "align.log" in file or file == "mapping_stats.txt" or "alignstats.pdf" in file):
			shutil.move(full_file, out_dir + "align")
		if ("cov" in file):
			shutil.move(full_file, out_dir + "cov")
		if ("counts".upper() in file.upper()):
			shutil.move(full_file, out_dir + "counts")
