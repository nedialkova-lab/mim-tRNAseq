#! /usr/bin/env python3

#   +-------------+
#   | mim-tRNAseq | 
#   +-------------+

####################################
# Main backbone and wrapper script #
####################################
# 
# author: Drew Behrens
# contact: aberens@biochem.mpg.de
# github: https://github.com/drewjbeh/mim-tRNAseq

import tRNAtools, tRNAmap, getCoverage, mmQuant, CCAanalysis, ssAlign
import sys, os, subprocess, logging, datetime
import argparse
from pyfiglet import figlet_format

def restrictedFloat(x):
## Method for restricting cluster_id argument to float between 0 and 1
	try:
		x = float(x)
		if x < 0.0 or x > 1.0:
			raise argparse.ArgumentTypeError('{} not in range 0.0 - 1.0'.format(x))
		return x
	except ValueError:
		raise argparse.ArgumentTypeError('{} not a real number'.format(x))

def mimseq(trnas, trnaout, name, out, cluster, cluster_id, posttrans, control_cond, threads, max_multi, snp_tolerance, \
	keep_temp, mode, cca, min_cov, mismatches, remap, misinc_thresh, mito_trnas, sample_data):
	
# Main wrapper

	# Integrity check for output folder argument...
	try:
		os.mkdir(out)
	except FileExistsError:
		raise FileExistsError("Output folder already exists!")

	if not out.endswith("/"):
		out = out + "/"

	###########
	# Logging #
	###########

	now = datetime.datetime.now()
	logging.basicConfig(
		format="%(asctime)s [%(levelname)-5.5s] %(message)s",
		level=logging.INFO,
		handlers=[
			logging.FileHandler(out + "mim-tRNAseq_{}.log".format(now.strftime("%H-%M-%S"))),
			logging.StreamHandler(sys.stdout)
		])
	log = logging.getLogger(__name__)
	log.info("mim-tRNAseq run with command:")
	log.info(" ".join(sys.argv))


	########
	# main #
	########

	# Parse tRNA and modifications, generate SNP index
	modifications = os.path.dirname(os.path.realpath(__file__))
	modifications += "/modifications"
	coverage_bed, snp_tolerance, mismatch_dict, mod_lists, tRNA_dict = tRNAtools.modsToSNPIndex(trnas, trnaout, mito_trnas, modifications, name, out, snp_tolerance, cluster, cluster_id, posttrans)
	ssAlign.structureParser()
	# Generate GSNAP indices
	genome_index_path, genome_index_name, snp_index_path, snp_index_name = tRNAtools.generateGSNAPIndices(name, out, snp_tolerance, cluster)

	# Align
	map_round = 1 #first round of mapping
	bams_list, coverageData = tRNAmap.mainAlign(sample_data, name, genome_index_path, genome_index_name, \
		snp_index_path, snp_index_name, out, threads, snp_tolerance, keep_temp, mismatches, map_round)

	# Coverage and plots
	cov_table, filtered_list = getCoverage.getCoverage(coverage_bed, coverageData, out, max_multi, min_cov)
	getCoverage.plotCoverage(out, mito_trnas)

	# if remap and snp_tolerance are enabled, skip further analyses, find new mods, and redo alignment and coverage
	if remap and snp_tolerance:
		new_mods = mmQuant.generateModsTable(coverageData, out, threads, cov_table, mismatch_dict, filtered_list, cca, remap, misinc_thresh, mod_lists)
		tRNAtools.newModsParser(out, name, new_mods, mod_lists, tRNA_dict)
		tRNAtools.generateSNPIndex(name, out, snp_tolerance)
		map_round = 2
		bams_list, coverageData = tRNAmap.mainAlign(sample_data, name, genome_index_path, genome_index_name, \
			snp_index_path, snp_index_name, out, threads, snp_tolerance, keep_temp, mismatches, map_round)
		cov_table, filtered_list = getCoverage.getCoverage(coverage_bed, coverageData, out, max_multi, min_cov)
		getCoverage.plotCoverage(out, mito_trnas)
		remap = False

	# featureCounts
	tRNAmap.countReads(bams_list, mode, threads, out)

	# DESeq2
	script_path = os.path.dirname(os.path.realpath(__file__))
	sample_data = os.path.abspath(coverageData)

	log.info("\n+----------------------------------------------+\
	\n| Differential expression analysis with DESeq2 |\
	\n+----------------------------------------------+")
	deseq_cmd = "Rscript " + script_path + "/deseq.R " + out + " " + sample_data + " " + control_cond
	subprocess.call(deseq_cmd, shell=True)
	deseq_out = out + "DESeq2"

	log.info("DESeq2 outputs located in: {}".format(deseq_out))

	# Misincorporation analysis
	mmQuant.generateModsTable(coverageData, out, threads, cov_table, mismatch_dict, filtered_list, cca, remap, misinc_thresh, mod_lists)
	# Output modification context file for plotting
	cons_mod_pos, mod_sites = ssAlign.modContext(out)
	# plot mods and stops
	log.info("Plotting modification and RT stop data...")
	modplot_cmd = "Rscript " + script_path + "/modPlot.R " + out + " " + str(cons_mod_pos) + " " + str(mod_sites)
	subprocess.call(modplot_cmd, shell=True)
	# CCA analysis (see mmQuant.generateModsTable and mmQuant.countMods_mp for initial counting of CCA vs CC ends)
	if cca:
		CCAanalysis.plotDinuc(out)

	# tidy files
	tRNAtools.tidyFiles(out, cca)


if __name__ == '__main__':
	
	################### 
	# Parse arguments #
	################### 
	
	parser = argparse.ArgumentParser(description = 'Custom high-throughput tRNA sequencing alignment and quantification pipeline\
		based on modifications and misincorporation.', add_help = True, usage = "%(prog)s [options] sample_data")

	inputs = parser.add_argument_group("Input files")
	inputs.add_argument('-t', '--trnas', metavar='genomic tRNAs', required = True, dest = 'trnas', help = \
		'Genomic tRNA fasta file, e.g. from gtRNAdb or tRNAscan-SE. Already avalable in data folder for a few model organisms. REQUIRED')
	inputs.add_argument('-o', '--trnaout', metavar = 'tRNA out file', required = True, dest = 'trnaout', help = \
		'tRNA.out file generated by tRNAscan-SE (also may be available on gtRNAdb). Contains information about tRNA features, including introns. REQUIRED')
	inputs.add_argument('-m', '--mito-trnas', metavar = 'mitochondrial tRNAs', required = False, dest = 'mito', \
		help = 'Mitochondrial tRNA fasta file. Should be downloaded from mitotRNAdb for species of interest. Already avaialable in data folder for a few model organisms.')
	
	options = parser.add_argument_group("Program options")
	options.add_argument('--cluster', required = False, dest = 'cluster', action = 'store_true',\
		help = 'Enable usearch sequence clustering of tRNAs by isodecoder - drastically reduces rate of multi-mapping reads.')
	options.add_argument('--cluster_id', metavar = 'clutering id cutoff', dest = 'cluster_id', type = restrictedFloat, nargs = '?', default = 0.95,\
		required = False, help = 'Identity cutoff for usearch clustering between 0 and 1. Default is 0.95.')
	options.add_argument('--threads', metavar = 'thread number', required = False, dest = 'threads', type = int, \
		help = 'Set processor threads to use during read alignment and read counting.')
	options.add_argument('--posttrans_mod_off', required = False, dest = 'posttrans', action = 'store_true', \
		help = "Disable post-transcriptional modification of tRNAs, i.e. addition of 3'-CCA and 5'-G (His) to mature sequences. Disable for certain \
		prokaryotes (e.g. E. coli) where this is genomically encoded. Leave enabled (default) for all eukaryotes.")
	options.add_argument('--control_condition', metavar = 'control condition', required = True, dest = 'control_cond', \
		help = 'Name of control/wild-type condition as per user defined group specified in sample data input. This must exactly match the group name \
		specified in sample data. This is used for differential expression analysis so that results are always in the form mutant/treatment vs WT/control. REQUIRED')
	options.add_argument('--cca_analysis', required = False, dest = 'cca', action = 'store_true',\
		help = "Enable analysis of 3'-CCA ends: Calculates proportions of CC vs CCA ending reads per cluster and performs DESeq2 analysis. \
		Useful for comparing functional to non-funtional mature tRNAs.")

	align = parser.add_argument_group("GSNAP alignment options")
	align.add_argument('--max-mismatches', metavar = 'allowed mismatches', required = False, dest = 'mismatches', type = float, \
		help = 'Maximum mismatches allowed. If specified between 0.0 and 1.0, then trated as a fraction of read length. Otherwise, treated as \
		integer number of mismatches. Default is an automatic ultrafast value calculated by GSNAP; see GSNAP help for more info.')
	align.add_argument('--snp_tolerance', required = False, dest = 'snp_tolerance', action = 'store_true',\
		help = 'Enable GSNAP SNP-tolerant read alignment, where known modifications from Modomics are mapped as SNPs.')


	outputs = parser.add_argument_group("Output options")
	outputs.add_argument('-n', '--name', metavar = 'experiment name', required = True, dest = 'name', help = \
		'Name of experiment. Note, output files and indeces will have this as a prefix. REQUIRED')
	outputs.add_argument('--out_dir', metavar = 'output directory', required = False, dest = 'out', help = \
		'Output directory. Default is current directory. Cannot be an exisiting directory.')
	outputs.add_argument('--keep-temp', required = False, dest='keep_temp', action = 'store_true', help = \
		'Keeps multi-mapping and unmapped bam files from GSNAP alignments. Default is false.')

	featurecounts = parser.add_argument_group("featureCounts options")
	featurecounts.add_argument('--count_mode', metavar = 'featureCounts mode', required = False, dest = 'mode', help = \
		"featureCounts mode to handle reads overlapping more than one feature. Choose from 'none' (multi-overlapping reads are not counted)\
		,'all' (reads are assigned and counted for all overlapping features), or 'fraction' (each overlapping feature receives a fractional count\
		of 1/y, where y is the number of features overlapping with the read). Default is 'none'",\
	 	choices = ['none','all','fraction'])

	bedtools = parser.add_argument_group("Bedtools coverage options")
	bedtools.add_argument('--min_cov', metavar = 'Minimum coverage per cluster', required = False, dest = 'min_cov', type = int, \
		help = "Minimum coverage per cluster to include this cluster in coverage plots, modification analysis, and 3'-CCA analysis. Clusters with \
		 less than this will be filtered out of these analyses. Note that all clusters are included for differential expression analysis with DESeq2.")
	bedtools.add_argument('--max_multi', metavar = 'Bedtools coverage multhreading', required = False, dest = 'max_multi', type = int, \
		help = 'Maximum number of bam files to run bedtools coverage on simultaneously. Increasing this number reduces processing time\
		by increasing number of files processed simultaneously. However, depending on the size of the bam files to process and\
		available memory, too many files processed at once can cause termination of mim-tRNAseq due to insufficient memory. If\
		mim-tRNAseq fails during coverage calculation, lower this number. Increase at your own discretion. Default is 3.')

	remapping = parser.add_argument_group("Analysis of unannotated modifications and realignment")
	remapping.add_argument('--remap', required = False, dest = 'remap', action = 'store_true',\
		help = 'Enable detection of unannotated (potential) modifications from misincorporation data. These are defined as having a total misincorporation rate\
		higher than the threshold set with --misinc_thresh. These modifications are then appended to already known ones, and read alignment is reperformed.\
		Very useful for poorly annotated species in Modomics. Due to realignment and misincorporation parsing, enabling this option slows the analysis down considerably.')
	remapping.add_argument('--misinc_thresh', metavar = 'threshold for unannotated mods', dest = 'misinc_thresh', type = restrictedFloat, nargs = '?', default = 0.1,\
		required = False, help = 'Threshold of total misincorporation rate at a position in a cluster used to call unannotated modifications. Value between 0 and 1, default is 0.1  (10%% misincorporation).')

	parser.add_argument('sample_data', help = 'Sample data sheet in text format, tab-separated. Column 1: full path to fastq (or fastq.gz). Column 2: condition/group.')

	parser.set_defaults(threads=1, out="./", mode = 'none', max_multi = 3, min_cov = 0, mito = '')


	#############################
	# Print help or run mim-seq #
	#############################

	if len(sys.argv[1:]) == 0:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print("     Modification-induced misincorporation sequencing of tRNAs\n")
		parser.print_help()
		parser.exit()
	if len(sys.argv) <= 1:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print("     Modification-induced misincorporation sequencing of tRNAs\n")
		parser.print_usage()
		sys.exit(1)
	else:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print("     Modification-induced misincorporation sequencing of tRNAs\n")
		args = parser.parse_args()
		# Check that control_cond exists in sample data
		conditions = list()
		with open(args.sample_data, "r") as sampleData:
			for line in sampleData:
				line = line.strip()
				if not line.startswith("#"):
					conditions.append(line.split("\t")[1])
		if args.control_cond not in conditions:
			raise argparse.ArgumentTypeError('{} not a valid condition in {}'.format(args.control_cond, args.sample_data))
		else:
			mimseq(args.trnas, args.trnaout, args.name, args.out, args.cluster, args.cluster_id, \
				args.posttrans, args.control_cond, args.threads, args.max_multi, args.snp_tolerance, \
				args.keep_temp, args.mode, args.cca, args.min_cov, args.mismatches, args.remap, args.misinc_thresh, args.mito, args.sample_data)
