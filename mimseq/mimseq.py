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
# github: https://github.com/nedialkova-lab/mim-tRNAseq

from __future__ import absolute_import
from . import version
from .tRNAtools import modsToSNPIndex, generateGSNAPIndices, newModsParser, tidyFiles
from .tRNAmap import mainAlign
from .getCoverage import getCoverage, plotCoverage
from .mmQuant import generateModsTable, plotCCA
from .ssAlign import structureParser, modContext 
from .splitClusters import splitIsodecoder, unsplitClusters, getIsodecoderSizes, writeIsodecoderTranscripts
import sys, os, subprocess, logging, datetime, copy
import argparse
from pyfiglet import figlet_format
from collections import defaultdict

log = logging.getLogger(__name__)

def restrictedFloat(x):
## Method for restricting cluster_id and cov_diff argument to float between 0 and 1
	try:
		x = float(x)
		if x < 0.0 or x > 1.0:
			raise argparse.ArgumentTypeError('{} not in range 0.0 - 1.0'.format(x))
		return x
	except ValueError:
		raise argparse.ArgumentTypeError('{} not a real number'.format(x))

def restrictedFloat2(x):
## Method for restricting min-cov argument to float between 0 and 1, or int greater than 1

	try:
		x = float(x)
		if x < 0.0:
			raise argparse.ArgumentTypeError('{} not greater than 0'.format(x))
		return x
	except ValueError:
		raise argparse.ArgumentTypeError('{} not a real number'.format(x))

def mimseq(trnas, trnaout, name, species, out, cluster, cluster_id, cov_diff, posttrans, control_cond, threads, max_multi, snp_tolerance, \
	keep_temp, cca, double_cca, min_cov, mismatches, remap, remap_mismatches, misinc_thresh, mito_trnas, pretrnas, local_mod, sample_data):
	
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
			logging.StreamHandler()
		])
	log.info("mim-tRNAseq v{} run with command:".format(version.__version__))
	log.info(" ".join(sys.argv))

	########
	# main #
	########

	map_round = 1 #first round of mapping

	# Parse tRNA and modifications, generate SNP index
	modifications = os.path.dirname(os.path.realpath(__file__))
	modifications += "/modifications"
	coverage_bed, snp_tolerance, mismatch_dict, insert_dict, del_dict, mod_lists, Inosine_lists, Inosine_clusters, tRNA_dict, cluster_dict, cluster_perPos_mismatchMembers \
	= modsToSNPIndex(trnas, trnaout, mito_trnas, modifications, name, out, double_cca, snp_tolerance, cluster, cluster_id, posttrans, pretrnas, local_mod)
	structureParser()
	# Generate GSNAP indices
	genome_index_path, genome_index_name, snp_index_path, snp_index_name = generateGSNAPIndices(species, name, out, map_round, snp_tolerance, cluster)

	# Align
	bams_list, coverageData = mainAlign(sample_data, name, genome_index_path, genome_index_name, \
		snp_index_path, snp_index_name, out, threads, snp_tolerance, keep_temp, mismatches, map_round)

	# define unique mismatches/insertions to assign reads to unique tRNA sequences
	unique_isodecoderMMs = defaultdict(dict)
	splitBool = list()
	newSplitBool = list()
	if cluster and cluster_id != 1:
		cluster_dict2 = copy.deepcopy(cluster_dict) # copy so splitReadsIsodecoder does not edit main cluster_dict
		unique_isodecoderMMs, splitBool, isodecoder_sizes = splitIsodecoder(cluster_perPos_mismatchMembers, insert_dict, del_dict, tRNA_dict, cluster_dict2, out, name)
		unsplit = unsplitClusters(coverageData, coverage_bed, unique_isodecoderMMs, threads, cov_diff)
		newSplitBool = list(set(splitBool).union(unsplit))
	elif cluster and cluster_id == 1:
		isodecoder_sizes = {iso:len(members) for iso, members in cluster_dict.items()}
		writeIsodecoderTranscripts(out, name, cluster_dict, tRNA_dict)
	elif not cluster:
		isodecoder_sizes = getIsodecoderSizes(out, name, tRNA_dict)

	# if remap and snp_tolerance are enabled, skip further analyses, find new mods, and redo alignment and coverage
	if remap and (snp_tolerance or not mismatches == 0.0):
		new_mods, new_Inosines, filtered_cov, filter_warning = generateModsTable(coverageData, out, name, threads, min_cov, mismatch_dict, insert_dict, del_dict, cluster_dict, cca, remap, misinc_thresh, mod_lists, Inosine_lists, tRNA_dict, Inosine_clusters, unique_isodecoderMMs, newSplitBool, isodecoder_sizes, cluster)
		Inosine_clusters, snp_tolerance, newtRNA_dict, new_mod_lists = newModsParser(out, name, new_mods, new_Inosines, mod_lists, Inosine_lists, tRNA_dict, cluster, remap, snp_tolerance)
		map_round = 2
		genome_index_path, genome_index_name, snp_index_path, snp_index_name = generateGSNAPIndices(species, name, out, map_round, snp_tolerance, cluster)
		bams_list, coverageData = mainAlign(sample_data, name, genome_index_path, genome_index_name, \
			snp_index_path, snp_index_name, out, threads, snp_tolerance, keep_temp, remap_mismatches, map_round)
		remap = False
	#else:
	#	log.info("\n*** New modifications not discovered as remap is not enabled ***\n")

	# redo checks for unsplit isodecoders based on coverage
	if map_round == 2 and cluster_id != 1:
		unsplit = unsplitClusters(coverageData, coverage_bed, unique_isodecoderMMs, threads, cov_diff)
		newSplitBool = list(set(splitBool).union(unsplit))

	# Misincorporation analysis
	filter_warning = False
	filtered_cov = list()
	if snp_tolerance or not mismatches == 0.0:
		if 'newtRNA_dict' in locals():
			new_mods, new_Inosines, filtered_cov, filter_warning = generateModsTable(coverageData, out, name, threads, min_cov, mismatch_dict, insert_dict, del_dict, cluster_dict, cca, remap, misinc_thresh, new_mod_lists, Inosine_lists, newtRNA_dict, Inosine_clusters, unique_isodecoderMMs, newSplitBool, isodecoder_sizes, cluster)
		else:
			new_mods, new_Inosines, filtered_cov, filter_warning = generateModsTable(coverageData, out, name, threads, min_cov, mismatch_dict, insert_dict, del_dict, cluster_dict, cca, remap, misinc_thresh, mod_lists, Inosine_lists, tRNA_dict, Inosine_clusters, unique_isodecoderMMs, newSplitBool, isodecoder_sizes, cluster)

	else:
		log.info("*** Misincorporation analysis not possible; either --snp-tolerance must be enabled, or --max-mismatches must not be 0! ***\n")

	# Output modification context file for plotting
	mod_sites, cons_pos_list, cons_pos_dict = modContext(out)

	script_path = os.path.dirname(os.path.realpath(__file__))
	
	if snp_tolerance or not mismatches == 0.0:
					# plot mods and stops, catch exception with command call and print log error if many clusters are filtered (known to cause issues with R code handling mods table)
		log.info("Plotting modification and RT stop data...")
		try:
			modplot_cmd = ["Rscript", script_path + "/modPlot.R", out, str(mod_sites), str(cons_pos_list), str(misinc_thresh), str(mito_trnas), control_cond]
			process = subprocess.Popen(modplot_cmd, stdout = subprocess.PIPE)
			while True:
				line = process.stdout.readline()
				if not line:
					break
				line = line.decode("utf-8")
				log.info(line.rstrip())
			exitcode = process.wait()
		except subprocess.CalledProcessError:
			if filter_warning:
				log.error("Error plotting modifications. Potentially caused by any clusters filtered by --min-cov: lower --min-cov or assess data quality and sequencing depth!")
				raise
		# CCA analysis (see mmQuant.generateModsTable and mmQuant.countMods_mp for initial counting of CCA vs CC ends)
		if cca:
			plotCCA(out, double_cca)

	# Coverage and plots
	sorted_aa = getCoverage(coverageData, out, control_cond, filtered_cov)
	plotCoverage(out, mito_trnas, sorted_aa)

	# DESeq2
	sample_data = os.path.abspath(coverageData)

	log.info("\n+----------------------------------------------+\
	\n| Differential expression analysis with DESeq2 |\
	\n+----------------------------------------------+")

	deseq_cmd = ["Rscript", script_path + "/deseq.R", out, sample_data, control_cond, str(cluster_id)]
	#subprocess.check_call(deseq_cmd)
	process = subprocess.Popen(deseq_cmd, stdout = subprocess.PIPE)
	while True:
		line = process.stdout.readline()
		if not line:
			break
		line = line.decode("utf-8")
		log.info(line.rstrip())
	exitcode = process.wait()
	deseq_out = out + "DESeq2"
	log.info("DESeq2 outputs located in: {}".format(deseq_out))

	# tidy files
	tidyFiles(out, cca)

def main():

	################### 
	# Parse arguments #
	################### 
	
	parser = argparse.ArgumentParser(description = 'Custom high-throughput tRNA sequencing alignment and quantification pipeline\
		based on modification induced misincorporation cDNA synthesis.', add_help = True, usage = "%(prog)s [options] sample data")

	inputs = parser.add_argument_group("Input files")
	inputs.add_argument('-s','--species', metavar='species', required = not ('-t' in sys.argv), dest = 'species', help = \
		'Species being analyzed for which to load pre-packaged data files (prioritized over -t, -o and -m). Options are: Hsap, Hsap38, Mmus, Scer, Spom, Dmel, Drer, Ecol', \
		choices = ['Hsap','Hsap38','Mmus','Scer','Spom','Dmel', 'Drer', 'Ecol'])
	inputs.add_argument('-t', '--trnas', metavar='genomic tRNAs', required = False, dest = 'trnas', help = \
		'Genomic tRNA fasta file, e.g. from gtRNAdb or tRNAscan-SE. Already avalable in data folder for a few model organisms.')
	inputs.add_argument('-o', '--trnaout', metavar = 'tRNA out file', required = (not '--species' or '-s' in sys.argv) or ('-t' in sys.argv), 
		dest = 'trnaout', help = 'tRNA.out file generated by tRNAscan-SE (also may be available on gtRNAdb). Contains information about tRNA features, including introns.')
	inputs.add_argument('-m', '--mito-trnas', metavar = 'mitochondrial tRNAs', required = False, dest = 'mito', \
		help = 'Mitochondrial tRNA fasta file. Should be downloaded from mitotRNAdb for species of interest. Already available in data folder for a few model organisms.')
	
	options = parser.add_argument_group("Program options")
	options.add_argument('--pretRNAs', required = False, dest = 'pretrnas', action = 'store_true',\
		help = "Input reference sequences are pretRNAs. Enabling this option will disable the removal of intron sequences and addition of 3'-CCA to generate \
		mature tRNA sequences. Useful for mapping and discovering pretRNA sequence reads.")
	options.add_argument('--no-cluster', required = False, dest = 'cluster', action = 'store_false',\
		help = 'Disable usearch sequence clustering of tRNAs by isodecoder which drastically reduces the rate of multi-mapping reads. Default is enabled.')
	options.add_argument('--cluster-id', metavar = 'clustering identity threshold', dest = 'cluster_id', type = restrictedFloat, nargs = '?', default = 0.97,\
		required = False, help = 'Identity cutoff for usearch clustering between 0 and 1. Default is 0.97.')
	options.add_argument('--deconv-cov-ratio', metavar='deconvolution coverage threshold', dest='cov_diff', type = restrictedFloat, nargs = '?', default=0.5,\
		required=False, help="Threshold for ratio between coverage at 3' end and mismatch used for deconvolution. Coverage reductions greater than the threshold will result in non-deconvoluted sequences. \
			Default is 0.5 (i.e. less than 50%% reduction required for deconvolution).")
	options.add_argument('--threads', metavar = 'thread number', required = False, dest = 'threads', type = int, \
		help = 'Set processor threads to use during read alignment and read counting.')
	options.add_argument('--posttrans-mod-off', required = False, dest = 'posttrans', action = 'store_true', \
		help = "Disable post-transcriptional modification of tRNAs, i.e. addition of 3'-CCA and 5'-G (His) to mature sequences. Disable for certain \
		prokaryotes (e.g. E. coli) where this is genomically encoded. Leave enabled (default) for all eukaryotes.")
	options.add_argument('--control-condition', metavar = 'control condition', required = True, dest = 'control_cond', \
		help = 'Name of control/wild-type condition as per user defined group specified in sample data input. This must exactly match the group name \
		specified in sample data. This is used for differential expression analysis so that results are always in the form mutant/treatment vs WT/control. REQUIRED')
	options.add_argument('--no-cca-analysis', required = False, dest = 'cca', action = 'store_false',\
		help = "Disable analysis of 3'-CCA ends. When enabled, this calculates proportions of CC vs CCA ending reads per cluster and performs DESeq2 analysis. \
		Useful for comparing functional to non-functional mature tRNAs. Default is enabled.")
	options.add_argument('--double-cca', required = False, dest = 'double_cca', action = "store_true",\
		help = "Enable analysis of 3'-CCACCA tagging for tRNA degradation pathway. Note that this will alter the output of the CCA analysis pipeline.")
	options.add_argument('--local-modomics', required=False, dest = 'local_mod', action='store_true',\
		help = "Disable retrieval of Modomics data from online. Instead use older locally stored data. Warning - this leads\
			to usage of older Modomics data!")

	align = parser.add_argument_group("GSNAP alignment options")
	align.add_argument('--max-mismatches', metavar = 'allowed mismatches', required = False, dest = 'mismatches', type = float, \
		help = 'Maximum mismatches allowed. If specified between 0.0 and 1.0, then treated as a fraction of read length. Otherwise, treated as \
		integer number of mismatches. Default is an automatic ultrafast value calculated by GSNAP; see GSNAP help for more info.')
	align.add_argument('--remap-mismatches', metavar = 'allowed mismatches for remap', required = False, dest = 'remap_mismatches', type = float,\
		help = 'Maximum number of mismatches allowed during remapping of all reads. Treated similarly to --max-mismatches. This is important to control misalignment of reads to similar clusters/tRNAs \
		Note that the SNP index will be updated with new SNPs from the first round of alignment and so this should be relatively small to prohibit misalignment.')
	align.add_argument('--no-snp-tolerance', required = False, dest = 'snp_tolerance', action = 'store_false',\
		help = 'Disable GSNAP SNP-tolerant read alignment, where known modifications from Modomics are mapped as SNPs. Default is enabled.')


	outputs = parser.add_argument_group("Output options")
	outputs.add_argument('-n', '--name', metavar = 'experiment name', required = True, dest = 'name', help = \
		'Name of experiment. Note, output files and indices will have this as a prefix. REQUIRED')
	outputs.add_argument('--out-dir', metavar = 'output directory', required = False, dest = 'out', help = \
		'Output directory. Default is current directory. Cannot be an existing directory.')
	outputs.add_argument('--keep-temp', required = False, dest='keep_temp', action = 'store_true', help = \
		'Keeps multi-mapping and unmapped bam files from GSNAP alignments. Default is false.')

	bedtools = parser.add_argument_group("Bedtools coverage options")
	bedtools.add_argument('--min-cov', metavar = 'Minimum coverage per cluster', required = False, dest = 'min_cov', type = restrictedFloat2, default=0.0005, \
		help = "Minimum coverage per cluster required to include this cluster in coverage plots, modification analysis, and 3'-CCA analysis. \
		Can be a fraction of total mapped reads between 0 and 1, or an integer of absolute coverage. Any cluster not meeting the threshold in 1 or more sample will be excluded. \
		Note that all clusters are included for differential expression analysis with DESeq2. Default = 0.0005 (0.05% mapped reads).")
	bedtools.add_argument('--max-multi', metavar = 'Bedtools coverage multithreading', required = False, dest = 'max_multi', type = int, \
		help = 'Maximum number of bam files to run bedtools coverage on simultaneously. Increasing this number reduces processing time\
		by increasing number of files processed simultaneously. However, depending on the size of the bam files to process and\
		available memory, too many files processed at once can cause termination of mim-tRNAseq due to insufficient memory. If\
		mim-tRNAseq fails during coverage calculation, lower this number. Increase at your own discretion. Default is 3.')

	remapping = parser.add_argument_group("Analysis of unannotated modifications and realignment")
	remapping.add_argument('--remap', required = False, dest = 'remap', action = 'store_true',\
		help = 'Enable detection of unannotated (potential) modifications from misincorporation data. These are defined as having a total misincorporation rate\
		higher than the threshold set with --misinc-thresh. These modifications are then appended to already known ones, and read alignment is reperformed.\
		Very useful for poorly annotated species in Modomics. Due to realignment and misincorporation parsing, enabling this option slows the analysis down considerably.')
	remapping.add_argument('--misinc-thresh', metavar = 'threshold for unannotated mods', dest = 'misinc_thresh', type = restrictedFloat, nargs = '?', default = 0.1,\
		required = False, help = 'Threshold of total misincorporation rate at a position in a cluster used to call unannotated modifications. Value between 0 and 1, default is 0.1  (10%% misincorporation).')

	parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version.__version__), help = 'Show version number and exit')
	parser.add_argument('sampledata', help = 'Sample data sheet in text format, tab-separated. Column 1: full path to fastq (or fastq.gz). Column 2: condition/group.')
	
	parser.set_defaults(threads=1, out="./", max_multi = 3, min_cov = 0, mito = '', cov_diff = 0.5)

	#########################################
	# Print help, check args or run mim-seq #
	#########################################

	if len(sys.argv[1:]) == 0:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print(" Modification-induced misincorporation analysis of tRNA sequencing data\n")
		parser.print_help()
		parser.exit()
	if len(sys.argv) <= 1:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print(" Modification-induced misincorporation analysis of tRNA sequencing data\n")
		parser.print_usage()
		sys.exit(1)
	else:
		print(figlet_format('mim-tRNAseq', font='standard'))
		print(" Modification-induced misincorporation analysis of tRNA sequencing data\n")
		args = parser.parse_args()
		if args.pretrnas:
			if args.cca:
				log.warning("Disabling CCA analysis in pre-tRNA mode...")
				args.cca = False
			if args.cluster:
				log.warning("Disabling tRNA clustering in pre-tRNA mode...")
				args.cluster = False
		# Check that control_cond exists in sample data
		conditions = list()
		with open(args.sampledata, "r") as sampleData:
			for line in sampleData:
				line = line.strip()
				if not line.startswith("#"):
					conditions.append(line.split("\t")[1])
		if args.control_cond not in conditions:
			raise argparse.ArgumentTypeError('{} not a valid condition in {}'.format(args.control_cond, args.sampledata))
		if not args.species and not (args.trnas or args.trnaout):
			parser.error('Must specify valid --species argument or supply -t (tRNA sequences) and -o (tRNAscan out file)!')						
		else:
			if args.species:
				if args.species == 'Hsap':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/hg19-eColitK/hg19_eColitK.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/hg19-eColitK/hg19_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/hg19-eColitK/hg19-mitotRNAs.fa"
				if args.species == 'Hsap38':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/hg38-eColitK/hg38-tRNAs-filtered.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/hg38-eColitK/hg38-tRNAs-detailed.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/hg38-eColitK/hg38-mitotRNAs.fa"
				if args.species == 'Scer':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/sacCer3-eColitK/sacCer3_eschColitK.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/sacCer3-eColitK/sacCer3_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/sacCer3-eColitK/sacCer3-mitotRNAs.fa"
				if args.species == 'Mmus':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/mm10-eColitK/mm10_eColitK-tRNAs.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/mm10-eColitK/mm10_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/mm10-eColitK/mm10-mitotRNAs.fa"
				if args.species == 'Spom':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/schiPomb-eColitK/schiPomb_972H-tRNAs.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/schiPomb-eColitK/schiPomb_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/schiPomb-eColitK/schiPomb-mitotRNAs.fa"
				if args.species == 'Dmel':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/dm6-eColitK/dm6_eColitK-tRNAs.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/dm6-eColitK/dm6_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/dm6-eColitK/dm6-mitotRNAs.fa"
				if args.species == 'Drer':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/danRer11-eColitK/danRer11_eColitK_filtered.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/danRer11-eColitK/danRer11_eschColi-tRNAs.out"
					args.mito = os.path.dirname(os.path.realpath(__file__)) + "/data/danRer11-eColitK/danRer11-mitotRNAs.fa"
				if args.species == 'Ecol':
					args.trnas = os.path.dirname(os.path.realpath(__file__)) + "/data/eschColi-K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.fa"
					args.trnaout = os.path.dirname(os.path.realpath(__file__)) + "/data/eschColi-K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.out"
					args.mito = ''
			else:
				args.species = args.trnas.split("/")[-1].split(".")[0]
			mimseq(args.trnas, args.trnaout, args.name, args.species, args.out, args.cluster, args.cluster_id, args.cov_diff, \
				args.posttrans, args.control_cond, args.threads, args.max_multi, args.snp_tolerance, \
				args.keep_temp, args.cca, args.double_cca, args.min_cov, args.mismatches, args.remap, args.remap_mismatches, \
				args.misinc_thresh, args.mito, args.pretrnas, args.local_mod, args.sampledata)

if __name__ == '__main__':
	main()
	
