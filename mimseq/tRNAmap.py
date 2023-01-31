#!/usr/bin/env python3

######################################################
# Wrapper functions for read aligmnent and placement #
######################################################

import subprocess, os, re, logging, copy
from pybedtools import BedTool
import pysam
from collections import defaultdict
from pathlib import Path
import glob
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import seaborn as sns

log = logging.getLogger(__name__)

def mainAlign(sampleData, experiment_name, genome_index_path, genome_index_name, snp_index_path, \
	snp_index_name, out_dir, threads, snp_tolerance, keep_temp, mismatches, map_round, remap):

	if map_round == 2:
		log.info("\n+------------------+ \
			\n| Realigning reads |\
			\n+------------------+")

	elif map_round == 1:
		log.info("\n+-----------+ \
			\n| Alignment |\
			\n+-----------+")

	# Read sampleData 
	unique_bam_list = list()
	alignstats_total = defaultdict(list)
	coverageData = open(out_dir + sampleData.split("/")[-1].split(".")[0] + "_cov." + sampleData.split(".")[-1], "w")
	with open(sampleData, "r") as sampleData:
		for line in sampleData:
			line = line.strip()
			if not line.startswith("#"):
				log.info("**** {} ****".format(line.split("\t")[0].split("/")[-1]))
				fq = line.split("\t")[0]
				group = line.split("\t")[1]

				# align
				if map_round == 1:
					unique_bam, librarySize, alignstats = mapReads(fq, genome_index_path, genome_index_name, snp_index_path, snp_index_name, threads, out_dir, snp_tolerance, keep_temp, mismatches, remap)
				elif map_round == 2:
					with open(out_dir + "mapping_stats.txt","a") as stats_out:
						stats_out.write("** NEW ALIGNMENT **\n\n")
					unique_bam, librarySize, alignstats = mapReads(fq, genome_index_path, genome_index_name, snp_index_path, snp_index_name, threads, out_dir, snp_tolerance, keep_temp, mismatches, remap)
				
				unique_bam_list.append(unique_bam)
				coverageData.write(unique_bam + "\t" + group + "\t" + str(librarySize) + "\n")

				if len(alignstats_total) == 0:
					alignstats_total = copy.deepcopy(alignstats)
				else:
					for key, value in alignstats.items():
						for item in value:
							alignstats_total[key].append(item)

	# alignstats plot and save
	alignstats_df = pd.DataFrame.from_dict(alignstats_total)
	alignstats_df["Proportion"] = alignstats_df.groupby(["Lib"])["Count"].apply(lambda x: x.astype(float)/x.sum())
	g = sns.barplot(data = alignstats_df, x = "Lib", y = "Proportion", hue = "Type", palette = "Set2")
	plt.xticks(rotation = 90)
	plt.tight_layout()
	fig = g.get_figure()
	if map_round ==1:
		alignstats_version = "Primary_"
	elif map_round == 2:
		alignstats_version = "Remap_"
	fig.savefig(out_dir + alignstats_version + "alignstats.pdf")
	plt.close(fig)

	log.info('Alignment statistics and plot saved to {}align/mapping_stats.txt'.format(out_dir))

	coverageData.close()

	return(unique_bam_list, coverageData.name)

def mapReads(fq, genome_index_path, genome_index_name, snp_index_path, snp_index_name, threads, \
	out_dir,snp_tolerance, keep_temp, mismatches, remap):
# map with or without SNP index and report initial map statistics

	# check zip status of input reads for command building
	zipped = ''
	if re.search(".gz",fq):
		zipped = '--gunzip'
		output_prefix = fq.split("/")[-1].split(".fastq.gz")[0]
	else:
		output_prefix = fq.split("/")[-1].split(".fastq")[0]
	
	if not mismatches == None:
		mismatch_list = ["--max-mismatches", str(mismatches)]
	elif mismatches == None:
		mismatch_list = ""

	if snp_tolerance:
		map_cmd = ["gsnap", zipped, "-D", genome_index_path, "-d", genome_index_name, "-V", snp_index_path, "-v", \
			snp_index_name, "-t", str(threads), "--split-output", out_dir + output_prefix, "--format", "sam", \
			"--genome-unk-mismatch", "0", "--md-lowercase-snp", "--ignore-trim-in-filtering", "1", fq]
		map_cmd = list(filter(None, map_cmd))
		map_cmd[-1:-1] = mismatch_list
	else:
		map_cmd = ["gsnap", zipped, "-D", genome_index_path, "-d", genome_index_name, "-t", str(threads), \
			"--split-output", out_dir + output_prefix, "--format", "sam", "--genome-unk-mismatch", "0", \
			"--md-lowercase-snp", "--ignore-trim-in-filtering", "1", fq]
		map_cmd = list(filter(None, map_cmd))
		map_cmd[-1:-1] = mismatch_list

	log.info("Aligning reads to {}...".format(genome_index_name))
	subprocess.check_call(map_cmd, stderr = open(out_dir + "align.log", "a"))
	
	# remove transloc sam output if no reads present (often the case)
	readcount = int(pysam.view("-@",str(threads),"-c", out_dir + output_prefix + ".unpaired_transloc").strip())
	if readcount == 0:
		os.remove(out_dir + output_prefix + ".unpaired_transloc")

	# write mapping stats and compress to bam - remove mutlimapping and unmapped unless keep_temp = True
	log.info("Compressing SAM files, sorting, and computing mapping stats...")
	with open(out_dir + "mapping_stats.txt","a") as stats_out:
		align_pathlist = Path(out_dir).glob(output_prefix + "*")
		for file in align_pathlist:
			if re.search("mult",file.name) and not re.search("bam", file.name):
				multi_count = int(pysam.view("-@",str(threads),"-F", "0x904", "-c", out_dir + file.name).strip())
				if keep_temp or remap:
					cmd = ["samtools", "view", "-@", str(threads), "-bh", "-o", out_dir + file.name + ".bam", out_dir + file.name]
					subprocess.check_call(cmd)
					os.remove(out_dir + file.name)
				elif not keep_temp or not remap:
					os.remove(out_dir + file.name)
			elif re.search("uniq",file.name) and not re.search("bam", file.name):
				unique_count = int(pysam.view("-@",str(threads),"-c", out_dir + file.name).strip())
				unique_bam = out_dir + file.name + ".bam"
				ps = subprocess.Popen(["samtools", "view", "-@", str(threads), "-bh", out_dir + file.name], stdout = subprocess.PIPE)
				cmd = ["samtools", "sort", "-@", str(threads), "-o", unique_bam]
				subprocess.check_call(cmd, stdin = ps.stdout)
				index_cmd = ["samtools", "index", "-@", str(threads), unique_bam]
				subprocess.check_call(index_cmd)
				os.remove(out_dir + file.name)
			elif re.search("nomapping",file.name) and not re.search("bam", file.name):
				unmapped_count = int(pysam.view("-@",str(threads),"-c", out_dir + file.name).strip())
				if keep_temp or remap:
					cmd = ["samtools", "view", "-@", str(threads), "-bh", "-o", out_dir + file.name + ".bam", out_dir + file.name]
					subprocess.check_call(cmd)
					os.remove(out_dir + file.name)
				elif not keep_temp or not remap:
					os.remove(out_dir + file.name)

		total_count = unique_count + multi_count + unmapped_count

		stats_out.write("{}\nUniquely mapped reads: {:d} ({:.0%}) \nMulti-mapping reads: {:d} ({:.0%}) \nUnmapped reads: {:d} ({:.0%}) \nTotal: {:d}\n\n"\
			.format(fq.split("/")[-1], unique_count, (unique_count/total_count),multi_count, (multi_count/total_count), unmapped_count, (unmapped_count/total_count), total_count))
	
	alignstats_dict = defaultdict(list)
	type_list = ["Uniquely mapped", "Multi-mapped", "Unmapped"]
	for i, count in enumerate([unique_count, multi_count, unmapped_count]):
		alignstats_dict["Lib"].append(fq.split("/")[-1].split(".")[0])
		alignstats_dict["Type"].append(type_list[i])
		alignstats_dict["Count"].append(count)

	return(unique_bam, unique_count, alignstats_dict)
