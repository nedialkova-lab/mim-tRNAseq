#!/usr/bin/env python3

######################################################
# Wrapper functions for read aligmnent and placement #
######################################################

import subprocess, os, sys, re, logging, copy
from pybedtools import BedTool
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
	snp_index_name, out_dir, threads, snp_tolerance, keep_temp, mismatches, map_round):

	if map_round == 2:
		log.info("\n+------------------+ \
			\n| Realigning reads |\
			\n+------------------+")

	elif map_round == 1:
		log.info("\n+-----------+ \
			\n| Alignment |\
			\n+-----------+")

	# Read sampleData 
	sampleDict = defaultdict()
	unique_bam_list = list()
	alignstats_total = defaultdict(list)
	coverageData = open(out_dir + sampleData.split("/")[-1].split(".")[0] + "_cov." + sampleData.split(".")[-1], "w")
	with open(sampleData, "r") as sampleData:
		for line in sampleData:
			line = line.strip()
			if not line.startswith("#"):
				log.info("**** {} ****".format(line.split("\t")[0].split("/")[-1]))
				fq = line.split("\t")[0]
				name = line.split("\t")[0].split("/")[-1].split(".")[0]
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

def remap(fq, genome_index_path, genome_index_name, snp_index_path, \
	snp_index_name, threads, out_dir, snp_tolerance, keep_temp, mismatches):
# remapping of multi and unmapped reads from round 1

	# generate fastq from unmaped reads
	nomap = out_dir + fq.split("/")[-1].split(".")[0] + ".nomapping.bam"
	nomap_bed = BedTool(nomap)
	nomap_fastq = "".join(nomap.split(".bam")[:-1]) + ".fastq"
	nomap_bed.bam_to_fastq(fq = nomap_fastq)

	# generate fastq from multi-mapping reads
	# first generate sam of primary alignments from multi-mappers
	multi = out_dir + fq.split("/")[-1].split(".")[0] + ".unpaired_mult.bam"
	cmd = "samtools view -@ " + str(threads) + " -F 0x904 -o " + out_dir + "temp_multi.bam " + multi
	subprocess.call(cmd, shell = True)
	multi_bed = BedTool(out_dir + "temp_multi.bam")
	multi_fastq = "".join(multi.split(".bam")[:-1]) + ".fastq"
	multi_bed.bam_to_fastq(fq = multi_fastq)
	os.remove(out_dir + "temp_multi.bam")

	if not keep_temp:
		os.remove(nomap)
		os.remove(multi)

	if not mismatches == None:
		mismatch_string = "--max-mismatches " + str(mismatches) + " "
	elif mismatches == None:
		mismatch_string = "" 

	output_prefix = fq.split("/")[-1].split(".")[0] + "_remap"
	if snp_tolerance:
 		map_cmd = "gsnap -D " + genome_index_path + " -d " + genome_index_name + " -V " + snp_index_path + " -v " \
 		+ snp_index_name + " -t " + str(threads) + " --split-output " + out_dir + output_prefix + \
 		" --format=sam --genome-unk-mismatch=0 --md-lowercase-snp  --ignore-trim-in-filtering 1 --force-single-end " + mismatch_string +\
 		nomap_fastq + " " + multi_fastq + " &>> " + out_dir + "remap_align.log"
	else:
 		map_cmd = "gsnap -D " + genome_index_path + " -d " + genome_index_name + " -t " + str(threads) + \
 		" --split-output " + out_dir + output_prefix + " --format=sam --genome-unk-mismatch=0 --md-lowercase-snp --ignore-trim-in-filtering 1 --force-single-end " + mismatch_string + " " + \
 		nomap_fastq + " " + multi_fastq + " &>> " + out_dir + "remap_align.log"

	log.info("Realigning multi-mapped and unmapped reads to {} with updated SNP index...".format(genome_index_name))
	subprocess.call(map_cmd, shell = True)

	os.remove(multi_fastq)
	os.remove(nomap_fastq)
	
	# remove transloc sam output if no reads present (often the case)
	cmd = "samtools view -c " + out_dir + output_prefix + ".unpaired_transloc"
	readcount = int(subprocess.check_output(cmd, shell = True))
	if readcount == 0:
		os.remove(out_dir + output_prefix + ".unpaired_transloc")

	# write mapping stats and compress to bam - remove mutlimapping and unmapped unless keep_temp = True
	log.info("Compressing SAM files, sorting, and computing mapping stats...")
	with open(out_dir + "mapping_stats.txt","a") as stats_out:
		align_pathlist = Path(out_dir).glob(output_prefix + "*")
		for file in align_pathlist:
			if re.search("mult",file.name) and not re.search("bam", file.name):
				cmd = "samtools view -@ " + str(threads) + " -F 0x904 -c " + out_dir + file.name
				multi_count = int(subprocess.check_output(cmd, shell = True))
				if keep_temp:
					cmd = "samtools view -@ " + str(threads) + " -bh -o " + out_dir + file.name + ".bam " + out_dir + file.name
					subprocess.call(cmd, shell = True)
					os.remove(out_dir + file.name)
				else:
					os.remove(out_dir + file.name)
			elif re.search("uniq",file.name) and not re.search("bam", file.name):
				cmd = "samtools view -@ " + str(threads) + " -bh " + out_dir + file.name + " | samtools sort -@ " + str(threads) + " -o " + out_dir + file.name + ".bam" + " - "
				subprocess.call(cmd, shell = True)
				os.remove(out_dir + file.name)
				# merge 1st run bam and remapped bam
				merged_bam = out_dir + fq.split("/")[-1].split(".")[0] + ".unpaired_uniq_remapMerge.bam"
				cmd = "samtools merge " + merged_bam + " " + out_dir + file.name + ".bam " + out_dir + fq.split("/")[-1].split(".")[0] + ".unpaired_uniq.bam"
				subprocess.call(cmd, shell = True) 
				#os.remove(out_dir + file.name + ".bam")  
				cmd = "samtools view -@ " + str(threads) + " -c " + merged_bam
				unique_count = int(subprocess.check_output(cmd, shell = True))
				unique_bam = merged_bam

			elif re.search("nomapping",file.name) and not re.search("bam", file.name):
				cmd = "samtools view -@ " + str(threads) + " -c " + out_dir + file.name
				unmapped_count = int(subprocess.check_output(cmd, shell = True))
				if keep_temp:
					cmd = "samtools view -@ " + str(threads) + " -bh -o " + out_dir + file.name + ".bam " + out_dir + file.name
					subprocess.call(cmd, shell = True)
					os.remove(out_dir + file.name)
				else:
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
		snp_index_name, "-t", str(threads), "--split-output", out_dir + output_prefix, \
		"--format", "sam", "--genome-unk-mismatch", "0", "--md-lowercase-snp", "--ignore-trim-in-filtering", "1", fq]
		map_cmd = list(filter(None, map_cmd))
		map_cmd[-1:-1] = mismatch_list
	else:
		map_cmd = ["gsnap", zipped, "-D", genome_index_path, "-d", genome_index_name, "-t", str(threads),\
		"--split-output", out_dir + output_prefix, "--format", "sam", "--genome-unk-mismatch", "0", "--md-lowercase-snp", "--ignore-trim-in-filtering", "1", fq]
		map_cmd = list(filter(None, map_cmd))
		map_cmd[-1:-1] = mismatch_list

	log.info("Aligning reads to {}...".format(genome_index_name))
	subprocess.check_call(map_cmd, stderr = open(out_dir + "align.log", "a"))
	
	# remove transloc sam output if no reads present (often the case)
	cmd = ["samtools", "view", "-c", out_dir + output_prefix + ".unpaired_transloc"]
	readcount = int(subprocess.check_output(cmd))
	if readcount == 0:
		os.remove(out_dir + output_prefix + ".unpaired_transloc")

	# write mapping stats and compress to bam - remove mutlimapping and unmapped unless keep_temp = True
	log.info("Compressing SAM files, sorting, and computing mapping stats...")
	with open(out_dir + "mapping_stats.txt","a") as stats_out:
		align_pathlist = Path(out_dir).glob(output_prefix + "*")
		for file in align_pathlist:
			if re.search("mult",file.name) and not re.search("bam", file.name):
				cmd = ["samtools", "view" ,"-@", str(threads), "-F", "0x904", "-c", out_dir + file.name]
				multi_count = int(subprocess.check_output(cmd))
				if keep_temp or remap:
					cmd = ["samtools", "view", "-@", str(threads), "-bh", "-o", out_dir + file.name + ".bam", out_dir + file.name]
					subprocess.check_call(cmd)
					os.remove(out_dir + file.name)
				elif not keep_temp or not remap:
					os.remove(out_dir + file.name)
			elif re.search("uniq",file.name) and not re.search("bam", file.name):
				cmd = ["samtools", "view", "-@", str(threads), "-c", out_dir + file.name]
				unique_count = int(subprocess.check_output(cmd))
				unique_bam = out_dir + file.name + ".bam"
				ps = subprocess.Popen(["samtools", "view", "-@", str(threads), "-bh", out_dir + file.name], stdout = subprocess.PIPE)
				cmd = ["samtools", "sort", "-@", str(threads), "-o", unique_bam]
				subprocess.check_call(cmd, stdin = ps.stdout)
				os.remove(out_dir + file.name)
			elif re.search("nomapping",file.name) and not re.search("bam", file.name):
				cmd = ["samtools", "view", "-@", str(threads), "-c", out_dir + file.name]
				unmapped_count = int(subprocess.check_output(cmd))
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

if __name__ == '__main__':
	mainAlign(sys.argv[1:])

