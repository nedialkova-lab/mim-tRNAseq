#!/usr/bin/env python3

######################################################
# Wrapper functions for read aligmnent and placement #
######################################################

import subprocess, os, sys, re, logging
from collections import defaultdict
from pathlib import Path
import glob

log = logging.getLogger(__name__)

def mainAlign(sampleData, experiment_name, genome_index_path, genome_index_name, snp_index_path, \
	snp_index_name, out_dir, threads, snp_tolerance, keep_temp):

	log.info("\n+-----------+ \
	\n| Alignment |\
	\n+-----------+")

	# Read sampleData 
	sampleDict = defaultdict()
	unique_bam_list = list()
	coverageData = open(out_dir + sampleData.split(".")[0] + "_cov." + sampleData.split(".")[-1], "w")
	with open(sampleData, "r") as sampleData:
		for line in sampleData:
			line = line.strip()
			if not line.startswith("#"):
				log.info("**** {} ****".format(line.split("\t")[0].split("/")[-1]))
				fq = line.split("\t")[0]
				name = line.split("\t")[0].split("/")[-1].split(".")[0]
				group = line.split("\t")[1]

				# align
				unique_bam, librarySize = mapReads(fq, genome_index_path, genome_index_name, snp_index_path, snp_index_name, threads, out_dir, snp_tolerance, keep_temp)
				unique_bam_list.append(unique_bam)
				coverageData.write(unique_bam + "\t" + group + "\t" + str(librarySize) + "\n")

	coverageData.close()

	return(unique_bam_list, coverageData.name)

def mapReads(fq, genome_index_path, genome_index_name, snp_index_path, snp_index_name, threads, \
	out_dir,snp_tolerance, keep_temp):
# map with or without SNP index and report initial map statistics

	# check zip status of input reads for command building
	zipped = ''
	if re.search(".gz",fq):
 		zipped = '--gunzip'

	if snp_tolerance:
		output_prefix = fq.split("/")[-1].split(".")[0] + "_SNP"
		map_cmd = "gsnap " + zipped + " -D " + genome_index_path + " -d " + genome_index_name + " -V " + snp_index_path + " -v " \
		+ snp_index_name + " -t " + str(threads) + " --split-output " + out_dir + output_prefix + \
		" --format=sam --genome-unk-mismatch=0 --md-lowercase-snp  " + \
		fq + " &>> " + out_dir + "align.log"
	else:
		output_prefix = fq.split("/")[-1].split(".")[0] + "_noSNP"
		map_cmd = "gsnap " + zipped + " -D " + genome_index_path + " -d " + genome_index_name + " -t " + str(threads) + \
		" --split-output " + out_dir + output_prefix + " --format=sam --genome-unk-mismatch=0 --md-lowercase-snp " + \
		fq + " &>> " + out_dir + "align.log"

	log.info("Aligning reads to {}...".format(genome_index_name))
	subprocess.call(map_cmd, shell = True)
	
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
			if re.search("mult",file.name):
				cmd = "samtools view -@ " + str(threads) + " -F 0x904 -c " + out_dir + file.name
				multi_count = int(subprocess.check_output(cmd, shell = True))
				if keep_temp:
					cmd = "samtools view -@ " + str(threads) + " -bh -o " + out_dir + file.name + ".bam " + out_dir + file.name
					subprocess.call(cmd, shell = True)
					os.remove(out_dir + file.name)
				elif not keep_temp:
					os.remove(out_dir + file.name)
			elif re.search("uniq",file.name):
				cmd = "samtools view -@ " + str(threads) + " -c " + out_dir + file.name
				unique_count = int(subprocess.check_output(cmd, shell = True))
				unique_bam = out_dir + file.name + ".bam"
				cmd = "samtools view -@ " + str(threads) + " -bh " + out_dir + file.name + " | samtools sort -@ " + str(threads) + " -o " + unique_bam + " - "
				subprocess.call(cmd, shell = True)
				os.remove(out_dir + file.name)
			elif re.search("nomapping",file.name):
				cmd = "samtools view -@ " + str(threads) + " -c " + out_dir + file.name
				unmapped_count = int(subprocess.check_output(cmd, shell = True))
				if keep_temp:
					cmd = "samtools view -@ " + str(threads) + " -bh -o " + out_dir + file.name + ".bam " + out_dir + file.name
					subprocess.call(cmd, shell = True)
					os.remove(out_dir + file.name)
				elif not keep_temp:
					os.remove(out_dir + file.name)

		total_count = unique_count + multi_count + unmapped_count

		stats_out.write("{}\nUniqely mapped reads: {:d} ({:.0%}) \nMulti-mapping reads: {:d} ({:.0%}) \nUnmapped reads: {:d} ({:.0%}) \nTotal: {:d}\n\n"\
			.format(fq.split("/")[-1], unique_count, (unique_count/total_count),multi_count, (multi_count/total_count), unmapped_count, (unmapped_count/total_count), total_count))
	
	return(unique_bam, total_count)

def countReads(unique_bam_list, mode, threads, out_dir):

	log.info('\n+-----------------------------------+\
		\n| Counting reads with featureCounts |\
		\n+-----------------------------------+')

	# Find tRNA gff produced earlier...
	gff_file = glob.glob(out_dir + '*.gff')[0]

	cmd = "featureCounts -T " + str(threads) 

	if mode == "none":
		cmd += " -a " + gff_file + " -o " + out_dir + "counts.txt " + " " .join(unique_bam_list) + " &>> " + out_dir + "featureCounts.log"
	elif mode == "all":
		cmd += " -O -a " + gff_file + " -o " + out_dir + "counts.txt " + " " .join(unique_bam_list) + " &>> " + out_dir + "featureCounts.log"
	elif mode == "fraction":
		cmd += " -O --fraction -a " + gff_file + " -o " + out_dir + "counts.txt " + " " .join(unique_bam_list) + " &>> " + out_dir + "featureCounts.log"

	subprocess.call(cmd, shell=True)

	log.info("Read counts per tRNA/cluster saved to " + out_dir + "counts/counts.txt")

if __name__ == '__main__':
	mainAlign(sys.argv[1:])

