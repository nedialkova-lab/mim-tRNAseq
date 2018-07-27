#! /usr/bin/env python3

import subprocess
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from pybedtools import BedTool

def getBamList (sampleGroups):
# reads sampleGroups file and creates dictionary of bam and groups
# sampleGroups text file contains bam file locations, group and library size in read number, tab-separated. 

	sampleGroups = open(sampleGroups,"r")
	bamfiles = defaultdict(list)
	for line in sampleGroups:
		line = line.strip()
		currbam = str(line.split("\t")[0])
		condition = line.split("\t")[1]
		librarySize = int(line.split("\t")[2])
		bamfiles[currbam] = [condition,librarySize]

	return(bamfiles)
	sampleGroups.close()

def getCoverage(tRNAbed, sampleGroups, out_dir):
# Uses bedtools coverage and pandas generate coverage in 5% intervals per gene and isoacceptor for plotting

	print("+-----------------------------------+\
		\n| Calculating coverage and plotting |\
		\n+-----------------------------------+\n")

	bamfiles = getBamList(sampleGroups)
	cov_mean = list()

	for bam, info in bamfiles.items():
		a = BedTool(tRNAbed)
		b = BedTool(bam)
		print("Running bedtools coverage on {}...".format(bam))
		coverage_out = a.coverage(b, d = True, s = True).saveas(out_dir + "coverage.txt")
		#cmd = "bedtools coverage -a " + tRNAbed + " -b " + bam + " -d -s > " + out_dir + "coverage.txt"
		#subprocess.call(cmd, shell = True)

		coverage = pd.read_table(out_dir + "coverage.txt", header = None, index_col = 0)[[6,7]]
		coverage['aa'] = coverage.index.format()
		coverage['aa'] = coverage['aa'].str.split("-").str[-4]
		coverage.columns = ['pos','cov','aa']
		coverage['condition'] = info[0]
		coverage['cov_norm'] = coverage['cov'] / info[1]

		# bin by grouping by each gene and cutting into 20 - i.e. 5% bins by gene length (pos)
		coverage['bin'] = coverage.groupby([coverage.index])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
		# append big coverage table
		cov_mean.append(coverage) 

	# concatenate all tables together, groupby + mean
	cov_mean = 	pd.concat(cov_mean, axis = 0)
	cov_mean_gene = cov_mean.groupby([cov_mean.index, 'bin', 'condition']).mean()
	cov_mean_gene.to_csv(out_dir + "coverage_bygene.txt", sep = "\t")
	cov_mean_aa = cov_mean.groupby(['aa', 'bin', 'condition']).mean()	
	cov_mean_aa.to_csv(out_dir + "coverage_byaa.txt", sep = "\t")

	os.remove(out_dir + "coverage.txt")	

def plotCoverage(out_dir):
	script_path = os.path.dirname(os.path.realpath(__file__))
	command = "Rscript " + script_path + "/coveragePlot.R " + out_dir + "coverage_bygene.txt " + out_dir + "coverage_byaa.txt " + out_dir 
	subprocess.call(command, shell = True)

#getCoverage("eColi_CCAvsCC/eColi_tK_clusters.bed", "eColi_CCAvsCC/sampleData_cov.txt","eColi_CCAvsCC/")
#plotCoverage("eColi_CCAvsCC/")
