#! /usr/bin/env python3

######################################################
# Coverage calculation and table output for plotting #
######################################################

import subprocess
import pandas as pd
import numpy as np
import os, logging
from functools import partial
from collections import defaultdict
from pybedtools import BedTool
from multiprocessing import Pool

log = logging.getLogger(__name__)

def bedtools_mp (tRNAbed, out_dir, inputs):
# runs pybedtools coverage for input data
# useful for multiprocessing of coverage tasks
	
	a = BedTool(tRNAbed)
	b = BedTool(inputs)
	log.info("Running bedtools coverage on {}...".format(inputs))
	cov = a.coverage(b, d = True, s = True).saveas(out_dir + inputs.split("/")[-1] + "_coverage.txt")
	log.info("Coverage calculation complete for {}".format(inputs))	

def getBamList (sampleGroups):
# reads sampleGroups file and creates dictionary of bam and groups
# sampleGroups text file contains bam file locations, group and library size in read number, tab-separated. 

	sampleGroups = open(sampleGroups,"r")
	baminfo = defaultdict(list)
	bamlist = list()
	for line in sampleGroups:
		line = line.strip()
		currbam = str(line.split("\t")[0])
		condition = line.split("\t")[1]
		librarySize = int(line.split("\t")[2])
		baminfo[currbam] = [condition,librarySize]
		bamlist.append(currbam)

	return(baminfo, bamlist)
	sampleGroups.close()

def getCoverage(tRNAbed, sampleGroups, out_dir, max_multi):
# Uses bedtools coverage and pandas generate coverage in 5% intervals per gene and isoacceptor for plotting

	log.info("\n+-----------------------------------+\
		\n| Calculating coverage and plotting |\
		\n+-----------------------------------+")

	baminfo, bamlist = getBamList(sampleGroups)
	cov_mean = list()

	# multiprocessing of bedtools coverage

	if max_multi > len(bamlist):
		max_multi = len(bamlist)

	# create process pool, make partial function of bedtools_mp (above) with some inputs, 
	# map partial function with bam inputs to pool of prcoesses for multithreaded coverage calculation
	pool = Pool(max_multi)
	func = partial(bedtools_mp, tRNAbed, out_dir)
	pool.map(func, bamlist)
	pool.close()
	pool.join()

	for bam, info in baminfo.items():

		coverage = pd.read_table(out_dir + bam.split("/")[-1] + "_coverage.txt", header = None, index_col = 0)[[6,7]]
		coverage['aa'] = coverage.index.format()
		coverage['aa'] = coverage['aa'].str.split("-").str[-4]
		coverage.columns = ['pos','cov','aa']
		coverage['condition'] = info[0]
		coverage['cov_norm'] = coverage['cov'] / info[1]

		# bin by grouping by each gene and cutting into 20 - i.e. 5% bins by gene length (pos)
		coverage['bin'] = coverage.groupby([coverage.index])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
		# append big coverage table and remove original coverage output
		cov_mean.append(coverage)
		os.remove(out_dir + bam.split("/")[-1] + "_coverage.txt")	 

	# concatenate all tables together, groupby + mean
	cov_mean = 	pd.concat(cov_mean, axis = 0)
	cov_mean_gene = cov_mean.groupby([cov_mean.index, 'bin', 'condition']).mean()
	cov_mean_gene.to_csv(out_dir + "coverage_bygene.txt", sep = "\t")
	cov_mean_aa = cov_mean.groupby(['aa', 'bin', 'condition']).mean()	
	cov_mean_aa.to_csv(out_dir + "coverage_byaa.txt", sep = "\t")

	return(cov_mean)

def plotCoverage(out_dir):
	script_path = os.path.dirname(os.path.realpath(__file__))
	command = "Rscript " + script_path + "/coveragePlot.R " + out_dir + "coverage_bygene.txt " + out_dir + "coverage_byaa.txt " + out_dir 
	subprocess.call(command, shell = True)

