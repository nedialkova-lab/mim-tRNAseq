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
from multiprocessing import Pool

log = logging.getLogger(__name__)

def filterCoverage (cov_table, min_cov):
# returns isodecoders as list from counts table with less than min_cov reads (excluding mito clusters)
	
	filtered_list = list(cov_table[(cov_table.values < min_cov).any(1) & (~cov_table.index.str.contains('mito'))].index)

	log.info("{} clusters filtered out according to minimum coverage threshold: {}".format(len(filtered_list), min_cov))

	return(filtered_list)

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

def getCoverage(sampleGroups, out_dir, min_cov, control_cond, filtered_cov):
# Uses bedtools coverage and pandas generate coverage in 5% intervals per gene and isoacceptor for plotting

	log.info("\n+-----------------------------------+\
		\n| Calculating coverage and plotting |\
		\n+-----------------------------------+")

	baminfo, bamlist = getBamList(sampleGroups)
	cov_mean = list()

	for bam, info in baminfo.items():

		coverage = pd.read_csv(out_dir + bam.split("/")[-1] + "_coverage.txt", index_col = 0, sep = "\t") 
		coverage['aa'] = coverage.index.format()
		coverage.loc[coverage.aa.str.contains('mito'), 'aa'] = "mito" + coverage[coverage.aa.str.contains('mito')].aa.str.split("-").str[-4]
		coverage.loc[coverage.aa.str.contains('nmt'), 'aa'] = "nmt" + coverage[coverage.aa.str.contains('nmt')].aa.str.split("-").str[-4]
		coverage.loc[~coverage.aa.str.contains('mito') & ~coverage.aa.str.contains('nmt'), 'aa'] = coverage[~coverage.aa.str.contains('mito') & ~coverage.aa.str.contains('nmt')].aa.str.split("-").str[-4]
		coverage = coverage[['pos','cov','aa','bam']]
		coverage['condition'] = info[0]
		coverage['cov'] = coverage['cov'].astype(float)
		coverage['cov_norm'] = coverage['cov'] / info[1]

		# bin by grouping by each gene and cutting into 20 - i.e. 5% bins by gene length (pos)
		coverage['bin'] = coverage.groupby([coverage.index])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
		# append big coverage table and remove original coverage output
		cov_mean.append(coverage)
		#os.remove(out_dir + bam.split("/")[-1] + "_coverage.txt")	 

	# concatenate all tables together, groupby + mean
	cov_mean = pd.concat(cov_mean, axis = 0)
	cov_mean = cov_mean[~cov_mean.index.isin(filtered_cov)]
	cov_mean_gene = cov_mean.copy()
	cov_mean_gene['Cluster'] = cov_mean_gene.index.format()
	cov_mean_gene.loc[cov_mean_gene.Cluster.str.contains("mito"), "Cluster"] = "mito" + cov_mean_gene[cov_mean_gene.Cluster.str.contains("mito")].Cluster.str.split("-").str[1:].str.join('-')
	cov_mean_gene.loc[cov_mean_gene.Cluster.str.contains("nmt"), "Cluster"] = "nmt" + cov_mean_gene[cov_mean_gene.Cluster.str.contains("nmt")].Cluster.str.split("-").str[1:].str.join('-')
	cov_mean_gene.loc[~cov_mean_gene.Cluster.str.contains("mito") & ~cov_mean_gene.Cluster.str.contains("nmt"), "Cluster"] = cov_mean_gene[~cov_mean_gene.Cluster.str.contains("mito") & ~cov_mean_gene.Cluster.str.contains("nmt")].Cluster.str.split("-").str[1:].str.join('-')
	cov_mean_gene = cov_mean_gene[['Cluster','pos','cov','aa','condition','cov_norm','bin','bam']]
	cov_mean_gene = cov_mean_gene.groupby(['Cluster', 'bin', 'condition', 'bam']).mean()
	cov_mean_gene = cov_mean_gene.dropna()
	cov_mean_gene.to_csv(out_dir + "coverage_bygene.txt", sep = "\t")

	# coverage per amino acid
	cov_mean_aa = cov_mean.groupby(['aa', 'condition', 'bam', 'pos']).sum() # sum coverages for all clusters for each amino acid
	cov_mean_aa = cov_mean_aa.reset_index()
	# remove last pos for each group (i.e. max pos) since this has very low coverage arising from members with different lengths
	cov_mean_aa = cov_mean_aa.groupby('aa').apply(lambda group: group.loc[group['pos'] != group['pos'].max()])
	cov_mean_aa = cov_mean_aa.drop(columns='aa')
	cov_mean_aa = cov_mean_aa.reset_index()
	cov_mean_aa = cov_mean_aa.drop(columns='level_1')
	cov_mean_aa['bin'] = cov_mean_aa.groupby(['aa','condition','bam'])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
	cov_mean_aa = cov_mean_aa.groupby(['aa', 'bin', 'condition', 'bam']).mean()
	cov_mean_aa = cov_mean_aa.dropna()
	cov_mean_aa	= cov_mean_aa.reset_index()
	cov_mean_aa.to_csv(out_dir + "coverage_byaa.txt", sep = "\t", index = False)
	
	# 5' to 3' coverage ratio calculation for ordering AAs on coverage plot
	cov_mean_aa_controlcond = cov_mean_aa[cov_mean_aa.condition == control_cond]
	bam = pd.unique(cov_mean_aa_controlcond['bam'])[0] 
	cov_mean_aa_controlcond = cov_mean_aa_controlcond[cov_mean_aa_controlcond.bam == bam]
	cov_ratios = dict()
	for aa, data in cov_mean_aa_controlcond.groupby('aa'):
		try:
			ratio = float(data[data.bin == 8]['cov_norm']) / float(data[data.bin == 92]['cov_norm'])
		except ZeroDivisionError:
			ratio = float(data[data.bin == 92]['cov_norm'])
		cov_ratios[aa] = ratio
	sorted_aa = sorted(cov_ratios, key = cov_ratios.get)
	sorted_aa = "_".join(str(e) for e in sorted_aa)

	return(sorted_aa)

def plotCoverage(out_dir, mito_trnas, sorted_aa):
	
	script_path = os.path.dirname(os.path.realpath(__file__))
	command = ["Rscript", script_path + "/coveragePlot.R", out_dir + "coverage_bygene.txt", out_dir + "coverage_byaa.txt", out_dir, sorted_aa, mito_trnas]
	try:
	#with StreamLogger.StreamLogger(logging.INFO) as out:
		subprocess.check_call(command)
	except Exception as e:
		logging.error("Error in {}".format(command), exc_info=e)
		raise
