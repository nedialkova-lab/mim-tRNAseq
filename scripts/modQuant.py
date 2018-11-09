#! /usr/bin/env python3

######################################################################################
# Analysis of modifications per tRNA position based on mismatches in bam files #
######################################################################################

import os, logging
import re
import pysam
from getCoverage import getBamList
from multiprocessing import Pool
import pandas as pd
import numpy as np
from collections import defaultdict

log = logging.getLogger(__name__)
modTable_list = list()

def countMods_mp(out_dir, inputs, cov_table):
# modification counting and table generation
	
	# generate DataFrame with the references (as index) and positions imported from cov_table
	modTable = defaultdict(lambda: defaultdict(int))
	bam_file = pysam.AlignmentFile(inputs, "rb")
	log.info('Analysing misincorporation at modified nucleotides for {}...'.format(inputs))
	for read in bam_file.fetch(until_eof=True):
		query = read.query_name
		reference = read.reference_name
		# get MD tags for mismatches, split into list of integers and characters
		md_tag = read.get_tag('MD')
		md_list = re.split('(.*?)([A-Z])', md_tag)
		md_list = list(filter(None, md_list))
		# get offset of read mapping position to reference start in order to adjust mismatch position - in the below code offset is also used
		# to keep track of the mismatch position relative to the reference
		# 20 is used to account for 20 upstream Ns in references
		offset = 20 - read.reference_start
		for interval in md_list:
			if interval.isdigit(): #stretch of matches
				interval = int(interval) - offset
				offset += int(interval)
			else: # is a mismatch
				offset += 1
				modTable[reference][offset] += 1

	log.info('Analysis complete for {}...'.format(inputs))

def generateModsTable(sampleGroups, out_dir, mod_lists, max_multi, cov_table):
# Wrapper function to call countMods_mp with multiprocessing

	log.info("\n+---------------------------------------------+\
	\n| Analysing misincorporations and stops to RT |\
	\n+---------------------------------------------+")

	# Get bam info using function in getCoverage
	baminfo, bamlist = getBamList(sampleGroups)

	if max_multi > len(bamlist):
		max_multi = len(bamlist)

	# initiate multiprocessing pool and run with bam names
	pool = Pool(max_multi)
	func = partial(countMods_mp, out_dir, bamlist)
	pool.map(func, cov_table)
	pool.close()
	pool.join()

	for bam, info in baminfo.items():
			