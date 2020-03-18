#! /usr/bin/env python3

####################################################
# 3' CCA analysis of aligned reads                 #
# plotting of proportions of 3' endings            #
# DESeq2 analysis on CC vs CCA counts for clusters #
####################################################

import os, subprocess
import logging

log = logging.getLogger(__name__)

def plotDinuc(out_dir):

	log.info("\n+-----------------+\
		\n| 3'-CCA analysis |\
		\n+-----------------+")

	out = out_dir + "CCAanalysis/"
	script_path = os.path.dirname(os.path.realpath(__file__))
	command = ["Rscript", script_path + "/ccaPlots.R", out + "AlignedDinucProportions.csv", out + "/CCAcounts.csv", out]
	subprocess.check_call(command)

	log.info("CCA analysis done and plots created. Located in {}".format(out))