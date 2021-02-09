#! /usr/bin/env python3

from __future__ import absolute_import
import logging
from collections import defaultdict
from itertools import chain, combinations as comb
from .ssAlign import aligntRNA
from .getCoverage import getBamList
import re
import pybedtools
from multiprocessing import Pool
from functools import partial
import pickle

####################################################################################################################
# Determine sets of distinguishing mismatches and insertions in cluster members to split reads to unique sequences #
####################################################################################################################

log = logging.getLogger(__name__)

def dd_set():
	return(defaultdict(set))

def dd():
	return(defaultdict(list))

def GetPowerset(s):
#Get the powerset of a given set (all subsets incl. empty and full set).
    
    return list(chain(*map(lambda x: comb(s, x), range(0, len(s)+1))))

def reformatInDelDict (dict, newDict, type):
# Code to reformat insert and delete dictionary information for processing to find unique sites

    for cluster, data in dict.items():
        for pos, members in data.items():
            for member in members:
                isodecoder = "-".join(member.split("-")[:4]) if not "chr" in member else member
                posIdentity = str(pos) + type
                newDict[cluster][isodecoder].add(posIdentity)

    return(newDict)

def updateMismatchPosDict (mismatchPosDict, inputDict):
    for cluster, data in inputDict.items():
        for isodecoder, posInfo in data.items():
            mismatchPosDict[cluster][isodecoder].update(posInfo)

    return(mismatchPosDict)

def atoi(text):
# string to int conversion of nuber component of mismatch
    return int(text) if text.isdigit() else text

def natural_keys(text):
# return keys to sort alphanumeric combinations present in mismatches according to human sorting (i.e numerically)
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def natural_keys_list(text):
# same as above but to sort tuples within a list of tuples
    l = list()
    for a in text:
        for c in re.split(r'(\d+)', a):
            l.append(atoi(c))
    
    return(l)

def findUniqueSubset (inputDict, outputDict):
# For dictionary of mismatches, insertions and deletions, find unique minimal distinguishing subset of positions and update outputDict

    for cluster, data in inputDict.items():
        powersets  = [GetPowerset(data[s]) for s in data]
        # temp dictionary of all isodecoder unique sets in the case that a subsequenct isodecoder unique set matches a previous one
        temp_uss_dict = defaultdict(list)
        # same for all total set of mismatches in the case that no subsets can be found for previous isodecoders that are unique
        temp_fullset_dict = defaultdict(list)
        for i in range(0, len(data)):
            # declare isodecoder name and powerset for current loop iteration (isodecoder)
            isodecoder, mismatch_powerset = list(data)[i], powersets[i]
            temp_fullset_dict[isodecoder] = tuple(sorted((max(mismatch_powerset, key=len))))
            # list the other powersets to compare against
            ops = [x for ind, x in enumerate(powersets) if ind != i]
            # find unique subsets: those that are only a subset of the current set
            # and not found in the powerset of any of the other sets
            uss = list(set(mismatch_powerset)-set([x for y in ops for x in y if x in mismatch_powerset]))
            # if the unique set is empty it means the the full set of mismatches is also found as some subset in another isodecoder
            # this is ok because this other isodecoder will have this full set subtracted - therefore use the full set as the unique set of mismatches
            if not uss:
                uss = [(max(mismatch_powerset, key=len))]
            # filter empty tuples
            uss = list(filter(None, uss))
            # sort descending within whole list
            uss = sorted(uss, key = natural_keys_list, reverse = True)
            # sort within individual subsets 
            for i, s in enumerate(uss):
                uss[i] = tuple(sorted(s, key=natural_keys))
            temp_uss_dict[isodecoder] = uss # add full uss set to temp dict
            # choose the shortest/minimal unique subset
            min_uss = tuple(sorted(min(uss, key = len)))
            index = 0
            while min_uss in outputDict[cluster].keys(): # if min_uss already in outputDict
                index += 1
                try:
                    min_uss = uss[index] # fetch next uss
                except IndexError: # if there is only one uss, change the uss of the conflicting isodecoder
                    conflictingIso = outputDict[cluster][min_uss][0]
                    conflictingIso_uss = None
                    iso_ss = temp_uss_dict[conflictingIso]
                    for ss in iso_ss: # find new unique uss for conflicting isodecoder
                        if not ss in outputDict[cluster].keys():
                            conflictingIso_uss = ss
                            outputDict[cluster][conflictingIso_uss].append(conflictingIso)
                            outputDict[cluster][conflictingIso_uss].append(temp_fullset_dict[conflictingIso])
                            del outputDict[cluster][min_uss]
                            break
                    if not conflictingIso_uss: # if no uss can be found, use the full set of mismatches
                        conflictingIso_uss = temp_fullset_dict[isodecoder]
                        outputDict[cluster][conflictingIso_uss].append(conflictingIso)
                        outputDict[cluster][conflictingIso_uss].append(conflictingIso_uss)
                        del outputDict[cluster][min_uss]
        
            outputDict[cluster][min_uss].append(isodecoder)
            outputDict[cluster][min_uss].append(temp_fullset_dict[isodecoder])

    # add "-1" so that referencing this isodecoder by its full gene name elsewhere in code is not an issue
    for cluster, data in outputDict.items():
        for uss, isodecoder in data.items():
            newIso = isodecoder[0] + "-1" if not "chr" in isodecoder[0] else isodecoder[0]
            newData = [newIso, isodecoder[1]]
            outputDict[cluster][uss] = newData

    return(outputDict)

def splitIsodecoder(cluster_perPos_mismatchMembers, insert_dict, del_dict, tRNA_dict, cluster_dict, out_dir, experiment_name):
# Determine minimal set of most 3' mismatches and/or insertions that characterise an isodecoder

    log.info("\n+------------------------------------------------------------------------------+\
        \n| Characterizing cluster mismatches for read splitting by unique tRNA sequence |\
       \n+------------------------------------------------------------------------------+")

    log.info("** Assessing mismatches between cluster members and parent... **")

    # Reformat cluster_perPos_mismatchMembers to contain cluster -> isodeocoder -> set of mismatches and identity
    cluster_MemberMismatchPos = defaultdict(dd_set)

    for cluster, data in cluster_perPos_mismatchMembers.items():
        for pos, members in data.items():
            for member in members:
                member_seq = tRNA_dict[member]['sequence']
                # find number of inserts before mismatch in question to ensure that the correct identity in the member is sliced by subtracting from the mismatch pos in the parent
                ins_num = len(set([ins for ins in insert_dict[cluster] if member in insert_dict[cluster][ins] and ins < pos]))
                # find number of deletions before mismatch in question to ensure that the correct identity in the member is sliced by subtracting from the mismatch pos in the parent
                del_num = len(set([deletion for deletion in del_dict[cluster] if member in del_dict[cluster][deletion] and deletion < pos]))
                identity = member_seq[pos-ins_num+del_num]

                isodecoder = "-".join(member.split("-")[:4]) if not "chr" in member else member
                posIdentity = str(pos) + identity
                cluster_MemberMismatchPos[cluster][isodecoder].add(posIdentity)
    
    # Reformat insert_dict similarly as above, and update cluster_MemberMismatchPos
    cluster_MemberInsertPos = defaultdict(dd_set)
    cluster_MemberInsertPos = reformatInDelDict(insert_dict, cluster_MemberInsertPos, "Ins")
    cluster_MemberMismatchPos =updateMismatchPosDict(cluster_MemberMismatchPos, cluster_MemberInsertPos)

    # And again for del_dict similarly as above
    cluster_MemberDeletePos = defaultdict(dd_set)
    cluster_MemberDeletePos = reformatInDelDict(del_dict, cluster_MemberDeletePos, "Del")
    cluster_MemberMismatchPos = updateMismatchPosDict(cluster_MemberMismatchPos, cluster_MemberDeletePos)
    
    # Build nested dictionary of unique minimal set of mismatches and insertions that distinguish an isodecoder from parent and all others in cluster
    unique_isodecoderMMs = defaultdict(dd)
    unique_isodecoderMMs = findUniqueSubset(cluster_MemberMismatchPos, unique_isodecoderMMs)
    isodecoder_sizes = defaultdict(int)

    # Check that all unique sequences can be deconvoluted
    # count clusters in cluster_dict that are composed of only one sequence and those that are composed of multiple sequences
    singleSeq_clusters = 0
    multiSeq_size = 0
    multiSeq_names = set()
    for cluster, members in cluster_dict.items():
        child_iso_set = {tRNA_dict[member]['sequence'].upper() for member in members}
        if len(child_iso_set) == 1:
            singleSeq_clusters += 1
            isodecoder_sizes[cluster] = 1
        else:
            multiSeq_size += len(child_iso_set)
            multiSeq_names.update({"-".join(member.split("-")[0:4]) for member in members} - {"-".join(cluster.split("-")[0:4])})

    # count deconvoluted sequences above
    deconv_sequences_num = 0
    deconv_names = set()
    for cluster, data in unique_isodecoderMMs.items():
        # count isodecoder sizes
        isodecoder_sizes[cluster] = len([info['sequence'].upper() for tRNA, info in tRNA_dict.items() if info['sequence'].upper() == tRNA_dict[cluster]['sequence'].upper()])
        for isodecoder in data.values():
            isodecoder_sizes[isodecoder[0]] = len([info['sequence'].upper() for tRNA, info in tRNA_dict.items() if info['sequence'].upper() == tRNA_dict[isodecoder[0]]['sequence'].upper()])
        deconv_sequences_num += 1
        unique_deconv = {member[0] for member in data.values()}
        deconv_sequences_num += len(unique_deconv)
        deconv_names.update({"-".join(member[0].split("-")[0:4]) for member in data.values()})

    log.info("Total unique sequences: {}".format(singleSeq_clusters + multiSeq_size)) # Total unique sequences 
    log.info("Single sequence clusters: {}".format(singleSeq_clusters))
    log.info("Deconvoluted sequences: {}".format(deconv_sequences_num))
    log.info("Total deconvoluted sequences {}".format(singleSeq_clusters + deconv_sequences_num))

    nondeconv_isos = list(multiSeq_names - deconv_names)
    # print warnings if not all clusters deconvoluted
    if (singleSeq_clusters + multiSeq_size) != (singleSeq_clusters + deconv_sequences_num):
        log.warning("*** Unique mismatches and indels not found for some unique tRNA sequences in input!!***\n \
Some clusters not fully deconvoluted. There can be various reasons for this, including errors in input sequences\n \
The following isodecoders were not distinguished uniquely from the cluster parent:")
        log.warning(nondeconv_isos)
        log.warning("Continuing with analysis - clusters with un-deconvoluted tRNAs will be excluded from differential expression analysis.")

    # nondeconv_isos above equivalent to old splitBool variable from earlier code. Copy and return for filtering these from DESeq2 analysis
    splitBool = nondeconv_isos

    # save isodecoder info
    with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
        isodecoderInfo.write("Isodecoder\tsize\n")
        for isodecoder, size in isodecoder_sizes.items():
            isodecoder = "-".join(isodecoder.split("-")[:-1]) if not "chr" in isodecoder else isodecoder
            isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	# write isodecoder fasta for alignment and context analysis
    with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
        for seq in isodecoder_sizes.keys():
            shortname = "-".join(seq.split("-")[:-1]) if not "chr" in seq else seq
            tempSeqs.write(">" + shortname + "\n" + tRNA_dict[seq]['sequence'] + "\n")
    aligntRNA(tempSeqs.name, out_dir)

    # debugging pickle output
    #with open("mismatch_dict.pkl","wb") as mismatch_dict_out, open("insert_dict.pkl","wb") as insert_dict_out, open("del_dict.pkl","wb") as del_dict_out, open("cluster_dict.pkl","wb") as cluster_dict_out, open("unique_mms.pkl", "wb") as unique_mms_out:
    #    pickle.dump(cluster_MemberMismatchPos, mismatch_dict_out)
    #    pickle.dump(insert_dict, insert_dict_out)
    #    pickle.dump(del_dict, del_dict_out)
    #    pickle.dump(cluster_dict, cluster_dict_out)
    #    pickle.dump(unique_isodecoderMMs, unique_mms_out) 

    return(unique_isodecoderMMs, splitBool, isodecoder_sizes)

def covCheck_mp(coverageBed, unique_isodecoderMMs, covDiff, input):
    # get positional coverage per cluster per bam file and check if 3':5' coverage is greater than covDiff
    unsplit = set()
    unsplit_isosOnly = set()
    log.info("Calculating nucleotide coverage for {}".format(input))
    bam = pybedtools.BedTool(input)
    bed = pybedtools.BedTool(coverageBed)
    cov = bed.coverage(bam, s = True, d = True)
    cov_df = cov.to_dataframe()

    # check coverage diff for each unique sequence in each cluster
    for cluster in unique_isodecoderMMs.keys():
        for mismatch in unique_isodecoderMMs[cluster]:
            # split positions from identities in distinguishing mismatches and find most 5' (i.e. min)
            mismatchNums = [int(re.search('(\d+\.|\d+)+', pos).group(0)) for pos in mismatch]
            minMismatch = min(mismatchNums) + 1 # 0 to 1 based numbering for cov_df thickStart

            # check coverage at position/coverage at 3' is >= covDiff: thickStart is pos, thickEnd is coverage
            # determine most 3' pos of current cluster
            end = max(cov_df.loc[cov_df.name == cluster,'thickStart'])
            if minMismatch < end:
                endCov = int(cov_df.loc[(cov_df.name == cluster) & (cov_df.thickStart == end - 5), 'thickEnd'])
                # use this end - 5 (to exlude variability at 3'-CCA coverage) as the 3' position to measure against
                ratio = cov_df.loc[(cov_df.name == cluster) & (cov_df.thickStart == minMismatch), 'thickEnd'].astype(float) / endCov
            #try:
                if float(ratio) < covDiff:
                    unsplit.add(cluster)
                    unsplit.add(unique_isodecoderMMs[cluster][mismatch][0])
                    unsplit_isosOnly.add(unique_isodecoderMMs[cluster][mismatch][0])
            #except TypeError:
            #   print(cov_df.head())
            #   print(mismatch)
            #    print(minMismatch)
            #    print(endCov)
            #    print(ratio)

    return(unsplit, unsplit_isosOnly)

def unsplitClusters(coverageData, coverageBed, unique_isodecoderMMs, threads, covDiff = 0.5):
    # Define unique sequences not able to be split based on significant drop in coverage before distinguishing mismatch
   
    log.info("\n+---------------------------------------------------------------------------------+\
        \n| Determining un-deconvoluted clusters due to insufficient coverage at mismatches |\
        \n+---------------------------------------------------------------------------------+")

    log.info("*** Un-deconvoluted sequences being defined based on coverage difference more than {}".format(covDiff))
    baminfo, bamlist = getBamList(coverageData)
    
    if len(baminfo) > threads:
	    multi = threads
    else:
	    multi = len(baminfo)

    # initiate custom non-daemonic multiprocessing pool and run with bam names
    log.info("Determining unsplittable sequences...")
    pool = Pool(multi)
    func = partial(covCheck_mp, coverageBed, unique_isodecoderMMs, covDiff)
    unsplit, unsplit_isosOnly = zip(*pool.map(func, bamlist))
    pool.close()
    pool.join()

    unsplit_all = list(set().union(*unsplit))
    unsplit_isosAll = len(list(set().union(*unsplit_isosOnly)))
    log.info("{} unique sequences excluded from deconvolution due to reductions in required coverage at mismatches of more than {:.2%}".format(unsplit_isosAll, covDiff))

    return(unsplit_all)

def writeIsodecoderTranscripts(out_dir, experiment_name, cluster_dict, tRNA_dict):
	# write isodecoderTransripts.fa when cluster_id == 1 to avoid issues with shortened isodecoder names and output files

	with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
		for seq in cluster_dict.keys():
			shortname = "-".join(seq.split("-")[:-1]) if not "chr" in seq else seq
			tempSeqs.write(">" + shortname + "\n" + tRNA_dict[seq]['sequence'] + "\n")
	aligntRNA(tempSeqs.name, out_dir)

def getIsodecoderSizes(out_dir, experiment_name, tRNAdict):
	# get isodecoder sizes for tRNA sequences - useful for when clustering is disabled and above function is not applicable

	isodecoder_sizes = defaultdict(int)
	already_added = set()
	for tRNA in tRNAdict:
		if tRNA not in already_added:
			sameSeq = [tRNAs for tRNAs, data in tRNAdict.items() if data['sequence'] == tRNAdict[tRNA]['sequence']]
			already_added.update(sameSeq)
			isodecoder_sizes[tRNA] = len(sameSeq)

	with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
		isodecoderInfo.write("Isodecoder\tsize\n")
		for isodecoder, size in isodecoder_sizes.items():
			isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	return(isodecoder_sizes)
