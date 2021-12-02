#! /usr/bin/env python3

from __future__ import absolute_import
import logging
from collections import defaultdict
from itertools import chain, combinations as comb

from numpy.core.fromnumeric import shape
from .ssAlign import aligntRNA, tRNAclassifier
from .getCoverage import getBamList
import re, copy
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
# Code to reformat insertion and deletion dictionary information for processing to find unique sites

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

def findUniqueSubset (inputDict, outputDict, tRNA_dict):
# For dictionary of mismatches, insertions and deletions, find unique minimal distinguishing subset of positions and update outputDict

    # define canonical position dictionary to cross check mismatches with canonical modified sites
    tRNA_ungap2canon = tRNAclassifier()[1]
    mods = {'9':['G', 'A'], 
            '20':['C', 'T'], 
            '20a':['T'], 
            '20b':['T'], 
            '26':['G'], 
            '32':['C'], 
            '37':['G', 'A'], 
            '58':['A']}

    # clusters and isodecoders not split due to canonical modified positions in all unique mismatch subsets
    notSplit_mods_clusterInfo = defaultdict(list)

    for cluster, data in inputDict.items():
        powersets  = [GetPowerset(data[s]) for s in data]
        # temp dictionary of all isodecoder unique sets in the case that a subsequent isodecoder unique set matches a previous one
        temp_uss_dict = defaultdict(list)
        # all full sets of mismatches in the case that no subsets can be found for previous isodecoders that are unique
        temp_fullset_dict = defaultdict(list)
        
        # filtered list of mods in cases where sites aren't present (e.g. 20a or b)
        modsFilter = {i:v for i, v in mods.items() if i in tRNA_ungap2canon[cluster].values()}
        # ungapped parent transcript locations of filtered canonical mod sites
        clusterMods_pos = {list(tRNA_ungap2canon[cluster].keys())[list(tRNA_ungap2canon[cluster].values()).index(i)]:v for i, v in modsFilter.items()}
        # concat. pos + nucleotide for each element of clusterMods_pos to match to uss tuples below
        clusterMods_list = [str(i) + str(x) for i, nucs in clusterMods_pos.items() for x in nucs]

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

            # sort within individual subsets and remove subsets that contain canonical modified positions
            for i, s in enumerate(uss):
                uss[i] = tuple(sorted(s, key=natural_keys))
                ussPosOnly = [int(re.search("([0-9]*)[A-Z]?", x).group(1)) for x in uss[i] if not re.search("Ins|Del", x)] # extract positions only from uss tuple (exclude Ins Del positions)
                ussClusterParent = [str(x) + tRNA_dict[cluster]['sequence'][x] for x in ussPosOnly] # build similar tuple to uss[i] except with cluster parent nucleotide identity instead of member
                # if any of member mismatches + identity, or cluster parent identity exists in transcript-specific list of canonical mod sites:
                # remove this subset
                if any(x in clusterMods_list for x in uss[i]) or any(y in clusterMods_list for y in ussClusterParent):
                    del uss[i]

            if not uss: # if empty after removing uss containing potential mod sites
                notSplit_mods_clusterInfo[cluster].append(isodecoder)
                continue

            temp_uss_dict[isodecoder] = uss # add full uss set to temp dict
            # choose the minimal unique subset
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

    # add "-1" to name so that referencing this isodecoder by its full gene name is possible
    for cluster, data in outputDict.items():
        for uss, isodecoder in data.items():
            isodecoder_list = [x for x in tRNA_dict.keys() if isodecoder[0] in x and not "chr" in isodecoder[0]]
            iso_min = min([x.split("-")[-1] for x in isodecoder_list])
            newIso = isodecoder[0] + "-" + str(iso_min) if not "chr" in isodecoder[0] else isodecoder[0]
            newData = [newIso, isodecoder[1]]
            outputDict[cluster][uss] = newData
    for cluster, data in notSplit_mods_clusterInfo.items():
        for i, isodecoder in enumerate(data):
            isodecoder_list = [x for x in tRNA_dict.keys() if isodecoder in x and not "chr" in isodecoder]
            iso_min = min([x.split("-")[-1] for x in isodecoder_list])
            newIso = isodecoder + "-" + str(iso_min) if not "chr" in isodecoder else isodecoder
            notSplit_mods_clusterInfo[cluster][i] = newIso

    return(outputDict, notSplit_mods_clusterInfo)

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

                isodecoder = "-".join(member.split("-")[:-1]) if not "chr" in member else member
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
    unique_isodecoderMMs, notSplit_mods_clusterInfo = findUniqueSubset(cluster_MemberMismatchPos, unique_isodecoderMMs, tRNA_dict)
    
    # Check that all unique sequences can be deconvoluted
    # count clusters in cluster_dict that are composed of only one sequence and those that are composed of multiple sequences
    singleSeq_clusters_num = 0
    multiSeq_size = 0
    multiSeq_names = set()
    for cluster, members in cluster_dict.items():
        child_iso_set = {tRNA_dict[member]['sequence'].upper() for member in members}
        if len(child_iso_set) == 1:
            singleSeq_clusters_num += 1
        else:
            multiSeq_size += len(child_iso_set)
            multiSeq_names.update({"-".join(member.split("-")[:-1]) for member in members})

    # count deconvoluted sequences above (stored in unique_isodecoderMMs)
    deconv_sequences_num = 0
    deconv_names = set()
    for cluster, data in unique_isodecoderMMs.items():
        unique_deconv = {member[0] for member in data.values()}
        deconv_sequences_num += len(unique_deconv)
        if len(unique_deconv) == len(cluster_MemberMismatchPos[cluster]):
            deconv_sequences_num += 1
            deconv_names.update({"-".join(cluster.split("-")[:-1])})
        deconv_names.update({"-".join(member[0].split("-")[:-1]) for member in data.values()})

    # convert to type defaultdict(set)    
    notSplit_mods_clusterInfo = defaultdict(set, ((k, set(v)) for k, v in notSplit_mods_clusterInfo.items()))
    
    log.info("Total unique sequences: {}".format(singleSeq_clusters_num + multiSeq_size)) # Total unique sequences 
    log.info("Single sequence clusters: {}".format(singleSeq_clusters_num))
    log.info("Deconvoluted sequences: {}".format(deconv_sequences_num))
    log.info("Single transcript resolution for {} sequences".format(singleSeq_clusters_num + deconv_sequences_num))
    log.info("Deconvolution not possible for {} sequences".format((multiSeq_size - deconv_sequences_num)))

 
    # print warnings if not all clusters deconvoluted
    nondeconv_isos = list(multiSeq_names - deconv_names)
    if ((singleSeq_clusters_num + multiSeq_size) != (singleSeq_clusters_num + deconv_sequences_num)) or (notSplit_mods_clusterInfo):
        log.warning("\n*** Deconvolution for some unique transcripts not possible ***\n\
*** due to overlap of distinguishing mismatches with canonical, potentially modified sites ***\n\
This might also be due to errors in input sequence naming in reference.\n\
Ensure all transcripts with identical isodecoder numbers are indeed unique in mature sequence (see *_tRNATranscripts.fa).")
        log.warning("Continuing with analysis - clusters with unsplit transcripts will be renamed accordingly in output.")

    # nondeconv_isos above equivalent to old splitBool variable from earlier code. 
    # unnest the dictionary into list and return for filtering these from DESeq2 analysis and metric reporting
    splitBool = notSplit_mods_clusterInfo

    return(unique_isodecoderMMs, splitBool)

def covCheck_mp(bedTool, unique_isodecoderMMs, splitBool, covDiff, input):
    # get positional coverage per cluster per bam file and check if 3':5' coverage is greater than covDiff
    unsplit = set()
    unsplit_isosOnly = set()
    log.info("Calculating nucleotide coverage for {}".format(input))
    bam = pybedtools.BedTool(input)
    cov = bedTool.coverage(bam, s = True, d = True)
    cov_df = cov.to_dataframe()

    # check coverage diff for each unique sequence in each cluster that is still able to be split
    # i.e. excluded unsplit clusters due to mismatches overlapping canonical mod positions defined in splitIsodecoder()
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

                if float(ratio) < covDiff:
                    unsplit.add(cluster)
                    unsplit.add(unique_isodecoderMMs[cluster][mismatch][0])
                    unsplit_isosOnly.add(unique_isodecoderMMs[cluster][mismatch][0])
                    splitBool[cluster].add(unique_isodecoderMMs[cluster][mismatch][0])


    return(unsplit, unsplit_isosOnly, splitBool)

def unsplitClusters(coverageData, coverageBed, unique_isodecoderMMs, splitBool, threads, covDiff = 0.5):
    # Define unique sequences not able to be split based on significant drop in coverage before distinguishing mismatch
   
    log.info("\n+---------------------------------------------------------------------------------+\
        \n| Determining un-deconvoluted clusters due to insufficient coverage at mismatches |\
        \n+---------------------------------------------------------------------------------+")

    log.info("*** Un-deconvoluted sequences being defined based on coverage difference more than {:.2%}".format(covDiff))
    baminfo, bamlist = getBamList(coverageData)
    
    if len(baminfo) > threads:
	    multi = threads
    else:
	    multi = len(baminfo)

    # initiate custom non-daemonic multiprocessing pool and run with bam names
    log.info("Determining unsplittable sequences...")
    pool = Pool(multi)
    bed = pybedtools.BedTool(coverageBed)
    func = partial(covCheck_mp, bed, unique_isodecoderMMs, splitBool, covDiff)
    unsplit, unsplit_isosOnly, splitBool = zip(*pool.map(func, bamlist))
    pool.close()
    pool.join()

    # create newSplitBool by updating with new cov unsplit transcripts
    newSplitBool = defaultdict(set)
    newSplitBool = {cluster: set().union(data) for proc in splitBool for cluster, data in proc.items()}

    # update unique_isodecoderMMs by removing new unsplit transcripts
    unique_isodecoderMMs_new = copy.deepcopy(unique_isodecoderMMs)
    for cluster, isos in newSplitBool.items():
        if cluster in unique_isodecoderMMs.keys():
            for mismatch, data in unique_isodecoderMMs[cluster].items():
                if data[0] in isos:
                    del unique_isodecoderMMs_new[cluster][mismatch]
    
    unsplit_all = list(set().union(*unsplit))
    unsplit_isosAll = len(list(set().union(*unsplit_isosOnly)))
    log.info("{} unique sequences not split from cluster parent".format(unsplit_isosAll))
    log.info("{} parents and isododecoders not deconvoluted due to insufficient coverage at mismatches".format(len(unsplit_all)))
    log.info("{} total unique sequences not deconvoluted due to mismatches at modified sites OR insufficient coverage".format(len(newSplitBool) + sum([len(data) for data in newSplitBool.values()])))

    return(newSplitBool, unique_isodecoderMMs_new)

def getDeconvSizes(splitBool, tRNA_dict, cluster_dict, unique_isodecoderMMs):
    # get isodecoder sizes in the case of deconvolution
    # this will be a mix of unique transcripts (and their numbers) and unsplit clusters with naming convention: tRNA-AA-Anti-IsoX/Yand their sizes
    
    isodecoder_sizes = defaultdict(int)

    unsplitIsos = [x for k in splitBool.values() for x in k]
    unsplitClusters = [x for x in splitBool.keys()]

    # add single seq cluster sizes
    for cluster, members in cluster_dict.items():
        child_iso_set = {tRNA_dict[member]['sequence'].upper() for member in members}
        if len(child_iso_set) == 1:
            isodecoder_sizes[cluster] = 1

    # generate new names for unsplit clusters (tRNA-AA-Anti-IsoX/Y) - write to lookup for later
    unsplitCluster_lookup = defaultdict()
    for cluster, isos in splitBool.items():
        member_IsoNums = tuple(int(iso.split("-")[-2]) for iso in isos)
        member_IsoNums = sorted(member_IsoNums)
        member_IsoString = "/" + "/".join([str(iso) for iso in member_IsoNums])
        newClusterName = "-".join(cluster.split("-")[:-1]) + member_IsoString
        shortClusterName = "-".join(cluster.split("-")[:-1])
        unsplitCluster_lookup[shortClusterName] = newClusterName

    # add unsplit cluster sizes to isodecoder_sizes with new names
    for cluster, isos in splitBool.items():
        shortClusterName = "-".join(cluster.split("-")[:-1])
        clusterSize = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[cluster]['sequence'].upper()])
        totalIsoSize = 0
        for isodecoder in isos:
            totalIsoSize += len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[isodecoder]['sequence'].upper()])
        isodecoder_sizes[unsplitCluster_lookup[shortClusterName]] = clusterSize + totalIsoSize

    # add split cluster sizes
    for cluster, data in unique_isodecoderMMs.items():
        if not cluster in unsplitClusters:
            isodecoder_sizes[cluster] = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[cluster]['sequence'].upper()])
        for isodecoder in data.values():
            isodecoder_sizes[isodecoder[0]] = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[isodecoder[0]]['sequence'].upper()])

    return(isodecoder_sizes, unsplitCluster_lookup)

def writeDeconvTranscripts(out_dir, experiment_name, tRNA_dict, isodecoder_sizes):
    # write unique transcript and unsplit cluster sequences and size info in the case of deconvolution 

    # save isodecoder info
    with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
        isodecoderInfo.write("Isodecoder\tsize\n")
        for isodecoder, size in isodecoder_sizes.items():
            isodecoder = "-".join(isodecoder.split("-")[:-1]) if not re.search("chr|/", isodecoder) else isodecoder
            isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	# write isodecoder fasta for alignment and context analysis
    with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
        for seq in isodecoder_sizes.keys():
            # new seq to match unsplit cluster parent to tRNA_dict
            if "/" in seq:
                tempSeq = seq.split("/")[0]
                isodecoder_list = [x for x in tRNA_dict.keys() if tempSeq in x and not "chr" in isodecoder[0]]
                iso_min = min([x.split("-")[-1] for x in isodecoder_list])
                newSeq = tempSeq + "-" + str(iso_min)
            else:
                newSeq = seq
            # make shortnames for fa file
            shortname = "-".join(seq.split("-")[:-1]) if not re.search("chr|/", seq) else seq.split("/")[0]
            tempSeqs.write(">" + shortname + "\n" + tRNA_dict[newSeq]['sequence'] + "\n")
    aligntRNA(tempSeqs.name, out_dir)

def writeSplitInfo(out, name, splitBool):
    # output info about unsplit clusters

    with open(out + name + "_unsplitClusterInfo.txt", "w") as unsplitOut:
        unsplitOut.write("Parent\tSize (total unsplit transctipts)\n")
        for cluster, members in splitBool.items():
            size = len(members) + 1 # add one to incude parent in size
            unsplitOut.write(cluster + "\t" + str(size) + "\n")
    log.info("Info on all clusters not deconvoluted written to {}".format(out + "annotation/" + name + "_unsplitClusterInfo.txt"))

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