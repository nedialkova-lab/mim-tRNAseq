# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:01:53 2021

@author: Xavier Hernandez-Alias
"""
import pandas as pd
import numpy as np
from os import listdir
from shutil import rmtree
from os.path import isdir,isfile,join
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import multiprocessing
from functools import partial

# Define function for parallelization
def analyze_1sample(s,dirpath,thres):
    outdf = pd.DataFrame(columns=["sample","ref","var1","var2","pval","odds_ratio","values"])
    n=0
    ref_files = [f for f in listdir(join(dirpath, s)) if isfile(join(dirpath, s, f))]
    for f in ref_files:
        # Load table
        ref = f[:-7].replace(".","/")
        ref = "-".join(ref.split("-")[:-1]) if not "chr" in ref and not "/" in ref else ref
        readsdf = pd.read_csv(join(dirpath, s, f),sep="\t",compression="gzip",index_col="READ",dtype="category")
        pairs = []
        for v1 in readsdf.columns:
            for v2 in readsdf.columns:
                if v1!=v2 and set([v1,v2]) not in pairs:
                    pairs.append(set([v1,v2]))
                    reads1 = readsdf[v1].dropna()
                    reads2 = readsdf[v2].dropna()
                    # Do test only if at least % positions are modified
                    if reads1.shape[0]>0 and reads2.shape[0]>0:
                        if sum(reads1!="0")/reads1.shape[0]>thres and sum(reads2!="0")/reads2.shape[0]>thres:
                            # Keep only reads that contain both v1 and v2
                            tempdf = readsdf[[v1,v2]].dropna(how="any")
                            # Build contingency table, var1 in rows and var 2 in columns
                            counts = (tempdf!="0").value_counts()
                            if counts.shape[0]==4:
                                cont_tab = np.array([[counts.loc[True,True], counts.loc[True,False]],
                                                     [counts.loc[False,True], counts.loc[False,False]]])
                                #p = chi2_contingency(cont_tab)[1]
                                oddsr, p = fisher_exact(cont_tab)
                                outdf.loc[n] = [s,ref,v1,v2,p,oddsr,counts]
                                n += 1
    return outdf

def crosstalks_wrapper(dirpath, thres, threads):
    # Multiprocessed function
    pool = multiprocessing.Pool(threads)
    frozen_fun = partial(analyze_1sample, dirpath=dirpath, thres=thres)
    samples = [f for f in listdir(dirpath) if isdir(join(dirpath, f))]
    dfs = pool.map(func=frozen_fun, iterable=samples, chunksize=1)
    pool.close()
    pool.join()
    # Concatenate tables
    outdf = pd.concat(dfs, ignore_index=True)
    # Correct multiple comparisons
    outdf["pval_corrected"] = multipletests(outdf["pval"].values,method="fdr_bh")[1]
    # Include canonical positions
    posinfo = pd.read_csv("/".join(dirpath.split("/")[:-1])+"/mods/mismatchTable.csv", sep="\t",
                          index_col=["isodecoder","pos"],usecols=["isodecoder","pos","canon_pos"],
                          dtype="category")
    posinfo = posinfo[~posinfo.index.duplicated(keep='first')]
    outdf["canon_var1"] = ["Charged" if s[1]=="Charged" else posinfo.loc[(s[0],s[1]),"canon_pos"] if s[0] in posinfo.index.get_level_values(0) else "NA" for s in outdf[["ref","var1"]].to_numpy()]
    outdf["canon_var2"] = ["Charged" if s[1]=="Charged" else posinfo.loc[(s[0],s[1]),"canon_pos"] if s[0] in posinfo.index.get_level_values(0) else "NA" for s in outdf[["ref","var2"]].to_numpy()]
    # Save table
    outdf.to_csv(dirpath+"/crosstalks.tsv",sep="\t",index=False)
