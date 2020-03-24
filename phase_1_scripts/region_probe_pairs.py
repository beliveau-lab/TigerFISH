"""
Robin Aguilar
Beliveau and Noble Labs
Compute probe pairs
Input: Probe filtered file
Output: file of all unique pairwise comparisons of probes
"""
from __future__ import print_function
import pandas as pd
from pandas import DataFrame
from Bio import pairwise2   #for pairwise alignments 
from Bio.pairwise2 import format_alignment   #for printing alignments out neatly 
import itertools
from Bio import Align
import csv
import argparse
from collections import defaultdict
from itertools import islice
from operator import itemgetter, attrgetter
import subprocess
import numpy as np
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
import time
from itertools import combinations
from nltk.metrics import *

start_time=time.time()

#include user arguments
userInput = argparse.ArgumentParser(description=\
        '%Requires a filtered probe file as input. This script will take unique probes and conduct pairwise local alignment. Probes with edit distance <=5 will be considered probes that are similar in sequence composition to one another and will be further filtered.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-p', '--probe_file', action='store', required=True,
                               help='The dataframe containing filtered probes')
requiredNamed.add_argument('-o', '--out_name', action='store', required=True,
                               help='The name of the output file')
args = userInput.parse_args()
test_probes=args.probe_file
out_f=args.out_name

def get_probe_list(probe_filter):
    """
    This function reads in the filtered probe file and appends all probes to region as a dictionary
    """
    #read in the probe file
    filtered=pd.read_csv(probe_filter, delimiter = '\t',names=["probe_seq",'region',"mer_length","p_start_i","r_start_i","r_end_i","p_start_idx","p_end_idx","hg_count","r_count","p_score"])

    #make a dictionary that contains all probes to one region. Region is the key, all probes are values in list
    d=defaultdict(list)
    probe_list = filtered['probe_seq'].tolist()
    region_list=filtered['region'].tolist()
    score_list=filtered['p_score'].tolist()

    #append probes to region dictionary
    for region,probe in zip(region_list,probe_list):
        d[region].append(probe)

    #keep track of probes and their scores
    probe_score_list = dict(zip(probe_list, score_list))

    return probe_list,region_list,d,probe_score_list

def generate_unique_pairs(probe_l,probe_dict):
    """
    This will generate all pairwise possibilities for each probe within each region
    """
    probes_list=[list(combinations(probe_dict[x],2)) for x in probe_dict.keys()]

    return probes_list

def align_probes(probes_list,probe_dict,score_dict):
    """
    This will map all the sequences to their regions so they can be appended to a dataframe
    """
    #append all probes to a particular region
    listed_dict = {k: oldk for oldk, oldv in probe_dict.items() for k in oldv}

    #for all pairwise combinations, add them to the DF as two columns, seq1 and seq2
    align_df=pd.DataFrame([t for lst in probes_list for t in lst],columns=["seq1","seq2"])

    #map the region to the probe and assign them as seq1_region and seq2_region
    align_df['seq1_region'] = align_df['seq1'].map(listed_dict)
    align_df['seq2_region'] = align_df['seq2'].map(listed_dict)
        
    #append the corresponding scores to each probe
    align_df['seq1_score']=align_df['seq1'].map(score_dict)
    align_df['seq2_score']=align_df['seq2'].map(score_dict)

    return align_df

def get_rev_comp(alignment_df):

    """
    This will take the RC of probe 2 and append to DF, this will also write out DF
    """
    comps=alignment_df['seq2'].tolist()
    #use biopython to obtain the RC
    sequences = [Seq(a_seq,generic_dna) for a_seq in comps]
    seqRC=[seq.reverse_complement() for seq in sequences]

    #append as DF column
    alignment_df['RC']=seqRC
    alignment_df=alignment_df.drop_duplicates()

    #write out the file
    alignment_df.to_csv(str(out_f)+'_pairs_scores.txt', header=True, index=False, sep="\t")

def main():

    probe_list,region_list,d,probe_score_list=get_probe_list(test_probes)

    print("---%s seconds ---"%(time.time()-start_time))

    probe_pairs=generate_unique_pairs(probe_list,d)

    print("---%s seconds ---"%(time.time()-start_time))

    align_df=align_probes(probe_pairs,d,probe_score_list)

    print("---%s seconds ---"%(time.time()-start_time))

    get_rev_comp(align_df)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
