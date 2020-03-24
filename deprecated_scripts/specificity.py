"""
Robin Aguilar
Beliveau and Noble Labs, 10/14/2019
Dependencies: Jellyfish, Bedtools
Input: Probe file
Output: A .bed file of probes that pass probe binding score in the form of a DF 
"""

#import all functions/modules needed
import time
import io
import sys
import re
import csv
import argparse
from collections import defaultdict
from itertools import islice
from operator import itemgetter, attrgetter
import subprocess
import numpy as np
import pandas as pd
import glob
import os
import refactoredBlockparse as bp
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
from scipy import signal
import time
from collections import OrderedDict
import collections
from operator import itemgetter
from itertools import groupby
import collections
from collections import Counter
import itertools

#declare a timer
start_time=time.time()

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-p', '--probe_file', action='store', required=True,
                               help='The FASTA file to find probes in')
requiredNamed.add_argument('-j', '--jf_file', action='store', required=True,
                               help='The kmer query file from jellyfish')
requiredNamed.add_argument('-k', '--idx_file', action='store', required=True,
                               help='The kmer query file from jellyfish')
requiredNamed.add_argument('-c', '--chr_file', action='store', required=True,
                               help='The chromosome name')
args = userInput.parse_args()
test_probes=args.probe_file
test_jf=args.jf_file
k_index_f=args.idx_file
chr_name=args.chr_file

MERLENGTH=18

def probe_df(probe_f):
    """
    This function will organize the format of the dataframe containing all probe information.
    """
    probe_data=pd.read_csv(probe_f, delimiter = '\t',names=["index","p_start","p_end","probe_seq","Tm","region"])
    probe_data['region'] = probe_data['region'].astype(str)

    probe_data[['chr','r_start','r_end']] = probe_data.region.str.replace(":","-").str.split("-",expand=True)

    probe_data['r_start']=probe_data['r_start'].astype(int)
    probe_data['r_end']=probe_data['r_end'].astype(int)

    #here you put all the probe sequences into a list
    probe_seqs = probe_data['probe_seq'].tolist()
    
    return probe_data,probe_seqs

indexList=[]
merWindows_lengths=[]

def probe_index(probe_seq_l):
    """
    This function will make an index of locations based on each probe length
    """
    for i in range(0, len(probe_seq_l), 1):
        probeLength = len(probe_seq_l[i])
        merWindows = probeLength - MERLENGTH + 1
        merWindows_lengths.append(merWindows)
    return merWindows_lengths

def append_probe_windows(probe_data_df,merlengths,k_index):
    """
    This function will add the lengths of the probe windows to the probe data DF
    """

    probe_data_df['probe_dist']=merlengths

    with open(k_index,"r") as index_file:
        kmer_indices=[int(line) for line in index_file]

    #make k-mer indices dictionary that keeps track of bases as indices
    kmer_ind_dict={k: v for v, k in enumerate(kmer_indices)}

    #make repeat start and end locations into a list
    repeat_start_indices = probe_data_df['r_start'].tolist()
    probe_data_df['r_end_mer']=probe_data_df['r_end']-18
    #repeat_end_indices = probe_data_df['r_end'].tolist()
    repeat_end_indices = probe_data_df['r_end_mer'].tolist()


    #make a new column that adds that probe start relative to the start of the repeat
    probe_data_df['probe_start_in_seq']=probe_data_df['r_start'].astype(int)+probe_data_df['p_start']
    probe_start_seq=probe_data_df['probe_start_in_seq'].tolist()
    
    #query the k-mer indices dictionary to find the key (sequence location of start) to return the index (val)
    starts=[kmer_ind_dict[item] for item in repeat_start_indices if item in kmer_ind_dict]
    ends=[kmer_ind_dict[item] for item in repeat_end_indices if item in kmer_ind_dict]
    probe_starts=[kmer_ind_dict[item] for item in probe_start_seq if item in kmer_ind_dict]

    #these indices correspond to the line in the file that corresponds to the correct probe locations
    probe_data_df['test_start_idx']=starts
    probe_data_df['test_end_idx']=ends
    probe_data_df['test_p_start']=probe_starts
    probe_data_df['test_p_end']=probe_data_df['test_p_start'].astype(int)+probe_data_df['probe_dist']

    #remove redundant columns
    probe_data_df=probe_data_df.drop(columns=['index', 'p_start','p_end','Tm','chr','r_start','r_end','r_end_mer'])

    return probe_data_df

def get_max_hg_repeat(jf_file,probe_data_df):
    """
    This function will iterate through the probe_data list items and get the max count for probes in the hg and in the repeat
    """
    max_list_region=[]
    max_list_all=[]
    k_mer=[]
    float_count=[]

    #read jellyfish file into two lists
    with open(jf_file,"r") as jf:
        for line in jf:
            k_mer.append(line.split()[0])
            float_count.append(float(line.split()[1])) 

    #convert the lists to numpy arrays
    k_mer = np.array(k_mer)
    float_count = np.array(float_count)

    #convert DF to numpy array
    probe_data_npy = probe_data_df.values

    #iterate over numpy array
    for row in probe_data_npy:
        probes=row[0]
        r_start=row[4]
        r_end=row[5]
        p_start=row[6]
        p_end=row[7]

        #obtain max value for each probe in all of hg38
        probe_kmers_all=(float_count[int(p_start)-1:int(p_end)-1])
        max_list_all.append(max(probe_kmers_all))

        #obtain max value a all 18mers exists in each repeat region/probe
        repeat_kmers=(k_mer[int(r_start):int(r_end)])
        probe_kmers_seqs=(k_mer[int(p_start)-1:int(p_end)-1])
        region_kmers=Counter(repeat_kmers)

        #check if probe k-mer is in repeat region dictionary of k-mers, append val if true
        region_vals=[region_kmers[item] for item in probe_kmers_seqs if item in region_kmers]

        largest=max(region_vals)
        max_list_region.append(largest)

    probe_data_df['max_mer_list']=max_list_all
    probe_data_df['max_repeat']=max_list_region
    probe_data_df["k-mer enrich"]=probe_data_df['max_repeat'].astype(int)/(probe_data_df['max_mer_list'].astype(int))

    probe_data_df=probe_data_df[(probe_data_df['k-mer enrich'] >= 0.50) & (probe_data_df['max_repeat'] >= 10)]
    probe_data_df.to_csv("initial_specificity_out/" + str(chr_name)+'_probes_final_score.txt', header=False, index=False, sep="\t")

    return probe_data_df

def main():

    probe_data,probe_seqs=probe_df(test_probes)
    
    print("---%s seconds ---"%(time.time()-start_time))

    merWindows_lengths=probe_index(probe_seqs)
    
    print("---%s seconds ---"%(time.time()-start_time))

    probe_data_df=append_probe_windows(probe_data,merWindows_lengths,k_index_f)
    
    print("---%s seconds ---"%(time.time()-start_time))

    get_max_hg_repeat(test_jf,probe_data_df)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
