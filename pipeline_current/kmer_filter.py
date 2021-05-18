"""
Created on Thu Oct 22 14:50:14 2020
Robin Aguilar
Beliveau and Noble Labs
Unit Testing Code Behavior
"""

#first you should load the libraries
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
import operator

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
requiredNamed.add_argument('-ch', '--chr_name', action='store', required=True,
                               help='The chromosome name')
requiredNamed.add_argument('-f','--fasta', action='store',required=True,
                               help = 'The fasta file for all repeat regions')          
userInput.add_argument('-m', '--merlength', action='store', default=18,
                           type=int,
                           help='The size of k-mers, should be same as jf build k-mer value; default is 18')
userInput.add_argument('-e', '--enrich_score', action='store', default=0.50,required = True, 
                           type=float,
                           help='The max # of times a given k-mer is found in a repeat region/entire HG, default is '
                                '0.50')
userInput.add_argument('-cn', '--copy_num', action='store', default=10,required = True,
                           type=int,
                           help='The minimum allowed number of times a k-mer is identified to be added to the final probe list, default is '
                                '10')
userInput.add_argument('-w', '--span_length', action='store', default=3000,
                           type=int,
                           help='The length of the scanning window for kmer enriched regions; default is '
                                '3000')
userInput.add_argument('-t', '--threshold', action='store', default=10,
                           type=int,
                           help='The minimum number of counts that defines a kmer enriched region; default is '
                                '10')
userInput.add_argument('-c', '--composition_score', action='store', default=0.5,
                           type=float,
                           help='The minimum percentage of kmers that pass threshold in window; default is '
                                '0.5')
requiredNamed.add_argument('-o', '--out_path', action='store', required=True,
                               help='The kmer query file from jellyfish')
args = userInput.parse_args()
probe_file=args.probe_file
jf_file=args.jf_file
chr_name=args.chr_name
fasta_file=args.fasta
MERLENGTH=args.merlength
ENRICH=args.enrich_score
COPY_NUM=args.copy_num
SPAN = args.span_length
THRESHOLD = args.threshold
COMPOSITION = args.composition_score
o_path = args.out_path

def read_probe_file(probe_file):
    #read probes csv
    colnames = ["chr", "p_start", "p_stop", "probe","Tm","regions"]
    probe_data=pd.read_csv(probe_file, delimiter = '\t',names=colnames)
 
    return probe_data

##############################################################################

def read_fasta_dict(fasta_file):   
    #parse out the probe file from blockparse
    seq_dict = {rec.id : rec.seq.upper() for rec in SeqIO.parse(fasta_file, "fasta")}

    return seq_dict

##############################################################################

def local_repeat_count(probe_data,seq_dict,jf_file):
    #now recreate the function you had to get the counts
    #now do a group by

    #read jellyfish file into two lists
    with open(jf_file) as jf:
        rows = (line.split() for line in jf)
        jf_dict = {row[0]:int(row[1]) for row in rows}

    #with dictionary
    all_max_hg_counts = []
    all_max_repeat_counts = []

    all_max_repeat_kmers = []

    grouped_regions = probe_data.groupby(['regions'],sort = False)
    for name,group in grouped_regions:
        if name in seq_dict:
            #make a list of the probes in this group
            probes_list = group["probe"].tolist()
            region_mers = generate_kmers(str(seq_dict[name]),int(MERLENGTH))
            region_mers = Counter(region_mers)
            for probe in probes_list:
                probe_region_count_list = []
                probe_genome_count_list = []
                region_mers_list = []
                probe_mers = generate_kmers(probe,int(MERLENGTH))
                for mer in probe_mers:
                    if mer in region_mers:
                        region_mers_list.append(mer)
                        probe_region_count_list.append(region_mers[mer])
                        region_mer_count_dict = dict(zip(region_mers_list,probe_region_count_list))
                        #return k-mer associated with max value
                        #max_mer_list = [k for k, v in region_mer_count_dict.items() if v == max(region_mer_count_dict.values())]
                        #store that k-mer and add to the all_max list type
                        max_mer_list = [k for k,v in region_mer_count_dict.items()]
                    if mer in jf_dict:
                        probe_genome_count_list.append(jf_dict[mer])
                all_max_repeat_counts.append(max(probe_region_count_list))
                all_max_hg_counts.append(max(probe_genome_count_list))
                all_max_repeat_kmers.append(max_mer_list)

    return all_max_repeat_counts,all_max_hg_counts,all_max_repeat_kmers

##############################################################################

def generate_kmers(sequence,k_size):

    #make a list to store the kmers
    kmers = []
    
    #compute size of kmer window
    n_kmers = len(sequence) - int(k_size) + 1

    #generate kmer within window
    for i in range(n_kmers):
        kmer = sequence[i:i + int(k_size)]
        kmers.append(kmer)

    return kmers

##############################################################################
    
def append_probe_df(probe_data,all_max_repeat_counts,all_max_hg_counts,all_max_repeat_kmers):
    
    probe_data["r_count"] = all_max_repeat_counts
    probe_data['r_count']=probe_data['r_count'].astype(int)

    probe_data["hg_count"] = all_max_hg_counts
    probe_data['hg_count']=probe_data['hg_count'].astype(int)

    probe_data['max_repeat_kmer']=all_max_repeat_kmers

    probe_data['k_score'] = probe_data["r_count"]/probe_data["hg_count"]

    probe_data.to_csv(str(o_path), header= False, index=False, sep="\t")

##############################################################################

def main():
    
    probe_data = read_probe_file(probe_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    seq_dict = read_fasta_dict(fasta_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    all_max_repeat_counts,all_max_genome_counts,all_max_repeat_kmers = local_repeat_count(probe_data,seq_dict,jf_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    append_probe_df(probe_data,all_max_repeat_counts,all_max_genome_counts,all_max_repeat_kmers)

    print("---%s seconds ---"%(time.time()-start_time))

if __name__== "__main__":
    main()
