"""
Created by Robin Aguilar
Beliveau and Noble Labs
Date Created: 11.02.2020
Purpose: Select top 24 probes one, from each chromosome
"""

#load libraries
import time
import io
import sys
import re
import csv
import os
import numpy as np
import pandas as pd
import glob
import itertools
from os import listdir
from os.path import isfile, join
import subprocess
from subprocess import check_output
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import combinations
import argparse
from itertools import product 
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pDups as pDuplex
import matplotlib.pyplot as plt
import seaborn as sns

start_time=time.time()

###################################################################################

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--probe_file', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')

args = userInput.parse_args()
probe_file = args.probe_file

###################################################################################

def read_probe_file(probe_file):

    colnames = ["chrom","probe_start","probe_end", "probe", "Tm", "region", "r_count", "h_count", "pb_score", "alpha_label", "sat_label"]

    probe_df = pd.read_csv(probe_file, names=colnames, header=None,
                        delim_whitespace=True)

    return probe_df

###################################################################################

def parse_selected_probes(probe_df):

    probe_df['pb_score'] = probe_df['pb_score'].astype(float)

    probe_df = probe_df.sort_values(by=['r_count','pb_score'], ascending=False)

    #do groupby

    grouped_df = probe_df.groupby('chrom')

    probe_list = []

    for name, group in grouped_df:
        probe_store_list = group['probe'].tolist()
        probe_list.append(probe_store_list[0])

    saved_probes = probe_df[probe_df['probe'].isin(probe_list)]

    return saved_probes

###################################################################################

def generate_pairs(saved_probes):

    saved_probes_list = saved_probes['probe'].tolist()

    combinations = list(product(saved_probes_list,repeat = 2))

    #split list of tuples into dataframe

    seq_pairs_df = pd.DataFrame(combinations, columns=['p1', 'p2'])

    return seq_pairs_df

###################################################################################

def generate_rc(seq_pairs_df):

    seq2_list=seq_pairs_df['p2'].tolist()
    sequences = [Seq(seq,generic_dna) for seq in seq2_list]
    seqRC=[str(seq.reverse_complement()) for seq in sequences]

    seq_pairs_df['p2_RC']=seqRC

    return seq_pairs_df

###################################################################################

def run_nupack(seq_pairs_df):

    seq1_list = seq_pairs_df['p1'].tolist()
    seq2RC_list = seq_pairs_df['p2_RC'].tolist()

    dup_vals_list=[pDuplex.calculateDup(p1,p2,'74.5','0.390') for p1,p2 in zip(seq1_list,seq2RC_list)]
            
    calc_dup_val=[((float(dup_val))/(float('1.0e-12'))) for dup_val in dup_vals_list]
    
    seq_pairs_df['pDups_seq1_seq2_RC']=calc_dup_val

    return seq_pairs_df

###################################################################################

def plot_heatmap(seq_pairs_df):

    seq_pairs_df = seq_pairs_df.drop('p2_RC',axis = 1)

    seq_pairs_df = seq_pairs_df.pivot("p1","p2","pDups_seq1_seq2_RC")

    sns.heatmap(seq_pairs_df)

    plt.show()

    plt.savefig('probe_pair_heatmap.png')

###################################################################################

def generate_probe_order(saved_probes):

    #split the region column into workable columns
    saved_probes[['chr','start_end']] = saved_probes['region'].astype(str).str.split(':', expand=True)
    saved_probes[['start', 'end']] = saved_probes['start_end'].astype(str).str.split('-', expand=True)

    saved_probes['probe_start'] = saved_probes['probe_start'].astype(int)
    saved_probes['probe_end'] = saved_probes['probe_end'].astype(int)
    saved_probes['start'] = saved_probes['start'].astype(int)
    saved_probes['end'] = saved_probes['end'].astype(int)

    saved_probes['p_start'] = saved_probes['start'] + saved_probes['probe_start']
    saved_probes['p_end'] = saved_probes['start'] + saved_probes['probe_end']

    #rename and adjust columns 
    
    saved_probes = saved_probes.drop(["probe_start","probe_end","chr","start_end","start","end"],axis = 1)

    saved_probes = saved_probes[['chrom', 'p_start', 'p_end', 'probe', 'Tm', 'r_count', 'h_count', 'pb_score', 'alpha_label', 'sat_label']]

    saved_probes.to_csv('tigerfish_alpha_sat_probes.bed',sep='\t', index=False) 
    print(saved_probes)

###################################################################################

def main():

    probe_df = read_probe_file(probe_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    saved_probes = parse_selected_probes(probe_df)

    print("---%s seconds ---"%(time.time()-start_time))

    seq_pairs_df = generate_pairs(saved_probes)

    print("---%s seconds ---"%(time.time()-start_time))

    seq_pairs_df = generate_rc(seq_pairs_df)

    print("---%s seconds ---"%(time.time()-start_time))

    seq_pairs_df = run_nupack(seq_pairs_df)

    print("---%s seconds ---"%(time.time()-start_time))

    plot_heatmap(seq_pairs_df)

    print("---%s seconds ---"%(time.time()-start_time))

    generate_probe_order(saved_probes)

    print("Done")

if __name__== "__main__":
    main()
