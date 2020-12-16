"""
Robin Aguilar
Beliveau and Noble Labs
Date Created: 02/25/2020
Purpose: To refactor the old specificity script to it runs at a much faster rate
"""

#import all the functions you need

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


#come up with a small test file that you can run through the whole pipeline and test speed

start_time=time.time()

# Define the command line usage.
USAGE = """Usage: specificity.py -p probe_file -j jf_file -k idx_file -c chr_name -f fasta_file [options]

  Options:

    -merlength     The size of k-mers used used to compute jellyfish query file (contains 18-mers, count in genome), default = 18

    -enrich_score  max k-mer found in repeat region/max k-mer found in genome for any given probe sequence, default = 0.50

    -copy_num      min number of times a k-mer must be found in a repeat region to be considered in the final list, default = 10

"""

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
requiredNamed.add_argument('-c', '--chr_name', action='store', required=True,
                               help='The chromosome name')
requiredNamed.add_argument('-f','--fasta', action='store',required=True,
                               help = 'The fasta file for all repeat regions')          
userInput.add_argument('-m', '--merlength', action='store', default=18,
                           type=int,
                           help='The size of k-mers, should be same as jf build k-mer value; default is 18')
userInput.add_argument('-e', '--enrich_score', action='store', default=0.50,
                           type=float,
                           help='The max # of times a given k-mer is found in a repeat region/entire HG, default is '
                                '0.50')
userInput.add_argument('-cn', '--copy_num', action='store', default=10,
                           type=int,
                           help='The minimum allowed number of times a k-mer is identified to be added to the final probe list, default is '
                                '10')
args = userInput.parse_args()
probe_file=args.probe_file
jf_file=args.jf_file
k_index_file=args.idx_file
chr_name=args.chr_name
fasta_file=args.fasta
MERLENGTH=args.merlength
ENRICH=args.enrich_score
COPY_NUM=args.copy_num


def probe_df(probe_f):
    """
    This function will organize the format of the dataframe containing all probe information.
    """
    #open probe file with probes generated from each repeat region
    probe_data=pd.read_csv(probe_f, delimiter = '\t',names=["index","p_start","p_end","probe_seq","Tm","region"])
    probe_data['region'] = probe_data['region'].astype(str)

    #split each repeat region into chr, region start, and region stop columns
    probe_data[['chr','r_start','r_end']] = probe_data.region.str.replace(":","-").str.split("-",expand=True)

    #cast the coordinates for region starts and stops as integers
    probe_data['r_start']=probe_data['r_start'].astype(int)
    probe_data['r_end']=probe_data['r_end'].astype(int)
    
    #here you put all the probe sequences into a list
    probe_seqs = probe_data['probe_seq'].tolist()
        
    return probe_data,probe_seqs

##############################################################################

def probe_index(probe_seq_l):
    """
    This function will make an index of locations based on each probe length
    """
    merWindows_lengths=[]

    #for all probes in the length of the list of repeat probes
    for i in range(0, len(probe_seq_l), 1):
        #this is a calculation for the number of 18mers in each probe essentially
        probeLength = len(probe_seq_l[i])
        merWindows = probeLength - MERLENGTH + 1
        merWindows_lengths.append(merWindows)
        
    #return this list to append it as a column in the overall probe dataframe
    return merWindows_lengths

##############################################################################

def append_probe_windows(probe_data_df,merWindows_lengths,k_index_file):
    """
    This function will add the lengths of the probe windows to the probe data DF
    """

    #add the total number of 18-mers in each probe as a column
    probe_data_df['probe_dist']=merWindows_lengths

    #read the k_mer index file and put idx into list
    with open(k_index_file,"r") as index_file:
        kmer_indices=[int(line) for line in index_file]

    #make k-mer indices dictionary that keeps track of bases as indices
    kmer_ind_dict={k: v for v, k in enumerate(kmer_indices)}

    #make repeat start and end locations into a list
    repeat_start_indices = probe_data_df['r_start'].tolist()
    
    #r_end is actually the final base in that repeat
    #to get the last k-mer in the repeat region, you must subtract the length of the k-mer
    probe_data_df['r_end_mer']=probe_data_df['r_end']-int(MERLENGTH)
    repeat_end_indices = probe_data_df['r_end_mer'].tolist()
    
    #make a new column that adds that probe start relative to the start of the repeat
    probe_data_df['probe_start_in_seq']=probe_data_df['r_start'].astype(int)+probe_data_df['p_start']
    probe_data_df['probe_end_in_seq']=probe_data_df['probe_start_in_seq'].astype(int)+probe_data_df['probe_seq'].str.len()-1

    #use this to query and identify the true start index of each probe in the jf file
    probe_start_seq=probe_data_df['probe_start_in_seq'].tolist()
        
    #query the k-mer indices dictionary to find the key (sequence location of start) to return the index (val)
    starts=[kmer_ind_dict[item] for item in repeat_start_indices if item in kmer_ind_dict]
    ends=[kmer_ind_dict[item] for item in repeat_end_indices if item in kmer_ind_dict]
    probe_start=[kmer_ind_dict[item] for item in probe_start_seq if item in kmer_ind_dict]

    #these indices correspond to the line in the file that corresponds to the correct probe locations
    probe_data_df['region_start_idx']=starts
    probe_data_df['region_end_idx']=ends
    probe_data_df['probe_start_idx']=probe_start
    probe_data_df['probe_end_idx']=probe_data_df['probe_start_idx'].astype(int)+probe_data_df['probe_dist']
    probe_end=probe_data_df['probe_end_idx'].tolist()
        
    #remove redundant columns
    probe_data_df=probe_data_df.drop(columns=['index', 'p_start','p_end','chr','r_start','r_end','r_end_mer'])
    
    return probe_data_df,starts,ends,probe_start,probe_end

##############################################################################

def get_max_hg_repeat(jf_file,probe_data_df,probe_start,probe_end):
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
            #add the kmers and counts into seperate lists
            k_mer.append(line.split()[0])
            float_count.append(float(line.split()[1])) 
        
    #parse the counts and obtain the max count for each probe in the whole genome
    probe_kmers_all=[max(float_count[int(p_start)-1:int(p_end)-1]) for p_start,p_end in zip(probe_start,probe_end)]
    
    return k_mer,float_count,probe_kmers_all, probe_data_df

##############################################################################

def build_kmers(sequence,k_size):
    """
    This function will split and string sequence into k_mers of a given size you specify
    """
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

def get_max_repeat(fasta_file,probe_data_df):
    """
    This function will return a list of the occurence of each probes k-mer within a given repeat region.
    From this list, the max (highest count kmer) is computed for each probe and that value is added to the dataframe.
    
    """
    #take the fasta file for the repeat regions designed against a scaffold
    
    #add distinct repeat region k-mers into a larger list
    region_kmers=[]
    
    #use biopython SeqIO to parse out the fasta sequences
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        #make the sequences uppercase so they match the case of the probes in the dataframe
        sequence = str(fasta.seq.upper())
        #break each region into 18mers
        region_mers=(build_kmers(sequence,int(MERLENGTH)))
        #apply a counter to creat a dict that counts the number of times each 18mer occurs in the repeat
        region_kmers.append(Counter(region_mers))
        
    #iterate over the grouped regions
    grouped_regions = probe_data_df.groupby(['region_start_idx','region_end_idx'])
    
    #add all probes withing each repeat region to a list that can be iterated over with the counter
    all_probe_seqs=[]
    for name,group in grouped_regions:
        probe_seqs=group["probe_seq"].tolist()
        all_probe_seqs.append(probe_seqs)

    #this list will store the max k-mer count for all 18mers in a given probe sequence
    all_probes_max=[]
    
    #iterate over the large probe list and the region counter
    for probes,regions in zip(all_probe_seqs,region_kmers): 
        #within each region, make a list of probes
        single_p_list=[]
        #for each of these probes for this region specific list, generate all 18-mers
        for p in probes:
            probe_mers=(build_kmers(p,int(MERLENGTH)))
            #add the 18-mers to the list
            single_p_list.append(probe_mers)
        
        #for each probe and it's 18-mers for each probe
        for l in single_p_list:
            #append the counts given intersecting information for each region Counter
            count_list=[]
            for mer in l: 
                #for each probe, count the number of times each kmer is identified within that repeat region
                count_list.append(regions[mer])
            #append the max of that list for each probe's 18-mer count
            max_val=(max(count_list))
            #add this to the overall list which will be added to the final dataframe
            all_probes_max.append(max_val)
            
    return probe_data_df,all_probes_max

##############################################################################

def compute_kmer_enrichment(probe_data_df,probe_kmers_all,all_probes_max):
            
    #append the two lists for max probe k-mer count in whole genome and repeat region
    probe_data_df['max_hg_count']=probe_kmers_all
    probe_data_df['max_repeat_count']=all_probes_max   
    
    #compute the k-mer enrichment score
    probe_data_df["k-mer enrich"]=probe_data_df['max_repeat_count'].astype(int)/(probe_data_df['max_hg_count'].astype(int))
    
    #filter based on these parameters, should be able to tune these as arguments
    probe_data_df=probe_data_df[(probe_data_df['k-mer enrich'] >= float(ENRICH)) & (probe_data_df['max_repeat_count'] >= int(COPY_NUM))]
        
    #drop columns that are no longer needed at this point
    probe_data_df=probe_data_df.drop(columns=['probe_dist','region_start_idx','region_end_idx',
                                              'probe_start_idx','probe_end_idx'])

    #reorder the columns so they are printed more neatly
    probe_data_df = probe_data_df[['region', 'probe_seq', 'Tm' , 'probe_start_in_seq','probe_end_in_seq',
                                  'max_hg_count','max_repeat_count','k-mer enrich']]
    
    #generate an output dataframe, will need to modify for snakemake
    probe_data_df.to_csv('initial_specificity_out/' +str(chr_name)+'_probes_final_score.txt', header=False, index=False, sep="\t")

##############################################################################
# MAIN
##############################################################################

def main():

    probe_data,probe_seqs=probe_df(probe_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    merWindows_lengths=probe_index(probe_seqs)
    
    print("---%s seconds ---"%(time.time()-start_time))

    probe_data_df,starts,ends,probe_start,probe_end=append_probe_windows(probe_data,merWindows_lengths,k_index_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    k_mer,float_count,probe_kmers_all,probe_data_df=get_max_hg_repeat(jf_file,probe_data_df,probe_start,probe_end)

    print("---%s seconds ---"%(time.time()-start_time))
    
    probe_data_df,all_probes_max=get_max_repeat(fasta_file,probe_data_df)

    print("---%s seconds ---"%(time.time()-start_time))
    
    compute_kmer_enrichment(probe_data_df,probe_kmers_all,all_probes_max)

    print("Done")

if __name__== "__main__":
    main()


