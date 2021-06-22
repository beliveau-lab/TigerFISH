"""
Created on Thu Oct 22 14:50:14 2020
Robin Aguilar
Beliveau and Noble Labs
Unit Testing Code Behavior
"""

#first you should load the libraries
import time
import argparse
import subprocess
import numpy as np
import pandas as pd
from itertools import groupby
import collections
from collections import Counter
from Bio import SeqIO
from Bio.Seq import reverse_complement as rev_comp

##############################################################################

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
requiredNamed.add_argument('-o', '--out_path', action='store', required=True,
                               help='The kmer query file from jellyfish')
requiredNamed.add_argument('-c1', '--c1_value', action='store', required=True,type=int,
                               help='The kmer query file from jellyfish')
requiredNamed.add_argument('-c2', '--c2_value', action='store', required=True,type=int,
                               help='The kmer query file from jellyfish')

args = userInput.parse_args()
probe_file=args.probe_file
jf_file=args.jf_file
chr_name=args.chr_name
fasta_file=args.fasta
MERLENGTH=args.merlength
o_path = args.out_path
c1_val = args.c1_value
c2_val = args.c2_value

##############################################################################

def read_probe_file(probe_file):

    #read probes csv
    colnames = ["chrom","p_start","p_stop","probe","Tm", "regions"]

    probe_data=pd.read_csv(probe_file, delimiter = '\t',names=colnames)
 
    probe_data = probe_data.drop_duplicates(subset=['probe'], keep='first')

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

    jf_dict = {}

    #read jellyfish file into two lists
    with open(jf_file) as jf:
        rows = (line.split() for line in jf)
        for row in rows:
            jf_dict[row[0]] = int(row[1])

    #with dictionary
    all_hg_counts = {}
    all_repeat_counts = {}

    grouped_regions = probe_data.groupby(['regions'],sort = False)
    for name,group in grouped_regions:
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
                #generate RC here
                rev_mer = rev_comp(mer)
                #then check if mer or mer RC in region
                if mer in region_mers:
                    probe_region_count_list.append(region_mers[mer])
                if rev_mer in region_mers:
                    probe_region_count_list.append(region_mers[rev_mer])
                if mer in jf_dict:
                    probe_genome_count_list.append(jf_dict[mer])
                if rev_mer in jf_dict:
                    probe_genome_count_list.append(jf_dict[rev_mer])
            all_repeat_counts[probe] = sum(probe_region_count_list)
            all_hg_counts[probe] = sum(probe_genome_count_list)

    return all_repeat_counts,all_hg_counts

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
    
def append_probe_df(probe_data,all_repeat_counts,all_hg_counts):

    repeat_count_list = []
    hg_count_list = []

    probe_list = probe_data['probe'].tolist()

    for probe in probe_list:
        repeat_count_list.append(all_repeat_counts[probe])
        hg_count_list.append(all_hg_counts[probe])

    #zip all 3 into df and then merge...
    probe_data["r_count_total"] = repeat_count_list
    probe_data['r_count_total']=probe_data['r_count_total'].astype(int)

    probe_data["hg_count_total"] = hg_count_list
    probe_data['hg_count_total']=probe_data['hg_count_total'].astype(int)

    probe_data['k_score_total'] = probe_data["r_count_total"]/probe_data["hg_count_total"]

    probe_data['k_score_total']=probe_data['k_score_total'].astype(float)

    return probe_data

##############################################################################

def compute_normalized_binding(probe_data):

    norm_val_list = []
    all_probe_norm_dict = {}
    #est weights, to be params later maybe?
    c1 = c1_val
    c2 = c2_val

    grouped_regions = probe_data.groupby(['regions'],sort = False)
    for name,group in grouped_regions:
        val_list = []
        #make a list for the probes
        probe_list = group['probe'].tolist()

        #make a list of the probes in this group
        repeat_counts = group["r_count_total"].tolist()
        k_binding = group["k_score_total"].tolist()

        max_r_sum = max(repeat_counts)
        max_k_binding = max(k_binding)

        #do math to produce normalized value
        for rc,kb in zip(repeat_counts,k_binding):
            norm_binding_val = (((rc/max_r_sum)*c1) + ((kb/max_k_binding)*c2))
            val_list.append(norm_binding_val)

        probe_norm_val_dict = dict(zip(probe_list,val_list))
        
        all_probe_norm_dict.update(probe_norm_val_dict)

    #get a total probe list
    probe_list = probe_data['probe'].tolist()
    for probe in probe_list:
        norm_val_list.append(all_probe_norm_dict[probe])

    probe_data['norm_vals'] = norm_val_list

    probe_data = probe_data.sort_values(by=['norm_vals'], ascending=False)

    return probe_data

##############################################################################

def write_file(probe_data):

    probe_data.to_csv(str(o_path), header= False, index=False, sep="\t")

##############################################################################

def main():
    
    probe_data = read_probe_file(probe_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    seq_dict = read_fasta_dict(fasta_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    all_repeat_counts,all_genome_counts = local_repeat_count(probe_data,seq_dict,jf_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    probe_data = append_probe_df(probe_data,all_repeat_counts,all_genome_counts)

    print("---%s seconds ---"%(time.time()-start_time))

    probe_data = compute_normalized_binding(probe_data)

    print("---%s seconds ---"%(time.time()-start_time))

    write_file(probe_data)

    print("---%s seconds ---"%(time.time()-start_time))

if __name__== "__main__":
    main()
