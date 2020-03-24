"""
Robin Aguilar
Beliveau and Noble Labs
Created: 02/18/20
Purpose: The goal is to decompose probes within RR into their kmers and filter probes with
the same kmer
"""

import pandas as pd
import numpy as np
import time
import sys
from nltk.metrics import *
from itertools import combinations 

from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn import svm, datasets
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics

from Bio import SeqIO
from Bio import pairwise2   
from Bio.pairwise2 import format_alignment 
from Bio.SeqUtils import GC
import itertools
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

"""
userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--probe_file', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')

requiredNamed.add_argument('-o', '--out_file', action='store', required=True,
                               help='The desired name of the output file')

requiredNamed.add_argument('-k', '--k_val', action='store', required=True,
                               help='The k-mer size to filter, default=15')

requiredNamed.add_argument('-n', '--n_k_pass', action='store', required=True,
                               help='Permissible number of shared k-mers to store probe, default=3')
args = userInput.parse_args()
p_file = args.probe_file
o_file=args.out_file
k_size=args.k_val
n_k_pass=args.n_k_pass

"""

p_file="chr19_probes_final_score.txt"
o_file="chr19_probes_out"
k_size=10
n_k_pass=0
    
kmer_set=set()
probe_set=set()
        
kmer_set_non=set()
probe_set_non=set()

##########

all_kmer_set=set()
all_probe_set=set()

all_kmer_set_non=set()
all_probe_set_non=set()

start_time=time.time()

##############################################################################

def build_kmers(sequence,k_size):
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

def read_filter_file(p_file):
    #read in as a dataframe
    colnames=['probe','region', 'window','probe_start','r_start_i',
              'r_end_i', 'p_start_i','p_end_i','hg_num','repeat_num',
              'enrichment']

    probe_file = pd.read_csv(p_file, names=colnames, header=None, sep='\t')

    #drop columns you don't need
    probe_file=probe_file.drop(['window','r_start_i','r_end_i',
                                'p_start_i','p_end_i'], axis=1)
    
    return probe_file

##############################################################################

def generate_ordered_df(probe_file):
    
    #I wonder if you want to use
    #let's sort in descending order based on the repeat number 
    probe_file=probe_file.sort_values(by='repeat_num', ascending=False)
    
    #remove probe duplicates
    probe_file = probe_file.drop_duplicates(subset=['probe'], keep='first')
    
    #let's add a columns for probe end
    probe_file["probe_length"]= probe_file["probe"].str.len() 
    probe_file['probe_end']=probe_file['probe_length']+probe_file['probe_start']
    
    #make the probes into a list
    probes_list=probe_file['probe'].tolist()
    
    idx_list=list(probe_file.index)
    
    #make dictionary of probe idxs
    idx_probe_dict = dict(zip(idx_list, probes_list))
    
    return probe_file,probes_list,idx_probe_dict

##############################################################################

def filter_region_probes(probe_file,probes_list):
    
    #group the regions within the dataframe
    grouped_regions = probe_file.groupby('region')

    #for each of the grouped regions in the dataframe
    for name, group in grouped_regions:
        
        #take each of the probes in that region and add it to a new list
        list_of_probes=group["probe"].tolist()
        
        #run filter probe function
        probe_set,probe_set_non = filter_probes(list_of_probes) 
                    
    return probe_set,probe_set_non

##############################################################################

def filter_probes(list_of_probes):  
    
    #add the first probe and it's kmers into the defined sets
    sequence=list_of_probes[0]
    probe_set.add(sequence)
    first_probe=(build_kmers(sequence,int(k_size)))
    kmer_set=set(first_probe)
        
    #take your range of probes in the list
    for i in range(1,len(list_of_probes)):
        #build the kmers for each one
        sequence=list_of_probes[i]
        probe_mers=build_kmers(sequence,int(k_size))
        #make the kmers into a set
        probe_mers_set=set(probe_mers)

        #perform a set intersection with the head probe
        c = kmer_set.intersection(probe_mers_set)
            
        if len(c) <= int(n_k_pass): 
            kmer_set.update(probe_mers_set)
            probe_set.add(sequence)
        else:
            kmer_set_non.update(probe_mers_set)
            probe_set_non.add(sequence)
    
    return probe_set, probe_set_non

##############################################################################

def parse_probe_sets(probe_set,probe_set_non,idx_probe_dict):
    
    kept_probe_df_idx=[]
    
    rm_probe_df_idx=[]
    
    #identify which indices correspond to what probe
    for key,val in idx_probe_dict.items():
        if val in probe_set:
            kept_probe_df_idx.append(key)
        if val in probe_set_non:
            rm_probe_df_idx.append(key)
            
    return kept_probe_df_idx,rm_probe_df_idx

##############################################################################

def subset_df_rows(probe_file,kept_probe_df_idx,rm_probe_df_idx):
    
    #generate probe dataframes according to kept and discarded
    kept_probe_df=probe_file.loc[kept_probe_df_idx,:]
        
    rm_probe_df=probe_file.loc[rm_probe_df_idx,:]
        
    return kept_probe_df,rm_probe_df

##############################################################################

def filter_between_rr(kept_probe_df,rm_probe_df):
        
    all_probe_list=kept_probe_df['probe'].tolist()
    
    #now take the remaining probes and filter between repeat regions
    sequence=all_probe_list[0]
    all_probe_set.add(sequence)
    first_probe=(build_kmers(sequence,int(k_size)))
    all_kmer_set=set(first_probe)
    
    for i in range(1,len(all_probe_list)):
        #build the kmers for each one
        sequence=all_probe_list[i]
        probe_mers=build_kmers(sequence,int(k_size))
        #make the kmers into a set
        all_probe_mers_set=set(probe_mers)

        #perform a set intersection with the head probe
        c = all_kmer_set.intersection(all_probe_mers_set)
            
        if len(c) <= int(n_k_pass): 
            all_kmer_set.update(all_probe_mers_set)
            all_probe_set.add(sequence)
        else:
            all_kmer_set_non.update(all_probe_mers_set)
            all_probe_set_non.add(sequence)
    
    return all_probe_set, all_probe_set_non

##############################################################################

def parse_probe_sets_global(all_probe_set,all_probe_set_non,idx_probe_dict):
    
    #generate the indices of the remaining seperated groups
    kept_probe_df_idx_global,rm_probe_df_idx_global=parse_probe_sets(all_probe_set,all_probe_set_non,
                                                       idx_probe_dict)
    
    return kept_probe_df_idx_global,rm_probe_df_idx_global

##############################################################################

def subset_row_df_global(probe_file,kept_probe_df_idx_global,rm_probe_df_idx_global,rm_probe_df):
    
    #generate the dataframes withr remaining information
    
    all_kept_probe_df,last_rm_probe_df=subset_df_rows(probe_file,kept_probe_df_idx_global,
                                             rm_probe_df_idx_global)

    #for the rm set, store and append to the last filtered set    
    all_rm_probe_df =pd.concat([rm_probe_df, last_rm_probe_df], axis=1)
    
    return all_kept_probe_df,all_rm_probe_df

##############################################################################

def write_probe_sets(all_kept_probe_df,all_rm_probe_df):
    
    all_kept_probe_df.to_csv(str(o_file)+ ".txt", sep='\t',index=False,header=False)
            
    all_rm_probe_df.to_csv(str(o_file)+ "_n_" + str(n_k_pass) + ".txt", sep='\t',index=False,header=False)

##############################################################################

def main():
    
    probe_file=read_filter_file(p_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    probe_file,probes_list,idx_probe_dict=generate_ordered_df(probe_file)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    probe_set,probe_set_non = filter_region_probes(probe_file,probes_list)
    
    print("---%s seconds ---"%(time.time()-start_time))

    kept_probe_df_idx,rm_probe_df_idx=parse_probe_sets(probe_set,probe_set_non,
                                                       idx_probe_dict)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    kept_probe_df,rm_probe_df=subset_df_rows(probe_file,kept_probe_df_idx,
                                             rm_probe_df_idx)

    print("---%s seconds ---"%(time.time()-start_time))

    filter_between_rr(kept_probe_df,rm_probe_df)
    
    print("---%s seconds ---"%(time.time()-start_time))

    kept_probe_df_idx_global,rm_probe_df_idx_global=parse_probe_sets_global(all_probe_set,
                                                                            all_probe_set_non,
                                                                            idx_probe_dict)
    print("---%s seconds ---"%(time.time()-start_time))

    
    all_kept_probe_df,last_rm_probe_df=subset_row_df_global(probe_file,kept_probe_df_idx_global,
                                                            rm_probe_df_idx_global,
                                                            rm_probe_df)
    
    print("---%s seconds ---"%(time.time()-start_time))

    write_probe_sets(all_kept_probe_df,last_rm_probe_df)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")
    
if __name__== "__main__":
    main()
