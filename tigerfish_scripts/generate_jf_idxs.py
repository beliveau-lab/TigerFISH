"""
Robin Aguilar
Date Created: 10.18.2020
Beliveau and Noble Labs
Purpose: To generate jellyfish index files and count files to avoid taxing
memory on cluster.
"""

#import all functions/modules needed
import time
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
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
from collections import OrderedDict
import collections
from collections import Counter
import itertools
import re

#start the time
start_time=time.time()

#write arguments so users can specify the commands they would like for each variable
userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing or optionally can be made into a  '
        '.bed file can be outputted instead if \'-b\' is flagged. Tm values '
        'are corrected for [Na+] and [formamide].')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--fasta_file', action='store', required=True,
                               help='The FASTA file to find probes in')
requiredNamed.add_argument('-j', '--jf_indexfile', action='store', required=True,
                               help='The jf file of a given genome')
requiredNamed.add_argument('-chr', '--chr_name', action='store', required=True,
                               help='Define the scaffold being queried')
requiredNamed.add_argument('-schr', '--scaffold_fasta', action='store',required=True,
                           help='Used to generate fasta file of kmer rich regions; default is ')
requiredNamed.add_argument('-st','--start',action='store',required=True,
                           help='The start sequence of the fasta region if not starting at beginning of scaffold; default is ' '0')

args = userInput.parse_args()
fasta_file = args.fasta_file
scaffold_fasta=args.scaffold_fasta
jf_indexfile = args.jf_indexfile
chrom= args.chr_name
START=args.start
MERLENGTH = 18
#the names of output files to be generated
chr_name_jf_out="results/jf_index_files/" + chrom + "_jf_temp.txt"
out_index="results/jf_index_files/" + chrom + "_index.txt"

##############################################################################

def jf_query(jf_indexfile,fasta_file):
    
    """
    Runs jellyfish to generate a query file from the jellyfish index provided
    """
    query_file=subprocess.call(['jellyfish', 'query', jf_indexfile, '-s',
                     fasta_file, '-o', chr_name_jf_out], stderr=None, shell=False)
    
    return chr_name_jf_out

##############################################################################
    
def map_coords(fa_file):
    
    """
    This function will take a the fasta file you pass it and identify where all N and non-N bases are located
    """
    
    bases_dist_start=[]
    n_bases_start=[]
    
    fa_seq = list(SeqIO.parse(open(fa_file),'fasta'))
    for fasta in fa_seq:
        sequence=str(fasta.seq).lower()
        #first you need to find where all of the regular bases are located
        for match in re.finditer('[atcg]',sequence):
            bases_dist_start.append(int(match.start())+int(START))
        #now you need to find where all of the N bases are located 
        for match in re.finditer('[n]', sequence):
            n_bases_start.append(int(match.start()))
    
    return bases_dist_start,n_bases_start

##############################################################################

def group_ranges(bases_dist_start,n_bases_start):
    """
    This function will return lists continous ranges of N and non-N bases.
    """
    ranges = []
    n_ranges=[]
    
    #this collapses the set of the regular bases to provide you with ranges in a list
    for k, g in groupby(enumerate(bases_dist_start), lambda x: x[1]-x[0]):
        group = list(map(itemgetter(1), g))
        ranges.append(str(group[0]) + "\t" + str(group[-1]+1))

    if n_bases_start:
        #this collapses the set of N bases to provide you with ranges in a list 
        for k, g in groupby(enumerate(n_bases_start), lambda x: x[1]-x[0]):
            group=list(map(itemgetter(1),g))
            n_ranges.append(str(group[0]) + "\t" + str(group[-1]+1))
        
    return ranges,n_ranges

##############################################################################

def compute_ranges(base_ranges,label):
    
    """
    This function will generate a dataframe with the start and stop of the ranges
    """
    
    #make a dataframe of bases
    ranges_v=pd.DataFrame(list(base_ranges), columns=['ranges'])
    #split the start and end to be seperate columns
    ranges_v[['start','end']] = ranges_v.ranges.str.split(expand=True)
    num_start=ranges_v['start'].astype(int)
    num_end=ranges_v['end'].astype(int)
    num_ranges = pd.concat([num_start, num_end], axis=1)
    #label whether the range is of a normal base "R" or an N, "N"
    num_ranges['type']=str(label)
    num_ranges['size']=num_ranges['end'].astype(int)-num_ranges['start'].astype(int)
    
    return num_ranges

##############################################################################

def merge_ranges(norm_ranges,n_ranges):
    
    """
    This function generated a dataframe of concat. ranges if both R and N type regions are present in the .fa
    """
    
    merged_ranges = pd.concat([norm_ranges, n_ranges])
    merged_ranges.index = range(len(merged_ranges.index))
    sorted_df = merged_ranges.sort_values(by=['start'], ascending=True)
    reorder_df=sorted_df.set_index(np.arange(len(sorted_df.index)))
        
    return reorder_df

##############################################################################
    
def create_df_ranges(ranges,n_ranges,n_bases_start):
    
    """
    This function is used to call on the compute ranges function and handle ranges if N bases are present in the .fa
    """

    #compute a dataframe for normal regions that are not N
    
    probes_ranges=compute_ranges(ranges,"R")
    
    #compute the dataframe for each respective range if N's are found in the file
    
    if n_bases_start:
        n_ranges_df=compute_ranges(n_ranges,"N")
        
        probes_ranges=merge_ranges(probes_ranges,n_ranges_df)
            
    #return that dataframe
    return probes_ranges
    
##############################################################################

def subtract_kmer_length(probes_ranges):
    
    """
    This function will take each row before each N type and deducts from it's end value the length of the kmer, all R regions are then put in a new DF
    """
    #for each row contianing an N
    for index, row in probes_ranges.iloc[1:].iterrows():
        if row['type']=="N":
            prev=probes_ranges.index.get_loc(index-1)
            parse_rows=(probes_ranges.iloc[prev])
            #take into account the different, the start of the next k-mer
            probes_ranges.loc[prev,['end']]=parse_rows["end"] - int(MERLENGTH)

    #subset all the R's into a seperate dataframe
    normal_ranges=probes_ranges.loc[(probes_ranges.type == "R")]
    
    return normal_ranges

##############################################################################
    
def generate_index_file(normal_ranges):
    """
    This function will take the two columns from the normal ranges DF and then append these ranges to a list. This list is then wri
tten out as a file
    """
    
    kmer_indices=[]

    #make the dataframes into lists
    start_ranges_list = normal_ranges['start'].tolist()
    end_ranges_list = normal_ranges['end'].tolist()

    #then zip the two lists together
    for start,end in zip(start_ranges_list,end_ranges_list):
        for x in range (start, end+1):
            kmer_indices.append(x)

    with open(out_index,"w") as k_file:
        for i in kmer_indices:
            k_file.write(str(i) + "\n")

    return out_index,kmer_indices

##############################################################################
    
def main():
    
    global USAGE 
    
    chr_name_jf_out = jf_query(jf_indexfile,fasta_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    bases_dist_start,n_bases_start = map_coords(fasta_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    ranges,n_ranges=group_ranges(bases_dist_start,n_bases_start)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    probes_ranges = create_df_ranges(ranges,n_ranges,n_bases_start)

    print("---%s seconds ---"%(time.time()-start_time))

    normal_ranges=subtract_kmer_length(probes_ranges)
    
    print("---%s seconds ---"%(time.time()-start_time))

    out_index,kmer_indices=generate_index_file(normal_ranges)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    print("Done")

if __name__== "__main__":
    main()
