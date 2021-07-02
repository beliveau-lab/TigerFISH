"""
Robin Aguilar
Beliveau and Noble labs
Date Modified: 03/01/2020
Modification: To make some functions more efficient and commented with user statements 
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

# Define the command line usage.
USAGE = """Usage: repeat_ID.py -f fasta_file -j jf_indexfile -c chr_name -schr scaffold_fasta -st start [options]

  Options:

    -s     The length of the scanning window to search for k-mer enriched regions, default = 3000

    -t     The minimum number of counts that defines a k-mer enriched region, default = 10

    -c     The minimum proportion of k-mers to consider a k-mer enriched region significant, default = 0.50
   

"""

#write arguments so users can specify the commands they would like for each variable
userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing or optionally can be made into a  '
        '.bed file can be outputted instead if \'-b\' is flagged. Tm values '
        'are corrected for [Na+] and [formamide].')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-j', '--jf_count', action='store', required=True,
                               help='The jf file of a given genome')
requiredNamed.add_argument('-i', '--index_file', action='store', required=True,
                               help='The index file of a given genome')
requiredNamed.add_argument('-chr', '--chr_name', action='store', required=True,
                               help='Define the scaffold being queried')
requiredNamed.add_argument('-st','--start',action='store',required=True,
                           help='The start sequence of the fasta region if not starting at beginning of scaffold; default is ' '0')
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
requiredNamed.add_argument('-schr', '--scaffold_fasta', action='store',required=True,
                           help='Used to generate fasta file of kmer rich regions; default is ')

requiredNamed.add_argument('-o_b', '--bed_file', action='store',required=True,
                           help='Used to generate fasta file of kmer rich regions; default is ')

#Import user-specified command line values.
args = userInput.parse_args()
index_file = args.index_file
jf_count = args.jf_count
chrom= args.chr_name
scaffold_fasta=args.scaffold_fasta
SPAN = args.span_length
THRESHOLD = args.threshold
COMPOSITION = args.composition_score
START=args.start
bed_file = args.bed_file

#a list for the regions and scaffold sequences 
name_list=[]
sequence_list=[]
zipped_list=[]
bases_dist_start=[]
MERLENGTH=18

##############################################################################

def open_index_file(index_file):
    """
    This function will take the two columns from the normal ranges DF and then 
    append these ranges to a list. This list is then written out as a file
    """
    
    with open(index_file, "r") as idx_f:
        kmer_indices = [line.rstrip('\n') for line in idx_f]

    return kmer_indices

##############################################################################

def generate_kmer_count_lists(jf_count):
    
    """
    This function will store the values from the jf query file into a list for 
    the k_mers and the corresponding count
    """
    
    k_mer=[]
    count=[]
    float_count=[]
    
    with open(jf_count, "r") as jf:
        for line in jf:
            k_mer.append(line.replace(" ","\t").split()[0])
            count.append(line.replace(" ","\t").split()[1])

    #This list is later used in the process of calling the max_count for the repeat region for each probe
    float_count=list(map(float,count))

    return k_mer,count,float_count

##############################################################################

def success_l(count):
    
    """
    This function will generate a binary list based on whether the count value 
    passes the defined threshold
    """
    
    success_list=[]
    
    for i in count:
        if int(i)>=THRESHOLD:
            success_list.append(1)
        else:
            success_list.append(0)

    return success_list

##############################################################################

def convolve_successes(success_list):
    """
    The purpose of this function is to convolve over the array in the size of 
    the window you define, then generate a dataframe of all passing ranges
    """
    #Contain an array of where all the successes are located
    iter_data = np.array(success_list)

    #iterate through each of the successes in a sum fashion
    iter_vals_convolve=np.convolve(iter_data,np.ones(SPAN,dtype=int),'valid')

    iter_list=[]
    [iter_list.append(float(element)) for element in iter_vals_convolve]

    #this is the sum of all the counts of the k-mers
    #the most important column here is the row_span_start, this will tell you which k-mer it is
    iterative_sum=pd.DataFrame(list(iter_list),columns=['iter_sum'])

    iterative_sum['row_span_start']=np.arange(len(iterative_sum))

    iterative_sum['row_span_end']=iterative_sum['row_span_start']+SPAN-1
    #these are all the k-mers that pass the range you specified
    #note that these are not their true indices in sequence (basically which number k-mer in order they appear)
    passing_range_iter = iterative_sum.loc[(iterative_sum['iter_sum']/SPAN>=COMPOSITION)]

    return passing_range_iter

##############################################################################

def obtain_repeat_indices(passing_range_iter):
    
    """
    This function will collapse the dataframe into continuous ranges, and then 
    append those to a seperate dataframe, then stored as seperate lists
    """
    
    min_ranges=[]
    max_ranges=[]

    start_index=[]
    end_index=[]
    
    for k,g in passing_range_iter.groupby(passing_range_iter['row_span_start'] - np.arange(passing_range_iter.shape[0])):
            min_ranges.append(min(g['row_span_start']))
            max_ranges.append(max(g['row_span_end']))

    #add these two values into an dataframe where we can scan the indices of jellyfish count
    indices_to_parse = pd.DataFrame(list(zip(min_ranges, max_ranges)), columns =['start_index_range', 'end_index_range'])

    start_index = indices_to_parse['start_index_range'].tolist()
    end_index = indices_to_parse['end_index_range'].tolist()

    return indices_to_parse

##############################################################################

def repeat_indices_in_seq(indices,kmer_indices,k_mer):
    
    """
    This function will refer to the indices to parse and append to a list, 
    then at that index, the start val in kmer_indices is appended, this is 
    also done for end sequences.
    """

    sequence_start=[]
    sequence_end=[]

    kmer_start_seq=[]
    start_list=[]

    kmer_end_seq=[]
    end_list=[]

    true_index_start=[]
    true_index_end=[]
    
    for item in indices['start_index_range']:
        true_index_start.append(int(item))
        start=kmer_indices[item]
        kmer_start_seq.append(int(start))

    repeat_starts_df=pd.DataFrame(list(zip(kmer_start_seq,true_index_start)),columns=['r_start','r_index_start'])

    for item in indices['end_index_range']:
        true_index_end.append(int(item))
        end=kmer_indices[item]
        kmer_end_seq.append(int(int(end) + len(str(k_mer[0]))))

    repeat_ends_df=pd.DataFrame(list(zip(kmer_end_seq,true_index_end)),columns=['r_end','r_index_end'])

    repeat_indices=pd.concat([repeat_starts_df, repeat_ends_df], axis=1)

    return kmer_start_seq,kmer_end_seq,true_index_start,true_index_end

##############################################################################

def nucleotide_range(start_l,end_l,start_index_l,end_index_l,bed_file):
    """
    This function will take the start and end lists and then map this back 
    to the regions in sequence
    """

    #zip the two lists together and this should get you the start and end nucleotides of where that 18mer pattern was
    nucleotide_range = pd.DataFrame(list(zip(start_l, end_l)), columns =['r_start', 'r_end'])

    nucleotide_index_range = pd.DataFrame(list(zip(start_index_l, end_index_l)), columns =['r_index_start', 'r_index_end'])

    #let's add a column to the dataframe that will contain the chromosome information
    nucleotide_range['chr']=chrom
    nucleotide_index_range['chr']=chrom

    #let's now merge by the chromosome to collapse overlapping regions
    collapsed_nucleotide_range=nucleotide_range.groupby((nucleotide_range.r_end.shift() - nucleotide_range.r_start).lt(0).cumsum()).agg(
        {'r_start': 'first','chr': 'first', 'r_end': 'last'})

    collapsed_nucleotide_index_range=nucleotide_index_range.groupby((nucleotide_index_range.r_index_end.shift() - nucleotide_index_range.r_index_start).lt(0).cumsum()).agg(
        {'r_index_start': 'first','chr': 'first', 'r_index_end': 'last'})

    #collapsed_nucleotide_range.columns=["chr","start","end"]
    collapsed_nucleotide_range=collapsed_nucleotide_range[["chr","r_start","r_end"]]

    collapsed_nucleotide_index_range=collapsed_nucleotide_index_range[['chr','r_index_start','r_index_end']]

    repeat_indices=pd.concat([collapsed_nucleotide_range, collapsed_nucleotide_index_range], axis=1)
    #repeat_indices.to_csv("repeat_indices.bed", header=None, index=None, sep='\t')

    collapsed_nucleotide_range.to_csv(bed_file, header=None, index=None, sep='\t')

    return repeat_indices

##############################################################################
# MAIN
##############################################################################

def main():
    
    global USAGE 

    #the names of output files to be generated
    
    kmer_indices=open_index_file(index_file)
    
    print("---%s seconds ---"%(time.time()-start_time))

    k_mer,count,float_count=generate_kmer_count_lists(jf_count)

    print("---%s seconds ---"%(time.time()-start_time))

    success_list=success_l(count)

    print("---%s seconds ---"%(time.time()-start_time))

    passing_range_iter=convolve_successes(success_list)

    print("---%s seconds ---"%(time.time()-start_time))

    indices_to_parse=obtain_repeat_indices(passing_range_iter)

    print("---%s seconds ---"%(time.time()-start_time))

    kmer_start_seq,kmer_end_seq,true_index_start,true_index_end=repeat_indices_in_seq(indices_to_parse,kmer_indices,k_mer)

    print("---%s seconds ---"%(time.time()-start_time))

    repeat_indices=nucleotide_range(kmer_start_seq,kmer_end_seq,true_index_start,true_index_end,bed_file)

    print("---%s seconds ---"%(time.time()-start_time))
    
    print("Done")

if __name__== "__main__":
    main()
