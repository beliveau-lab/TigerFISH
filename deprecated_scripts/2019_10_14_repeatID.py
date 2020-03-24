"""
Robin Aguilar
Beliveau and Noble Labs
Dependencies: Jellyfish, Bedtools
Input: Fasta for region of interest, JF index of genome of interest
Output: A .bed file of the elevated kmer regions, a .fa of said regions,
and a DF containing information about all the probes generated in those
regions. 
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
userInput.add_argument('-s', '--span_length', action='store', default=3000,
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
requiredNamed.add_argument('-st','--start',action='store',required=True,
                           help='The start sequence of the fasta region if not starting at beginning of scaffold; default is ' '0')

#Import user-specified command line values.
args = userInput.parse_args()
test_fasta = args.fasta_file
scaffold_fasta=args.scaffold_fasta
index = args.jf_indexfile
chrom= args.chr_name
SPAN = args.span_length
THRESHOLD = args.threshold
COMPOSITION = args.composition_score
START=args.start

#the names of output files to be generated
chr_name_jf_out=chrom+"_jf_temp.txt"
out_index=chrom+"_index.txt"
bed_file=chrom + "_regions.bed"
out_fasta=chrom + "_regions.fa"


#a list for the regions and scaffold sequences 
name_list=[]
sequence_list=[]
zipped_list=[]
bases_dist_start=[]
MERLENGTH=18

def jf_query(jf_index,fa_query):
    """
    Runs jellyfish to generate a query file from the jellyfish index provided
    """
    query_file=subprocess.call(['jellyfish', 'query', jf_index, '-s',
                     fa_query, '-o', chr_name_jf_out], stderr=None, shell=False)
    return chr_name_jf_out

#these are lists that will be used to find the indices of the regions that contain N's and non-N bases 
bases_dist_start=[]
n_bases_start=[]

def map_coords(fa_file):
    """
    This function will take a the fasta file you pass it and identify where all N and non-N bases are located
    """
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

#these are the lists that will hold the continuous ranges of N and non-N bases
ranges = []
n_ranges=[]

def group_ranges(bases_start,n_bases):
    """
    This function will return continous ranges of N and non-N bases.
    """
    #this collapses the set of the regular bases to provide you with ranges in a list
    for k, g in groupby(enumerate(bases_start), lambda x: x[1]-x[0]):
        group = list(map(itemgetter(1), g))
        ranges.append(str(group[0]) + "\t" + str(group[-1]+1))

    #this collapses the set of N bases to provide you with ranges in a list 
    for k, g in groupby(enumerate(n_bases), lambda x: x[1]-x[0]):
        group=list(map(itemgetter(1),g))
        n_ranges.append(str(group[0]) + "\t" + str(group[-1]+1))

    return ranges,n_ranges

def df_ranges(ranges_list,n_ranges_list):
    """
    This function will generate a dataframe with the start and end ranges of non-N bases
    """
    #dataframe with the regular bases
    normal_ranges=pd.DataFrame(list(ranges_list), columns=['ranges'])
    normal_ranges[['start','end']] = normal_ranges.ranges.str.split(expand=True)
    num_start=normal_ranges['start'].astype(int)
    num_end=normal_ranges['end'].astype(int)
    num_ranges = pd.concat([num_start, num_end], axis=1)
    #the R represents that the column is composed of regular bases
    num_ranges['type']="R"
    num_ranges['size']=num_ranges['end'].astype(int)-num_ranges['start'].astype(int)

    #dataframe for the n bases
    n_ranges=pd.DataFrame(list(n_ranges_list), columns=['ranges'])
    n_ranges[['start','end']] = n_ranges.ranges.str.split(expand=True)
    num_n_start=n_ranges['start'].astype(int)
    num_n_end=n_ranges['end'].astype(int)
    num_n_ranges = pd.concat([num_n_start, num_n_end], axis=1)
    #the N represents that these are N bases
    num_n_ranges['type']="N"
    num_n_ranges['size']=num_n_ranges['end'].astype(int)-num_n_ranges['start'].astype(int)

    return num_ranges,num_n_ranges

def merge_range_dfs(regular_ranges,n_regular_ranges):
    """
    This function merges the two N base and non-N base dataframes into one and sorts them accordingly in numerical order
    """ 
    merged_ranges = pd.concat([regular_ranges, n_regular_ranges])
    merged_ranges.index = range(len(merged_ranges.index))

    sorted_df = merged_ranges.sort_values(by=['start'], ascending=True)
    reorder_df=sorted_df.set_index(np.arange(len(sorted_df.index)))

    return reorder_df
def subtract_kmer_length(ordered_df):
    """
    This function will take each row before each N type and deducts from it's end value the length of the kmer, all R regions are then put in a new DF
    """
    for index, row in ordered_df.iloc[1:].iterrows():
        if row['type']=="N":
            prev=ordered_df.index.get_loc(index-1)
            parse_rows=(ordered_df.iloc[prev])
            ordered_df.loc[prev,['end']]=parse_rows["end"] - int(MERLENGTH)

    #subset all the R's into a seperate dataframe
    normal_ranges=ordered_df.loc[(ordered_df.type == "R")]

    return normal_ranges

kmer_indices=[]

def generate_index_file(normal_ranges_df):
    """
    This function will take the two columns from the normal ranges DF and then append these ranges to a list. This list is then written out as a file
    """
    #make the dataframes into lists
    start_ranges_list = normal_ranges_df['start'].tolist()
    end_ranges_list = normal_ranges_df['end'].tolist()

    #then zip the two lists together
    for start,end in zip(start_ranges_list,end_ranges_list):
        for x in range (start, end+1):
            kmer_indices.append(x)

    with open(out_index,"w") as k_file:
        for i in kmer_indices:
            k_file.write(str(i) + "\n")

    return out_index

k_mer=[]
count=[]
float_count=[]
def generate_kmer_count_lists(jf_file):
    """
    This function will store the values from the jf query file into a list for the k_mers and the corresponding count
    """ 
    with open(jf_file, "r") as jf:
        for line in jf:
            k_mer.append(line.replace(" ","\t").split()[0])
            count.append(line.replace(" ","\t").split()[1])

    #This list is later used in the process of calling the max_count for the repeat region for each probe
    float_count=list(map(float,count))

    return k_mer,count,float_count

success_list=[]

def success_l(count_l):
    """
    This function will generate a binary list based on whether the count value passes the defined threshold
    """
    for i in count_l:
        if int(i)>=THRESHOLD:
            success_list.append(1)
        else:
            success_list.append(0)

    return success_list

def convolve_successes(s_list):
    """
    The purpose of this function is to convolve over the array in the size of the window you define, then generate a dataframe of all passing ranges
    """
    #Contain an array of where all the successes are located
    iter_data = np.array(s_list)

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

min_ranges=[]
max_ranges=[]

start_index=[]

end_index=[]

def obtain_repeat_indices(pass_iter):
    """
    This function will collapse the dataframe into continuous ranges, and then append those to a seperate dataframe, then stored as seperate lists
    """
    for k,g in pass_iter.groupby(pass_iter['row_span_start'] - np.arange(pass_iter.shape[0])):
            min_ranges.append(min(g['row_span_start']))
            max_ranges.append(max(g['row_span_end']))

    #add these two values into an dataframe where we can scan the indices of jellyfish count
    indices_to_parse = pd.DataFrame(list(zip(min_ranges, max_ranges)), columns =['start_index_range', 'end_index_range'])

    start_index = indices_to_parse['start_index_range'].tolist()
    end_index = indices_to_parse['end_index_range'].tolist()

    return indices_to_parse
sequence_start=[]
sequence_end=[]

kmer_start_seq=[]
start_list=[]

kmer_end_seq=[]
end_list=[]

true_index_start=[]
true_index_end=[]

def repeat_indices_in_seq(indices):
    """
    This function will refer to the indices to parse and append to a list, then at that index, the start val in kmer_indices is appended, this is also done for end sequences.
    """

    for item in indices['start_index_range']:
        true_index_start.append(int(item))
        start=kmer_indices[item]
        kmer_start_seq.append(int(start))

    repeat_starts_df=pd.DataFrame(list(zip(kmer_start_seq,true_index_start)),columns=['r_start','r_index_start'])

    for item in indices['end_index_range']:
        true_index_end.append(int(item))
        end=kmer_indices[item]
        kmer_end_seq.append(int(end+len(k_mer[0])))

    repeat_ends_df=pd.DataFrame(list(zip(kmer_end_seq,true_index_end)),columns=['r_end','r_index_end'])

    repeat_indices=pd.concat([repeat_starts_df, repeat_ends_df], axis=1)

    return kmer_start_seq,kmer_end_seq,true_index_start,true_index_end
def nucleotide_range(start_l,end_l,start_index_l,end_index_l):
    """
    This function will take the start and end lists and then map this back to the regions in sequence
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
    repeat_indices.to_csv("repeat_indices.bed", header=None, index=None, sep='\t')

    collapsed_nucleotide_range.to_csv(bed_file, header=None, index=None, sep='\t')

    return repeat_indices

def generate_fasta():
    """
    This function will take the dataframe generated about as a bed file and then generate a fasta from those coords
    """

    #call bedtools to generate a fasta from the sequence you have generated
    subprocess.call(['bedtools', 'getfasta', '-fi', scaffold_fasta, '-bed', bed_file, '-fo', out_fasta], stderr=None, shell=False)

    return out_fasta

def main():

    returned_jf_file = jf_query(index,test_fasta)

    print("---%s seconds ---"%(time.time()-start_time))

    bases_start,n_bases=map_coords(test_fasta)

    print("---%s seconds ---"%(time.time()-start_time))

    normal_ranges,non_normal_ranges=group_ranges(bases_start,n_bases)

    print("---%s seconds ---"%(time.time()-start_time))

    num_ranges,num_n=df_ranges(normal_ranges,non_normal_ranges)

    print("---%s seconds ---"%(time.time()-start_time))

    out_df=merge_range_dfs(num_ranges,num_n)

    print("---%s seconds ---"%(time.time()-start_time))

    range_df=subtract_kmer_length(out_df)

    print("---%s seconds ---"%(time.time()-start_time))

    index_file_lists=generate_index_file(range_df)

    print("---%s seconds ---"%(time.time()-start_time))

    kmer_list,kmer_list_count,kmer_list_float=generate_kmer_count_lists(returned_jf_file)

    print("---%s seconds ---"%(time.time()-start_time))

    to_convolve=success_l(kmer_list_count)

    print("---%s seconds ---"%(time.time()-start_time))

    iter_range=convolve_successes(to_convolve)

    print("---%s seconds ---"%(time.time()-start_time))

    indices_parse=obtain_repeat_indices(iter_range)

    print("---%s seconds ---"%(time.time()-start_time))

    kmer_start,kmer_end,kmer_i_start,kmer_i_end=repeat_indices_in_seq(indices_parse)

    print("---%s seconds ---"%(time.time()-start_time))

    print("---%s seconds ---"%(time.time()-start_time))

    nuc_range=nucleotide_range(kmer_start,kmer_end,kmer_i_start,kmer_i_end)

    print("---%s seconds ---"%(time.time()-start_time))

    generate_fasta()

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
