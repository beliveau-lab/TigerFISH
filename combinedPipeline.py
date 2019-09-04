"""
Author: Robin Aguilar
Create date: 07/06/2019, Python3.7
Date last modified: 08/15/2019
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
from collections import OrderedDict
import collections

#modules needed to run refactoredBlockParse
import refactoredBlockparse as bp
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
from scipy import signal

#declare a timer
start_time=time.time()


#define command line usage
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

#Import user-specified command line values.
args = userInput.parse_args()
test_fasta = args.fasta_file
scaffold_fasta=args.scaffold_fasta
index = args.jf_indexfile
chrom= args.chr_name
SPAN = args.span_length
THRESHOLD = args.threshold
COMPOSITION = args.composition_score

########################################################################################################################

"""
Write a function that will run jellyfish and then will later remove the query file after the script is run
"""
#this function reads the jellyfish index file (hg38_index) and the fasta sequence to query 

def jf_query(jf_index,fa_query):
    query_file=subprocess.call(['jellyfish', 'query', jf_index, '-s',
                     fa_query, '-o', chr_name_jf_out], stderr=None, shell=False)
    return chr_name_jf_out

########################################################################################################################
"""
Write a function that will handle the jellyfish and fasta processing to find map repeat coordinates. 
"""

#the jf query file from the jf_query function in addition to the fasta being queried is used in this function
def map_coords(fa_file,jf_file):
    fa_seq = list(SeqIO.parse(open(fa_file),'fasta'))
    for fasta in fa_seq:
        sequence=str(fasta.seq).lower()
        
    #will append the k_mers and the corresponding counts to seperate lists
    k_mer=[]
    count=[]
    with open(jf_file, "r") as jf:
        for line in jf:
            k_mer.append(line.replace(" ","\t").split()[0])
            count.append(line.replace(" ","\t").split()[1])

    #create a list that is dependent on the values in count according to the threshold
    success_list=[]
    for i in count:
        if int(i)>=THRESHOLD:
            success_list.append(1)
        else:
            success_list.append(0)


    #make a list that will store the indices of each k-mer count
    indexList = []
    for i in range(0, (len(count))):
        indexList.append(i)

    #Contain an array of where all the successes are located
    iter_data = np.array(success_list)

    #iterate through each of the successes in a sum fashion using convolve, used to identify elevated k-mer enrich. regions
    iter_vals_convolve=np.convolve(iter_data,np.ones(SPAN,dtype=int),'valid')

    iter_list=[]
    [iter_list.append(float(element)) for element in iter_vals_convolve]

    #creates a DF from the convolve function and calls the column iter_sum (to represent the convolve iterative sum)
    iterative_sum=pd.DataFrame(list(iter_list),columns=['iter_sum'])

    #keeps track of the span of where the k-mer starts
    iterative_sum['row_span_start']=np.arange(len(iterative_sum))

    #keeps track of where the k-mer stops
    iterative_sum['row_span_end']=iterative_sum['row_span_start']+SPAN-1

    #keeps track of elevated sums that match the composition criteria
    passing_range_iter = iterative_sum.loc[(iterative_sum['iter_sum']/SPAN>=COMPOSITION)]

    #make two lists to take the min ranges and max ranges
    min_ranges=[]
    max_ranges=[]

    #collapse continuous regions by their indices
    for k,g in passing_range_iter.groupby(passing_range_iter['row_span_start'] - np.arange(passing_range_iter.shape[0])):
            min_ranges.append(min(g['row_span_start']))
            max_ranges.append(max(g['row_span_end']))

    #condenses dataframe into continuous index regions
    indices_to_parse = pd.DataFrame(list(zip(min_ranges, max_ranges)), columns =['start_index_range', 'end_index_range'])

    #now you need to return the nucleotide slice at which the ranges start and end
    sequence_start=[]
    sequence_end=[]

    #this will find the indices of that specific kmer that you're looking for (i.e the start and end)
    kmer_start_seq=[]
    start_list=[]
    for item in indices_to_parse['start_index_range']:
        start=k_mer[item]
        kmer_start_seq.append(start)
    for item in kmer_start_seq:
        start_list.append(sequence.lower().find(str(item).lower()))

    #do the same for the end indices as you did above
    kmer_end_seq=[]
    end_list=[]
    for item in indices_to_parse['end_index_range']:
        end=k_mer[item]
        kmer_end_seq.append(end)
    for item in kmer_end_seq:
        end_list.append(sequence.lower().rfind(str(item).lower())+len(k_mer[0]))

    #zip the two lists together and this should get you the start and end nucleotides of where that 18mer pattern was
    nucleotide_range = pd.DataFrame(list(zip(start_list, end_list)), columns =['start', 'end'])

    #let's add a column to the dataframe that will contain the chromosome information
    nucleotide_range['chr']=chrom

    #let's now merge by the chromosome to collapse overlapping regions
    collapsed_nucleotide_range=nucleotide_range.groupby((nucleotide_range.end.shift() - nucleotide_range.start).lt(0).cumsum()).agg(
        {'start': 'first','chr': 'first', 'end': 'last'})

    collapsed_nucleotide_range=collapsed_nucleotide_range[["chr","start","end"]]

    collapsed_nucleotide_range.to_csv(bed_file, header=None, index=None, sep='\t')

    #call bedtools to generate a fasta from the sequence you have generated
    subprocess.call(['bedtools', 'getfasta', '-fi', scaffold_fasta, '-bed', bed_file, '-fo', out_fasta], stderr=None, shell=False)

    return out_fasta

########################################################################################################################
"""
Write a function that will now run blockParse as a module to output the df of probes designed against the regions.
"""

#this function takes the fasta that you returned from the last function
def blockParse_run(output_fasta):
    fasta_sequences = list(SeqIO.parse(open(output_fasta),'fasta'))
    #you want to parse each fasta sequence as a string and append to a sequence list
    for fasta in fasta_sequences:
        name_list.append(fasta.id)
        sequence=(str(fasta.seq))
        sequence_list.append(sequence)
    #zip the names (headers of the fasta) and the string sequence into a list, then make into a dict
    zipped_list=zip(name_list,sequence_list)
    dict_name_seq=dict(zipped_list)
    #then run each item in the dict into the refactored blockParse script which can now handle multi-lined fastas
    for names,sequences in dict_name_seq.items():   
        print(bp.runSequenceCrawler(sequences,names, 36, 41, 20, 80,mt.DNA_NN3, 
                                      42, 47, 'AAAAA,TTTTT,CCCCC,GGGGG', 390, 50, 0,25, 
                                      25, None, True ,False, False,True, False, False))
            
########################################################################################################################
# MAIN
########################################################################################################################

def main():    
    global USAGE
    
    returned_jf_file = jf_query(index,test_fasta)  
        
    #the names of output files to be generated
    chr_name_jf_out=chrom+"_jf_temp.txt"
    bed_file=chrom + "_regions.bed"
    out_fasta=chrom + "_regions.fa"

    #a list for the regions and scaffold sequences 
    name_list=[]
    sequence_list=[]
    zipped_list=[]
 
    mapped_fasta=map_coords(test_fasta,returned_jf_file)
    
    blockParse_run(mapped_fasta)
    
    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")
    
if __name__== "__main__":
    main()


