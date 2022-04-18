#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##############################################################################

"""
Created on Fri Jun 25 12:40:14 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#specific script name
script_name = "repeat_ID"

#import libraries
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

##############################################################################

def test_open_index_file(index_file):
    """
    This function will take the two columns from the normal ranges DF and then 
    append these ranges to a list. This list is then written out as a file    
    Parameters
    ----------
    index_file : file
        file containing the index location of ATCG bases in the
        genome
        
    Returns
    -------
    kmer_indices : list
        file is read into a list

    Assertions tested
    ----------------
    - Opens valid index file
    """
    
    #open index file and add to list
    with open(index_file, "r") as idx_f:
        kmer_indices = [line.rstrip('\n') for line in idx_f]
    
    if len(kmer_indices) == 0:
        print("Invalid kmer file, exiting ...")
        exit()

    #check that the length if not equal to 0
    assert(len(kmer_indices) != 0)

    return kmer_indices

##############################################################################

def test_generate_kmer_count_lists(jf_count):
    """
    This function will store the values from the jf query file into a list for 
    the k_mers and the corresponding count
    
    Parameters
    ----------
    jf_count : file
        file derived from jellyfish count providing count of k-mers in whole
        genome
    Returns
    -------
    k_mer : list
        list containing specific k-mers
    count : list
        list containing k-mer count values

    Assertions tested
    -----------------
    - Valid strings are seperated into k_mer and count lists        
    """
    
    k_mer=[]
    count=[]
    
    with open(jf_count, "r") as jf:
        for line in jf:
            k_mer.append(line.replace(" ","\t").split()[0])
            count.append(line.replace(" ","\t").split()[1])

    #checks if k_mer list contains items
    assert(len(k_mer) != 0)

    #checks if count list contains items
    assert(len(count) != 0)

    return k_mer,count

##############################################################################

def test_check_threshold():
    """
    This function will generate a binary list based on whether the count value 
    passes the defined threshold
    
    Parameters
    ----------
    count : list
        list containing k-mer count values
    THRESHOLD : user argument (int)
        min val of k-mer count to be flagged as enrchiched
    Returns
    -------
    pass_thresh_list : list
        binary list where (1) means k-mer count value is >= threshold,
        else assigned 0

    Assertions tested
    -----------------
    - Tests logic if a list is correctly binarized based on threshold.
    """
    
    THRESHOLD = 5
    count = [1, 1, 1, 10, 12, 15, 20, 2, 4, 3]
    expected_output = [0, 0, 0, 1, 1, 1, 1, 0, 0, 0]

    pass_thresh_list=[]
    
    for i in count:
        if int(i)>=THRESHOLD:
            pass_thresh_list.append(1)
        else:
            pass_thresh_list.append(0)

    assert(expected_output == pass_thresh_list)

    return pass_thresh_list

##############################################################################

def test_convolve_successes():
    """
    The purpose of this function is to convolve over the array in the size of 
    the window you define, then generate a dataframe of all passing ranges   
    Parameters
    ----------
    pass_thresh_list : list
        binary list where (1) means k-mer count value is >= threshold,
        else assigned 0
    WINDOW : user argument (int)
        the length of the window to be searched in the pass_thresh_list
    COMPOSITION : user argument (float)
        the percentage of k-mers in the list that must be flagged as >=
        composition score
    Returns
    -------
    passing_range_iter : TYPE
        DESCRIPTION

    Assertions tested
    -----------------
    - Checks that array of data is appropriate sums
    - See if array appears as list
    - Map location of sliding window
    - Composition score selects windows
    """

    pass_thresh_list = [0, 0, 0, 1, 1, 1, 1, 0, 0, 0]

    convolve_check = [0, 1, 2, 3, 3, 2, 1, 0]

    WINDOW = 3

    COMPOSITION = 0.5

    #convert list of binary k-mer count vals into array
    iter_data = np.array(pass_thresh_list)
 
    #iterate through array over length of specific WINDOW
    iter_vals_convolve=np.convolve(iter_data,np.ones(WINDOW,dtype=int),
                                   'valid')

    #checks if convolve sum was conducted appropriately over window size
    assert(list(iter_vals_convolve) == convolve_check)

    #list of values in sum values in sliding window casted as float
    iter_list=[(float(element)) for element in iter_vals_convolve]

    #this is the sum of all the counts of the k-mers
    #the most important column here is the row_span_start,
    #this will tell you which k-mer it is
    #this casts to dataframe
    w_sum=pd.DataFrame(list(iter_list),columns=['iter_sum'])

    #the length of the list represents the start values
    w_sum['sliding_win_start']=np.arange(len(w_sum))

    #the location of the end of each sliding window
    w_sum['sliding_win_end']=w_sum['sliding_win_start']+WINDOW-1
    
    #these are all the k-mers that pass the range you specified
    #note that these are not their true indices in sequence 
    #(basically which number k-mer in order they appear)
    pass_w = w_sum.loc[(w_sum['iter_sum']/WINDOW>=COMPOSITION)]

    #checks logic for values that pass composition score
    comp_check = []
    comp_idx = []
    for i in range (len(iter_list)):
        if iter_list[i]/WINDOW>=COMPOSITION:
            comp_check.append(iter_list[i])
            comp_idx.append(i)

    #checks takes the sums and start values that pass
    pass_w_sum_list = pass_w['iter_sum'].tolist()
    pass_w_start_list = pass_w['sliding_win_start'].tolist()

    #checks if the passed sums match both checks
    assert(comp_check == pass_w_sum_list)

    #checks if the indices of passed sum vals also match
    assert(comp_idx == pass_w_start_list)

    return pass_w

##############################################################################

def test_obtain_repeat_indices():
    """
    This function will collapse the dataframe into continuous ranges, and then 
    append those to a seperate dataframe, then stored as seperate lists
    
    Parameters
    ----------
    pass_w : dataframe
        subset of the dataframe containing only the rows that pass the 
        composition score from surveying the sliding window
    Returns
    -------
    indices_to_parse : dataframe

    Assertions tested
    -----------------
    - Collapses dataframe into expected ranges 
        
    """

    pass_w_list = [2.0, 3.0, 3.0, 2.0]
    sliding_window_start = [2,3,4,5]
    sliding_window_end = [4,5,6,7]

    #generates expected to test logic
    #takes first of start and last of end since continuous
    expected_start = [2]
    expected_end = [7]

    pass_w = pd.DataFrame(
        {'iter_sum': pass_w_list,
         'sliding_win_start': sliding_window_start,
         'sliding_win_end': sliding_window_end
        })

    #lists for the continuous ranges of start and end values that will
    #be used to identify the repeat regions
    min_ranges=[]
    max_ranges=[]

    start_index=[]
    end_index=[]
    
    for k,g in pass_w.groupby(pass_w['sliding_win_start'] -
                              np.arange(pass_w.shape[0])):
            min_ranges.append(min(g['sliding_win_start']))
            max_ranges.append(max(g['sliding_win_end']))

    #add these two values into a dataframe where we can scan the
    #indices of jellyfish count
    indices_to_parse = pd.DataFrame(list(zip(min_ranges, max_ranges)),
                                    columns =['start_index_range',
                                              'end_index_range'])

    start_index = indices_to_parse['start_index_range'].tolist()
    end_index = indices_to_parse['end_index_range'].tolist()

    #checks that logic above sucessfully collapses start and end
    assert(start_index == expected_start)
    assert(end_index == expected_end)

    return indices_to_parse

##############################################################################

def test_find_repeats():
    """
    This function will refer to the indices to parse and append to a list, 
    then at that index, the start val in kmer_indices is appended, this is 
    also done for end sequences.
    
    Parameters
    ----------
    indices_to_parse : datafrane
        contains continous start and end coordinates
    kmer_indices : list
        contains the index derived from jellyfish files
        
    Returns
    -------
    kmer_start_seq : list
        list containing start coords of k-mer
    kmer_end_seq : list
        list containing ned coords of k-mer
    index_start : list
        list of where start indices were noted from parsing sliding window
    index_end : list
        list of where end indices were noted from parsing sliding window

    Assertions tested
    -----------------
    - Checks that the end val (last mer base) accounts for k-mer size.
    """

    #establish test lists to make a dataframe to test assertion
    expected_start = [2]
    expected_end = [7]

    indices_to_parse = pd.DataFrame(
        {'start_index_range': expected_start,
         'end_index_range': expected_end
        })

    #produce the kmer_indices list
    kmer_indices = [0,1,2,3,4,5,6,7,8,9]

    #assigned merlength
    MERLENGTH = 18

    #list of indices of all valid ATCG bases that represent k-mer location
    kmer_start=[]
    kmer_end=[]

    #index of regions found as valid in k-mer sliding window
    index_start=[]
    index_end=[]
    
    #obtains start ranges for k-mers
    for item in indices_to_parse['start_index_range']:
        index_start.append(int(item))
        start=kmer_indices[item]
        kmer_start.append(int(start))

    #generates start and end vals
    repeat_starts_df=pd.DataFrame(list(zip(kmer_start,index_start)),
                                  columns=['r_start','r_index_start'])

    #generates end that accounts for MERLENGTH
    for item in indices_to_parse['end_index_range']:
        index_end.append(int(item))
        end=kmer_indices[item]
        kmer_end.append(int(int(end) + MERLENGTH))
        
    repeat_ends_df=pd.DataFrame(list(zip(kmer_end,index_end)),
                                columns=['r_end','r_index_end'])

    #makes repeat indices
    repeat_indices=pd.concat([repeat_starts_df, repeat_ends_df], axis=1)

    assert(len(repeat_indices) != 0)

    return kmer_start,kmer_end,index_start,index_end

##############################################################################

def test_nucleotide_range(bed_file,chrom):
    """
    This function will take the start and end lists and then map this back 
    to the regions in sequence
    
    Parameters
    ----------
    kmer_start : list
        location of start of k-mer
    kmer_end : list
        location of end of k-mer
    index_start : list
        start index of k-mer sliding window
    index_end : list
        end index of k-mer sliding window
    bed_file : file
        path and name of output bed file
    chrom : user arg string
        given scaffold to report bed file
    Returns
    -------
    None. The output bed file is written at specified path

    Assertions tested
    -----------------
    - Checks if collapse range features works appropriately
    - Validates BED file created
    """

    #adds additional repeat to what was generated
    #here 30 is the largest val continuous and 2 is the smallest
    kmer_start = [2,4]
    kmer_end = [25,30]
    index_start = [2,5]
    index_end = [7,10]

    #fill expected dataframe to test collapsed behavior
    expected_chrom = ["chrX"]
    expected_start = [2]
    expected_end = [30]

    expected_df = pd.DataFrame(
        {'chr': expected_chrom,
         'r_start': expected_start,
         'r_end': expected_end
        })

    #zip the two lists together and this should get you the start and end 
    #nucleotides of where that 18mer pattern was
    seq_range = pd.DataFrame(list(zip(kmer_start, kmer_end)), 
                                    columns =['r_start', 'r_end'])

    print(seq_range)

    seq_index_range = pd.DataFrame(list(zip(index_start,index_end)),
                                   columns =['r_index_start','r_index_end'])

    print(seq_index_range)

    #let's add a column to the dataframe that will contain the
    #chromosome information
    seq_range['chr']=chrom
    seq_index_range['chr']=chrom

    #let's now merge by the chromosome to collapse overlapping regions
    collapsed_nucleotide_range=seq_range.groupby((seq_range.r_end.shift() - 
                                                  seq_range.r_start).lt(0).cumsum()).agg(
        {'r_start': 'first','chr': 'first', 'r_end': 'last'})

    collapsed_nucleotide_range=collapsed_nucleotide_range[["chr","r_start","r_end"]]

    #checks that collapsed function works
    assert(collapsed_nucleotide_range.equals(expected_df))

    #generates bed file in tab seperated format
    collapsed_nucleotide_range.to_csv(bed_file, header=None, index=None, sep='\t')

##############################################################################







