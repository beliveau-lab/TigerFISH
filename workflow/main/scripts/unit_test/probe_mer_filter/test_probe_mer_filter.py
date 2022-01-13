#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Mon Jun 28 15:47:32 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#specific script name
script_name = "probe_mer_filter"

#load libraries used
import time
import argparse
import pandas as pd
from itertools import groupby

##############################################################################

def test_read_region(file_path, enrich_score, copy_num):
    """
    
    Function reads in probe file with normalized filter scores that are
    pre-filter. This takes the ENRICH and COPY_NUM values to filter probes.
    
    Parameters
    ----------
    file_path : file
        probe file that contains pre-filter k_score and k_norm values.
    Returns
    -------
    region_df : dataframe
        filtered dataframe that contains the filtered probes based on user
        specified arguments

    Assertions tested
    -----------------
    - Checks that file is read and opened.
    - Checks that parameters filter file as expected.
    """
    
    #read in the file with the appropriate column names
    colnames = ["chrom","p_start","p_end","probe","Tm","region",
                "r_count_total","h_count_total","k_score","k_norm"]
    
    #read as a dataframe
    region_df = pd.read_csv(file_path, delimiter = '\t', names = colnames)
        
    region_df['k_score']=region_df['k_score'].astype(float)

    #filters dataframe based on user parameters
    region_df=region_df[(region_df['k_score'] >= float(enrich_score)) &
                        (region_df['r_count_total'] >= int(copy_num))]
    
    #sorts the pandas dataframe based on descending order
    region_df=region_df.sort_values(by='k_norm', ascending=False)

    #two probes were provided in test file
    #only one should remain based on filter specs
    assert(len(region_df) == 2)

    return region_df

##############################################################################

def test_split_mers(file_path, enrich_score, copy_num, merlength):
    """
    Implements the generate_kmers function to add a column to the probe
    dataframe that decomposes each probe into it's k-mers (length specified
    by user)'
    
    Parameters
    ----------
    region_df : dataframe
        dataframe that contains the probe sequences
    Returns
    -------
    region_df : dataframe
        dataframe that now contains a column that is a list of k-mers
        for the probe sequence
    Assertions tested
    -----------------
    - That k-mers are stored in a list in order as a column.
    """

    #read in the test dataframe
    region_df = test_read_region(file_path, enrich_score, copy_num)

    #make a list of the probes in the dataframe that have been filtered
    probe_list = region_df['probe'].tolist()

    #make a list for the k-mers
    all_mers_list = []
    
    #iterate over each probe in the list
    for probe in probe_list:
        #decompose the probe into k-mers of given length
        #function for kmers was already tested in prior scripts
        mer_list = generate_kmers(probe,merlength)
        #append list to list
        all_mers_list.append(mer_list)

    #all_mers_list should equal the number of probes provided
    assert(len(all_mers_list) == len(probe_list))
 
    #add list of lists to appropriate dataframe rows
    region_df['mers_list'] = all_mers_list

    return region_df

##############################################################################

def generate_kmers(sequence,k_size):
    """
    Function take a string sequence and splits it into k-mers of a specified
    size in the form of a list.
    Parameters
    ----------
    sequence : string
        sequence to be decomposed into k-mers
    k_size : int
        length of k-mers to be decomposed
    Returns
    -------
    kmers : list
        list of k-mers of specified size
    Assertions tested
    -----------------
    - Script was tested in kmer_filter. No assertions checked here.
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

def test_rm_shared_mer_probes(file_path, enrich_score, copy_num, merlength, 
        mer_cutoff):
    """
    For all repeat regions that have more than one probe in them, this
    function will take the list of k-mers and add the top ranking probe
    based on k_norm into a set. For all probes within the same repeat region,
    if, it shares at least >= specified mer_cutoff (user specified prop), 
    the probe is removed from the final probe set. This list column of each
    probes k-mers is then dropped.
    Parameters
    ----------
    region_df : dataframe
        contains relevant probe information
    Returns
    -------
    region_df : dataframe
        contains relevant probe information
    """
    
    region_df = test_split_mers(file_path, enrich_score, copy_num, merlength)

    #make a list for probes to be flagged that will be removed
    probe_to_remove = []
    
    #do a groupby over the repeat regions
    grouped_regions=region_df.groupby('region',sort=False)
    
    for name,group in grouped_regions:
        #if there is more than one probe in a given repeat region
        if len(group) > 1:
            
            #gets a list of all k-mers in each probe
            max_mers_list = group['mers_list'].tolist()
    
            #defines top ranking "keep set of 18-mers"
            keep_mers_set = set(max_mers_list[0])
 
            #get a list of all probes
            probes_list = group['probe'].tolist()
            
            for i in range(1,len(max_mers_list)):
                #will keep track of number of 18-mers within a given mers
                #list in the keep set 
                in_count = 0
                #iterates starting with max_mers_list[1]
                for mer in max_mers_list[i]:
                    #checks if a mer is in the keep set
                    if mer in keep_mers_set:
                        #if so tally count
                        in_count += 1
                    if mer not in keep_mers_set:
                        #if mer is unique add it to the keep set
                        keep_mers_set.add(mer)

                #take the lengh of the given max_mer_list[i]
                list_len = len(max_mers_list[i])

                #compute proportion of shared 18-mers
                #if proportion is >= cutoff, then you want to cull
                if float(in_count/list_len) >= float(mer_cutoff):
                    probe_to_remove.append(probes_list[i])
    
    #remove rows of the dataframe that have probes that should be filtered
    region_df = region_df[~region_df['probe'].isin(probe_to_remove)]

    #drop the mers column
    region_df.drop(['mers_list'], axis=1, inplace=True)
   
    #validates that if a particular prop of k-mers are shared
    #2/2 probes remains
    assert(len(region_df) == 2)

    return region_df

##############################################################################

def test_writee_file(file_path, enrich_score, copy_num, merlength,
        mer_cutoff,out_path):
    """
    This function writes the filtered probe dataframe into a file
    Parameters
    ----------
    region_df : dataframe
        dataframe containing filtered probes
    Returns
    -------
    None. Probe file is generated from dataframe

    Assertions tested
    -----------------
    - That file specified is written
    """

    region_df = test_rm_shared_mer_probes(file_path, enrich_score, copy_num,
         merlength, mer_cutoff)

    region_df.to_csv(str(out_path), header=False, index=False, sep="\t")

##############################################################################

