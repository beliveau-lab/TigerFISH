#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################

"""
Created on Mon Jun 28 11:09:32 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#specific script name
script_name = "kmer_filter"

#first you should load the libraries
import time
import argparse
import pandas as pd
from itertools import groupby
from collections import Counter
from Bio import SeqIO
from Bio.Seq import reverse_complement as rev_comp

##############################################################################

def test_read_probe_file(probe_file):
    """
    This function reads the probe file containing designed probes into a 
    pandas dataframe
    
    Parameters
    ----------
    probe_file : file
        file containing probes designed from the design_probes script
    Returns
    -------
    probe_data : dataframe
        dataframe of probes designed with labeled columns

    Assertions tested
    -----------------
    - Function reads probe file and is of expected length without duplicates
    """

    expected_length = 26

    #read probes csv
    colnames = ["chrom","p_start","p_stop","probe","Tm", "regions"]

    probe_df=pd.read_csv(probe_file, delimiter = '\t',names=colnames)
 
    probe_df = probe_df.drop_duplicates(subset=['probe'], keep='first')

    #checks if expected output matches test file
    assert(len(probe_df) == expected_length)

    return probe_df

##############################################################################

def test_read_fasta_dict(fasta_file):
    """
    This script reads the fasta file generated for all repeat regions that
    underwent probe design and adds them into a dictionary for easy region
    lookup.
    
    Parameters
    ----------
    fasta_file : fasta file
        multi entry fasta contains sequence of each repeat region, generated
        in the design probes script
    Returns
    -------
    seq_dict : dictionary
        contains the probe region header as key and sequence as value

    Assertions tested
    -----------------
    - Checks that seq_dict being generated is not empty
    """

    #parse out the probe file from blockparse
    seq_dict = {rec.id : rec.seq.upper() for rec in SeqIO.parse(fasta_file,
                                                                "fasta")}


    assert(len(seq_dict) != 0)

    return seq_dict

##############################################################################

def test_repeat_count(jf_file,fasta_file,probe_file):
    
    """
    This function takes the probes from probes data, and generates indep.
    dictionaries accountaing for the total sum of k-mer count is found 
    within a probes target repeat vs the entire human genome.
    
    Parameters
    ----------
    probe_data : dataframe
        dataframe of designed probes
    seq_dict : dictionary
        dictionary of repeat region and repeat sequence
    jf_file : file
        file containing jellyfish k-mer and k-mer counts
    Returns
    -------
    all_repeat_counts : dictionary
        stores the probe seq as key, and sum of all k-mer counts in repeat
        as val
    all_hg_counts : dictionary
        stores the probe seq as key, and sum of all k-mer counts in genome
        as val

    Assertions tested
    -----------------
    - Checks if jellyfish count file provided populates a dict
    - Validates that probe k-mer sums are generated correctly
    """    


    #establish testing parameters by obtaining needed files and vals
    MERLENGTH = 18
    seq_dict = test_read_fasta_dict(fasta_file)
    probe_df = test_read_probe_file(probe_file)

    #just take the first few probes to test
    probe_df = probe_df[:1]

    #make dictionary ot read jellyfish file into k-mer (key) and count (val)
    jf_dict = {}

    #read jellyfish file into two lists
    with open(jf_file) as jf:
        rows = (line.split() for line in jf)
        for row in rows:
            jf_dict[row[0]] = int(row[1])

    #checks that the dictionary is populated
    assert(len(jf_dict) != 0)

    #make dictionaries for repeat counts and genome counts
    #probe (key) and sum of all k-mer counts (val)
    all_hg_counts = {}
    all_repeat_counts = {}

    #group over each repeat region
    grouped_regions = probe_df.groupby(['regions'],sort = False)
    
    #iterate over each repeat region(name)
    for name,group in grouped_regions:
        
        #make a list of the probes in this group
        probes_list = group["probe"].tolist()
        
        #implement function to split region from the fasta dict(seq_dict)
        #into mers of a given size
        #region_mers is a list
        region_mers = generate_kmers(str(seq_dict[name]),int(MERLENGTH))
        
        #Counter collapses list into dictionary where mer (key) and Counter
        #appends count (val)
        region_mers = Counter(region_mers)
        
        #iterates over probes within a repeat region
        for probe in probes_list:
            
            #creates lists for each probe's k-mer counts in the list
            probe_region_count_list = []
            probe_genome_count_list = []
            
            #implement function to split probe into its k-mers of specified
            #length, returns the mers as a list
            #kmers generated checked in assertion in generate_kmers function
            probe_mers = generate_kmers(probe,int(MERLENGTH))
            
            #iterate over each mer in the probe_mers list for each probe
            for mer in probe_mers:
                
                #implement function to generate RC
                #if there are RC matches in repeat or genome, those counts
                #are also considered
                rev_mer = rev_comp(mer)
                
                #then check if mer or mer RC in region
                #append that count value
                if mer in region_mers:
                    probe_region_count_list.append(region_mers[mer])

                if rev_mer in region_mers:
                    probe_region_count_list.append(region_mers[rev_mer])
                    
                #then check if the mer or mer RC is in genome
                #append that count value
                if mer in jf_dict:
                    probe_genome_count_list.append(jf_dict[mer])

                if rev_mer in jf_dict:
                    probe_genome_count_list.append(jf_dict[rev_mer])
            
            #takes the sum of all k-mer counts identified for each probe
            #adds probe (key) and total sum of k-mer counts (val)

            all_repeat_counts[probe] = sum(probe_region_count_list)
            all_hg_counts[probe] = sum(probe_genome_count_list)

            #checks that the total k-mer sums match expected values
            assert(sum(probe_region_count_list) == 24)

            assert(sum(probe_genome_count_list) == 24)

    return probe_df,all_repeat_counts,all_hg_counts

##############################################################################

def generate_kmers(sequence,k_size):
    """
    This function takes a sequence, splits it into mers based on user k_size,
    returns mers as a list
    Parameters
    ----------
    sequence : string
        sequence to be split into mers of specified size
    k_size : int
        the size of k-mer desired
    Returns
    -------
    kmers : list
        list of k-mers depending on specified size
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

def test_generate_kmers():
    """
    This function takes a sequence, splits it into mers based on user k_size,
    returns mers as a list
    Parameters
    ----------
    sequence : string
        sequence to be split into mers of specified size
    k_size : int
        the size of k-mer desired
    Returns
    -------
    kmers : list
        list of k-mers depending on specified size

    Assertions tested
    -----------------
    - A provided sequence can be split into 18-mers
    """

    sequence = "ATCGAATCTGGTTATCGATC"
    k_size = 18
    expected_kmers = ['ATCGAATCTGGTTATCGA', 'TCGAATCTGGTTATCGAT', 'CGAATCTGGTTATCGATC']
    #make a list to store the kmers
    kmers = []
    
    #compute size of kmer window
    n_kmers = len(sequence) - int(k_size) + 1

    #generate kmer within window
    for i in range(n_kmers):
        kmer = sequence[i:i + int(k_size)]
        kmers.append(kmer)

    assert(expected_kmers == kmers)

    return kmers

##############################################################################

def test_append_probe_df(jf_file,fasta_file,probe_file):
    """
    This function takes a probe dataframe, and the dictionary for repeat
    and genome counts and appends the counts according to probes in cols.
    Parameters
    ----------
    probe_df : dataframe 
        dataframe that contains probe information and sequences
    all_repeat_counts : dictionary
        contains probe seq as key and total sum of all k-mer counts found in
        target repeat region as val
    all_hg_counts : dictionary
        contains probe seq as key and total sum of all k-mer counts found in
        human genome.
    Returns
    -------
    probe_df : dataframe
        contains probe infornation with addition of repeat count, hg count, 
        and the k_score which is the r_count/hg_count which determines how
        many of the probes k-mers exclusively target the repeat region of
        interest

    Assertions tested
    -----------------
    - See that dataframe appended is not empty with appropriate cols.

    """

    probe_df,all_repeat_counts,all_hg_counts = test_repeat_count(jf_file,fasta_file,probe_file)

    #list to use to append r_count and hg_count as columns in the probe
    #dataframe
    repeat_count_list = []
    hg_count_list = []

    #list of probes
    probe_list = probe_df['probe'].tolist()
    
    #iterate that for probe in the list, append the correct total k-mer sum
    #associated with that particular probe
    for probe in probe_list:
        repeat_count_list.append(all_repeat_counts[probe])
        hg_count_list.append(all_hg_counts[probe])

    #add the lists as columns and cast as ints
    probe_df["r_count"] = repeat_count_list
    probe_df['r_count']=probe_df['r_count'].astype(int)
    probe_df["hg_count"] = hg_count_list
    probe_df['hg_count']=probe_df['hg_count'].astype(int)

    #compute k_score, which is the k-mer binding proportion of the total
    #sum of k-mer counts within the target repeat vs the whole genome
    probe_df['k_score'] = probe_df["r_count"]/probe_df["hg_count"]
 
    #cast k_score column as float
    probe_df['k_score']=probe_df['k_score'].astype(float)

    #assert that the k_score column is as expected, making other cols valid
    assert(probe_df['k_score'].values[0] == 1.0)

    return probe_df

##############################################################################

def test_compute_normalized_binding():
    """
    Function groups probes by repeat region and implements weights to rank
    order probes. Higher c1_values indicate that probes are more likely to
    be ranked by higher repeat count, higher c2 values are to implement 
    ranking mostly drivern by k-mer binding proportion.
    
    Parameters
    ----------
    probe_df : datframe
        contains probe information in a dataframe
    Returns
    -------
    probe_df : dataframe
        contains probe information with updated ranking scheme based on 
        user input of c1 and c2 values
    Assertions tested
    -----------------
    - See that the weights are as expected.
    - Validates that probes are ranked correctly.

    """

    #make the probe dataframe to match expected output, containing two seqs
    #that have different values of k_score to test ranking system


    chrom_list = ["chrX","chrX"]
    start_list = [118,159]
    stop_list = [158,198]
    seq_list = ["AAATGATTGGTTCAAAACCAGGGTGACAAGTTGATGAAGTC",
         "ATAATAGGTTTATGATGCTGATCACCCTCTCTCACCTGGT"]
    tm_list = [42.07,42.19]
    region_list = ["chrX:1-300","chrX:1-300"]
    r_count = [24,24]
    h_count = [24,24]
    k_score = [1.0,1.0]

    #load into dataframe
    probe_df = pd.DataFrame({'chrom':chrom_list,
        'start':start_list,'stop':stop_list,
        'probe':seq_list,'tm':tm_list,'regions':region_list,
        'r_count':r_count,'h_count':h_count,'k_score':k_score})

    #list for the normalized binding value
    norm_val_list = []
    
    #dictionary to contain probe sequence as key and normalized value as val
    all_probe_norm_dict = {}

    #here in test case weights are provided
    c1 = 1
    c2 = 2

    #iterate by repeat region
    grouped_regions = probe_df.groupby(['regions'],sort = False)
    for name,group in grouped_regions:
        
        #store the value to update and add to the dictionary
        val_list = []
        
        #make a list for the probes
        probe_list = group['probe'].tolist()

        #make a list of the repeat count total and k_binding score
        repeat_counts = group["r_count"].tolist()
        k_binding = group["k_score"].tolist()

        #take the max of the group for the repeat count and k_binding score
        max_r_sum = max(repeat_counts)
        max_k_binding = max(k_binding)

        #do math to produce normalized value
        #in order of probes list
        for rc,kb in zip(repeat_counts,k_binding):
            
            #takes k-mer count in repeat and divides by max k-mer count in
            #entire repeat region and multiplies by the c1 weight.
            #This value is added to the k-mer binding score which is divided
            #by the max k-mer binding score in the repeat region and is
            #multiplied by the c2 weight. The sum of these two normalized
            #values drives ranking of each of these probes within each
            #region
            
            norm_binding_val = (((rc/max_r_sum)*c1) +
                                ((kb/max_k_binding)*c2))
            
            #append this normalized value to the list
            val_list.append(norm_binding_val)

        #zip the probe to the value as a dictionary
        probe_norm_val_dict = dict(zip(probe_list,val_list))
                
        #update this dictionary for the whole df
        all_probe_norm_dict.update(probe_norm_val_dict)

    #get a total probe list of normalized values
    probe_list = probe_df['probe'].tolist()
    for probe in probe_list:
        norm_val_list.append(all_probe_norm_dict[probe])

    #validates that the test probes provided match dict values
    check_probe_dict = {'AAATGATTGGTTCAAAACCAGGGTGACAAGTTGATGAAGTC': 1.0, 
    'ATAATAGGTTTATGATGCTGATCACCCTCTCTCACCTGGT': 1.0}

    #checks that probe was appended with correct rank value
    assert(set(all_probe_norm_dict) == set(check_probe_dict))

    #append this list to the aggregate df of probe information
    probe_df['norm_vals'] = norm_val_list

    #sort probes by the normalized values
    probe_df = probe_df.sort_values(by=['norm_vals'], ascending=False)

    #assert that the values match based on the constants provided
    check_val_list = [3.0, 3.0]

    assert(norm_val_list == check_val_list)

    return probe_df

##############################################################################

def test_write_file(o_path):
    """
    Function takes the dataframe with probe data and writes to a specified
    user output file
    
    Parameters
    ----------
    probe_data : dataframe
        dataframe containing probe information in a ranked order
    o_path : file 
        path for the dataframe to be written as an output file
        
    Returns
    -------
    None. Returns a file with the dataframe information

    Assertions tested
    -----------------
    - See that the end file is made correctly
    """

    probe_df = test_compute_normalized_binding()

    #write dataframe to file
    probe_df.to_csv(str(o_path), header= False, index=False, sep="\t")

##############################################################################

