#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Mon Jun 28 10:41:18 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################
#specific script name
script_name = "design_probes"

#import libraries
import time
import argparse
import subprocess
import numpy as np
import pandas as pd
import refactoredBlockparse as bp

#import biopython libraries
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO

##############################################################################

def test_make_fasta_from_bed(bed, region_fa, genome_fa):
    """
    This function will run a bedtools process on a genomic fasta provided,
    from a bed file to return a fasta of the bed file regions listed.

    Parameters
    ----------
    bed: bed file
    File containing the genomic coordinates to return a multi-entry fasta.
    These genomic coordinates are the identified repeat regions from the
    repeat identification script.

    region_fa: multi-entry fasta file
    The output fasta file of sequences from the generated bed file from
    the repeat identification script.

    genome_fa: fasta file name
    Reference fasta file to be used to create fasta seqs against repeats.

    Returns
    -------
    genome_fa described above

    Assertions tested
    -----------------
    - Subprocess opens fasta reference and bed file
    - Output successfully is formatted as fasta file
    - Validates file arguments passed into function
    """

    #sequence of corresponding bed coordinates chr1:2-30
    #aassertion equates if valid bases are called
    test_seq = "CCCTAAACCCTAACCCCTAANNNNNNNN"

    #takes a genomic fasta reference and bed coordinates as input
    #output is fasta file derived from bed coordinates
    subprocess.call(['bedtools', 'getfasta', '-fi', genome_fa, '-bed',
                     bed, '-fo', region_fa], stderr=None, shell=False)

    #opens the fasta file
    fasta_sequences = SeqIO.parse(open(region_fa),'fasta')
    
    #parses the chromosome of interest from file
    for fasta in fasta_sequences:
        sequence = (fasta.seq)

    assert(sequence == test_seq)

##############################################################################

def test_blockParse_run(test_regions,name,probe_out):
    """
    This function takes the provided fasta seqs derived from the bed file
    and passes each repeat seq into the refactoredBlockParse script from
    Oligominer. Here, probes are written and appended to an output dataframe.
    runSequenceCrawler is run using parameters described in Oligominer for 
    default probe generation. 

    Parameters
    ----------
    region_fa: fasta file
    Contains fasta sequences from derived repeat regions

    name: string
    The chromosome that repeat regions are derived from

    Returns
    -------
    probe_out: tsv file
    The output file name containing all designed probes for all provided
    repeat region fasta sequences.

    Assertions tested
    -----------------
    - Tests if fasta names are appended correctly to dataframes
    - Using secondary test file provided probes are designed.
    """
    name_list = []
    sequence_list = []

    fasta_sequences = list(SeqIO.parse(open(test_regions),'fasta'))

    #you want to parse each fasta sequence as a string and append
    #to a sequence list
    
    for fasta in fasta_sequences:
        name_list.append(fasta.id)
        sequence=(str(fasta.seq))
        sequence_list.append(sequence)

    #zip the names (headers of the fasta) and the string sequence into a list,
    #then make into a dict
    zipped_list=zip(name_list,sequence_list)
    dict_name_seq=dict(zipped_list)
 
    #then run each item in the dict into the refactored blockParse script
    #which can now handle multi-lined fastas
    for names,sequences in dict_name_seq.items():
        
        probes=bp.runSequenceCrawler(sequences,names,name, 36, 41, 20, 80, 
                                     mt.DNA_NN3, 42, 47,
                                     'AAAAA,TTTTT,CCCCC,GGGGG', 390, 50, 0,
                                     25, 25, None, True ,False, False, False,
                                     False, False,(str(probe_out)))

    #reads in probes that have been designed
    colnames = ["chrom","p_start","p_stop","probe","Tm", "regions"]

    probe_df=pd.read_csv(probe_out, delimiter = '\t',names=colnames)

    repeat_list = probe_df['regions'].tolist()

    #collapse uniques
    repeat_list = (set(repeat_list))

    #checks if all repeats provided in file have probes designed 
    assert(repeat_list == set(name_list))

##############################################################################

