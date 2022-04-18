#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##############################################################################

"""
Created on Fri Jun 25 09:54:15 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""

##############################################################################

#specific script name
scriptName = "split_bed"

#import libraries
import time
import pandas as pd
import argparse

##############################################################################

def test_read_bed(bed_file):
    
    """
    Function takes the bed file and reads it into a pandas df.
    Parameters
    ----------
    bed_file : bed file
        file with region coords chr \t start \t stop
    Returns
    -------
    bed_df : dataframe
        dataframe of coordinates in file

    Assertions tested
    -----------------
    - Can open and read bed file into a dataframe.

    """
    #make a dataframe that matches the provide test  file
    chrom_list = ["chrX","chr2"]
    start_list = [1,1]
    stop_list = [50,100]

    test_df = pd.DataFrame(
        {'chrom': chrom_list,
         'start': start_list,
         'stop': stop_list
        })
    
    colnames = ["chrom","start","stop"]
    
    bed_df = pd.read_csv(bed_file, names=colnames,header=None,
                            sep='\t')

    #check if tested file matches expected result of test_df    
    assert(bed_df.equals(test_df))

    return bed_df
    
##############################################################################

def test_split_write_bed(bed_file,chrom,bed_out):
    """
    Function takes dataframe and subsets specific chrom into an output bed
    file with regions only corresponding to the chrom of choice.
    Parameters
    ----------
    bed_df : dataframe
        dataframe of region coordinates from file
    Returns
    -------
    None. Writes file with chromosome name as bed.
    """    

    bed_df = test_read_bed(bed_file)

    #generate new test dataframe to probe correct chrom is subset
    chrom_list = ["chrX"]
    start_list = [1]
    stop_list = [50]

    test_df = pd.DataFrame(
        {'chrom': chrom_list,
         'start': start_list,
         'stop': stop_list
        })

    #subsets rows containing specific chromosome
    chrom_df = bed_df[bed_df['chrom'].str.contains(str(chrom))]

    #checks that expected df matches output    
    assert(chrom_df.equals(test_df))

    #writes row(s) to output bed file
    chrom_df.to_csv(bed_out,index = False,
                               index_label=False,
                               header = None,
                               sep = '\t')
    
    
##############################################################################

