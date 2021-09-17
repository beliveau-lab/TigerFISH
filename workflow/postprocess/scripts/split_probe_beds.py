#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Thu Jul  1 10:43:36 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#specific script name
script_name = "split_probe_beds"

#load libraries used
import time
import argparse
import pandas as pd
from itertools import groupby

##############################################################################

def read_probe_bed(file_path):
    """
    Function takes the bed file provided and reads it into a dataframe

    Parameters
    ----------
    file_path : file
        bed file with probe coordinates.

    Returns
    -------
    probe_df : dataframe
        contains probe information in dataframe

    """

    #names of columns in file
    colnames = ["probe_coords","repeat_coords","probe"]
    
    #read in as dataframe
    probe_df = pd.read_csv(file_path, delimiter = '\t', names = colnames)

    return probe_df
    
##############################################################################

def split_file(probe_df,out_path):
    """
    Function takes the probe dataframe and splits each probe coord into a 
    distince bed file for downstream pipeline parsing.

    Parameters
    ----------
    probe_df : dataframe
        pandas dataframe containing probe information
    out_path : path to output file
        contains the path to an output probe file

    Returns
    -------
    None. creates output file

    """
    
    #do a groupby over the distict coord start values
    grouped_regions=probe_df.groupby("probe_coords",sort=False)

    #write a probe coord to a seperate file
    for name,group in grouped_regions:
        #use dict to call on full coord of probe 
        group.to_csv(str(out_path) + str(name) + ".bed",
                     header=False,index=False, sep="\t")
        
##############################################################################

def main():
    
    start_time=time.time()

    userInput = argparse.ArgumentParser(description=\
        '%Requires a probe bed file to be split into seperate files'
        'based on probe coordinates')

    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file_path', action='store', 
                               required=True, help='The post probe file'
                               'that contains all sequences and regions')
    requiredNamed.add_argument('-o', '--out_path', action='store',
                               required=True, help='The directory where'
                               'the split files will be located by repeat')
    
    args = userInput.parse_args()
    file_path = args.file_path
    out_path = args.out_path
    
    probe_df = read_probe_bed(file_path)

    print("---%s seconds ---"%(time.time()-start_time))
    
    split_file(probe_df,out_path)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    print("done")

if __name__== "__main__":
    main()
