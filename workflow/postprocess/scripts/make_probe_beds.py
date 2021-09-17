#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Wed Jun 30 18:17:36 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#import libraries

import argparse
import time
import pandas as pd

##############################################################################

def read_probe_file(p_file):
    """
    Function reads in provided probe file following analysis via Tigerfish
    pipeline.
    
    Parameters
    ----------
    p_file : file
        probe file after Tigerfish pipeline is completed, contains subset
        of probes of interest to survey in silico binding predictions

    Returns
    -------
    probe_df : dataframe
        probe information is stored in a pandas dataframe

    """
    
    #add columns to dataframe
    colnames = ['probe_coords','region','probe','Tm','r_count',
                'h_count','k_score','k_norm','on_target_binding',
                'off_target_binding','on_target_binding_prop']
    
    #load in file as dataframe
    probe_df = pd.read_csv(p_file, delimiter = '\t', names = colnames)

    return probe_df
    
##############################################################################

def write_clean_probe_df(probe_df,or_file):
    """
    Function will remove columns not needed for probe to repeat mapping
    
    Parameters
    ----------
    probe_df : dataframe
        dataframe storing probe information
    or_file : file
        file containing probe coords, repeat region, and probe sequence

    Returns
    -------
    probe_df : dataframe
        dataframe containing probe coord, repeat region, and probe sequence
        columns

    """
    
    #drop columns so only the probe coords and region remain so that way
    #it's possible to map back which alignment sequences fall within the 
    #probe's target region
    
    probe_df = probe_df.drop(['Tm', 'r_count', 'h_count', 'k_score',
                              'k_norm', 'on_target_binding',
                              'off_target_binding', 'on_target_binding_prop'],
                             axis = 1)
    
    #write as output file
    probe_df.to_csv(or_file, header=False, index=False, sep="\t")

    return probe_df
    
##############################################################################

def write_probe_bed(probe_df,ob_file):
    """
    Function containing probe dataframe and output bed file

    Parameters
    ----------
    probe_df : dataframe
        contains remaining probe information
    ob_file : file
        bed file containing coordinates of probes

    Returns
    -------
    None. Creates file with bed coordinates of selected probes

    """
    
    #split into two columns for chrom and coords
    probe_df[['chrom','coords']] = probe_df['probe_coords'].str.split(':',
                                                                      expand=True)

    #split coord column into start and stop vals
    probe_df[['p_start','p_end']] = probe_df['coords'].str.split('-',
                expand=True)

    #cast coord columns as ints
    probe_df[['p_start','p_end']] = probe_df[['p_start',
                'p_end']].astype(int)
    
    
    #remove additional columns
    columns = ["region","probe","probe_coords","coords"]
    probe_df = probe_df.drop(columns,axis=1)

    #rearrange cols
    probe_df = probe_df[['chrom','p_start', 'p_end']]
    
    #write as output file
    probe_df.to_csv(ob_file, header=False, index=False, sep="\t")
    
##############################################################################

def main():

    start_time=time.time()

    userInput = argparse.ArgumentParser(description=\
        'Requires subset of probes of interest derived after alignment'
        'filter has been implemented in pipeline. This script will'
        'clean probe coords into bed files and store a simplified'
        'dataframe of probe coords and corresponding repeat coords')
    
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--probe_file', action='store',
                               required=True, help='The filtered probe file'
                               'that contains all sequences and regions')
    requiredNamed.add_argument('-o_b', '--out_bed', action='store',
                               required=True, help='The filtered probe file'
                               'that contains all sequences and regions')
    requiredNamed.add_argument('-o_r', '--out_regions', action='store',
                               required=True, help='The filtered probe file'
                               'that contains all sequences and regions')
    
    args = userInput.parse_args()
    p_file = args.probe_file
    ob_file = args.out_bed
    or_file = args.out_regions


    probe_df = read_probe_file(p_file)
    print("---%s seconds ---"%(time.time()-start_time))

    probe_df = write_clean_probe_df(probe_df,or_file)
    print("---%s seconds ---"%(time.time()-start_time))

    write_probe_bed(probe_df,ob_file)
    print("---%s seconds ---"%(time.time()-start_time))

    print("done")    
    
if __name__== "__main__":
    main()
