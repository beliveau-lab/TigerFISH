#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Thu Jul  1 13:13:48 2021
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
import subprocess
from subprocess import Popen

##############################################################################

def generate_bins(file_path,out_path,genome_name):
    
    genome_bin_file = str(out_path) + str(genome_name) + ".txt"
    
    #use subprocess to call bowtie
    subprocess.Popen(['bedtools', 'makewindows', '-g', str(file_path),
                     '-w', '1000000', genome_bin_file], \
    stdout=subprocess.PIPE)
    
##############################################################################

def main():
    
    start_time=time.time()

    userInput = argparse.ArgumentParser(description=\
        '%Requires a chrom sizes file to generate a genome bin based'
        'on the size of the desired window')
        
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file_path', action='store', 
                               required=True, help='The genome file'
                               'containing chrom sizes')
    requiredNamed.add_argument('-o', '--out_path', action='store',
                               required=True, help='The directory where'
                               'the bin file should be kept')
    requiredNamed.add_argument('-g', '--genome_name', action='store',
                               required=True, help='The name of the genome')
    
    args = userInput.parse_args()
    file_path = args.file_path
    out_path = args.out_path
    genome_name = args.genome_name
    
    
    generate_bins(file_path,out_path,genome_name)
    
    print("---%s seconds ---"%(time.time()-start_time))
    
    print("done")

if __name__== "__main__":
    main()
