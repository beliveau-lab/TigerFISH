"""
Robin Aguilar
Beliveau and Noble Labs, 08/2/2019
Dependencies: Jellyfish, Bedtools
Input: Fasta for region of interest, JF index of genome of interest
Output: A .bed file of the elevated kmer regions, a .fa of said regions,
and a DF containing information about all the probes generated in those
regions. 
"""

#import all functions/modules needed
import time
import argparse
import subprocess
import numpy as np
import pandas as pd
import refactoredBlockparse as bp
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
import tempfile

##############################################################################

#declare a timer
start_time=time.time()

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-b', '--bed_name', action='store', required=True,
                               help='The chromosome corresponding to the probes being generated')
requiredNamed.add_argument('-r_o', '--region_out', action='store', required=True,
                               help='The chromosome corresponding to the probes being generated')
requiredNamed.add_argument('-p_o', '--probes_out', action='store', required=True,
                               help='The chromosome corresponding to the probes being generated')
requiredNamed.add_argument('-g', '--genome_fasta', action='store', required=True,
                               help='The chromosome corresponding to the probes being generated')
requiredNamed.add_argument('-c', '--chrom_name', action='store', required=True,
                               help='The chromosome corresponding to the probes being generated')
args = userInput.parse_args()
bed = args.bed_name
region_fa = args.region_out
probe_out = args.probes_out
genome_fa = args.genome_fasta
name = args.chrom_name

##############################################################################

def make_fasta_from_bed(bed,region_fa):

    subprocess.call(['bedtools', 'getfasta', '-fi', genome_fa, '-bed', bed, '-fo', region_fa], stderr=None, shell=False)

##############################################################################

#this function takes the fasta that you returned from the last function
def blockParse_run(region_fa,name):

    name_list = []
    sequence_list = []

    fasta_sequences = list(SeqIO.parse(open(region_fa),'fasta'))
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
        probes=bp.runSequenceCrawler(sequences,names,name, 36, 41, 20, 80,mt.DNA_NN3,42, 47, 'AAAAA,TTTTT,CCCCC,GGGGG', 390, 50, 0,25, 25, None, True ,False, False,False, False, False,(str(probe_out)))

##############################################################################

def main():

    make_fasta_from_bed(bed,region_fa)

    print("---%s seconds ---"%(time.time()-start_time))

    blockParse_run(region_fa,name)
    
    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
