"""
Robin Aguilar
Beliveau and Noble Labs, 10/14/2019
Dependencies: Jellyfish, Bedtools
Input: Fasta of k-mer enriched regions. 
Output: A DF containing probes generated in each region.
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
import refactoredBlockparse as bp
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
from scipy import signal
import time
from collections import OrderedDict
import collections
from operator import itemgetter
from itertools import groupby
import collections
from collections import Counter
import itertools

#declare a timer
start_time=time.time()

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--fasta_file', action='store', required=True,
                               help='The FASTA file to find probes in')

args = userInput.parse_args()
test_fasta = args.fasta_file

name_list=[]
sequence_list=[]

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
        probes=(bp.runSequenceCrawler(sequences,names, 36, 41, 20, 80,mt.DNA_NN3,
                                      42, 47, 'AAAAA,TTTTT,CCCCC,GGGGG', 390, 50, 0,25,
                                      25, None, True ,False, False,True, False, False))

    return probes

def main():

    blockParse_run(test_fasta)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()

