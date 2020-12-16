"""
Robin Aguilar
Beliveau and Noble labs
Date Modified: 03/01/2020
Modification: To make some functions more efficient and commented with user 
statements 
"""
#import libraries needed
from os import listdir
from os.path import isfile, join
import time
import argparse
import pandas as pd
###################################################################################

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--file_path', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')
requiredNamed.add_argument('-o', '--out_path', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')
userInput.add_argument('-e', '--enrich_score', action='store', default=0.50,required = True, 
                           type=float,
                           help='The max # of times a given k-mer is found in a repeat region/entire HG, default is '
                                '0.50')
userInput.add_argument('-cn', '--copy_num', action='store', default=10,required = True,
                           type=int,
                           help='The minimum allowed number of times a k-mer is identified to be added to the final probe list, default is '
                                '10')

args = userInput.parse_args()
file_path = args.file_path
out_path = args.out_path
ENRICH=args.enrich_score
COPY_NUM=args.copy_num

start_time=time.time()

###################################################################################

def filter_file(file_path):

    colnames = ["chrom","p_start","p_end","probe","Tm","region","r_count","h_count","k_score"]

    probe_data = pd.read_csv(file_path, delimiter = '\t', names = colnames)

    probe_data['k_score']=probe_data['k_score'].astype(float)

    probe_data=probe_data[(probe_data['k_score'] >= float(ENRICH)) & (probe_data['r_count'] >= int(COPY_NUM))]

    probe_data = probe_data.drop_duplicates(subset=['probe'], keep='first')

    probe_data.to_csv(str(out_path), header=False, index=False, sep="\t")

###################################################################################

def main():

    filter_file(file_path)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()

