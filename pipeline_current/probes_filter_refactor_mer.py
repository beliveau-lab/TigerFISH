"""
Robin Aguilar
Created:05.09.2021
Purpose:refactor this script so it filters more efficiently by RR
"""

import argparse
import time
import pandas as pd
import nupack
from nupack import *

from Bio.Seq import reverse_complement as rev_comp

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

MER_CUTOFF = 0.95

start_time=time.time()

###################################################################################

def read_region(file_path):
    
    #read in the file with the appropriate column names
    colnames = ["chrom","p_start","p_end","probe","Tm","region","r_count","h_count","max_mer","k_score"]
    
    #read as a dataframe
    region_df = pd.read_csv(file_path, delimiter = '\t', names = colnames)
        
    region_df['k_score']=region_df['k_score'].astype(float)

    region_df=region_df[(region_df['k_score'] >= float(ENRICH)) & (region_df['r_count'] >= int(COPY_NUM))]

    region_df = region_df.drop_duplicates(subset=['probe'], keep='first')
    
    region_df=region_df.sort_values(by='r_count', ascending=False)

    return region_df

def rm_shared_mer_probes(region_df):
            
    probe_to_remove = []
    
    #do a groupby over the repeat regions
    grouped_regions=region_df.groupby('region',sort=False)
    
    for name,group in grouped_regions:
        if len(group) > 1:
            
            #gets a list of all max_mers
            max_mers_list = group['max_mer'].tolist()
            max_mers_list = [i[1 : -1].split(', ') for i in max_mers_list]
    
            #define your top ranking "keep set of 18-mers"
            keep_mers_set = set(max_mers_list[0])
        
            #get a list of all probes
            probes_list = group['probe'].tolist()
            
            for i in range(1,len(max_mers_list)):
                #will keep track of number of 18-mers within a given mers list in the keep set
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
                if float(in_count/list_len) >= MER_CUTOFF:
                    probe_to_remove.append(probes_list[i])
    
    #remove rows of the dataframe that have probes that should be filtered
    region_df = region_df[~region_df['probe'].isin(probe_to_remove)]
            
    return region_df
                    
def write_file(region_df):

    region_df.to_csv(str(out_path), header=False, index=False, sep="\t")

def main():
    
    region_df = read_region(file_path)
    print("---%s seconds ---"%(time.time()-start_time))

    region_df = rm_shared_mer_probes(region_df)
    print("---%s seconds ---"%(time.time()-start_time))
        
    write_file(region_df)
    print("---%s seconds ---"%(time.time()-start_time))
    
    print("done")
    
if __name__== "__main__":
    main()

