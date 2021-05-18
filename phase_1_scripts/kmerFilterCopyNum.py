#########################
#!/usr/bin/env kmerFilter
###5/22/2019
###Robin Aguilar
###Beliveau lab
###Input:The output of the repeatFinder.py script
###Dependencies:kmerFilter.py from Oligominer
###Output:A dataframe containing the following information about the kmers for each of the probes designed:
###'chr','start','end','start_base','end_base','probe_seq','temp','m_length','monomer','copy_num','r_length','probes','ciel_copy_num','new_k','found','total','Kmer_Enrich'
###Note that the output of this should be depdendent on the K_adjust value.
###Purpose: For each probe designed, it will give the number of kmers depending on the argument you give for k_adjust. K_adjust is what is used in the kmerFilter.py algorithm
###The output of this script is to be used in the kmer_max_mer.py script 
########################

import time
import io
import sys
import re
import csv
from collections import defaultdict
from itertools import islice
from operator import itemgetter, attrgetter
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
import numpy as np
import pandas as pd
import glob
import os
import kmerFilter
import timeit
import itertools


rootdir = '.'
extensions = glob.glob(rootdir + "/*.bed")
K_ADJUST=2.0

file_list=[]
text_file=glob.glob('*.bed')
for i in text_file:
    file_list.append(i)
#print(file_list)

TRF_output="2019_04_26_candidate_probes_list.txt"
#open the output text file and define the start stop, etc values
start_list=[]
stop_list=[]
with open(TRF_output,'r') as file:
    for i, line in enumerate(file):
        if i >= 1:
            #format the candidate probe file contents
            replace_space=line.replace(' ','\t')
            elements=replace_space.split()
            chr=elements[0]
            start=elements[1]
            start_list.append(start)
            #make a list and store all the start positions
            stop=elements[2]
            stop_list.append(stop)
            #make a list and store all the stop positions
            m_length=elements[3]
            monomer=elements[4]
            copynum=elements[5]
            r_length=elements[6]
            probes=elements[7]

#this will open all files in the current directory that end in *bed and print their contents
bed_file_contents=[]
concat_probe_file=[]
for file in file_list:
    with open(file, 'r') as open_file:
    #Since the file name is chr_start_stop.bed
    #element 0 should be chr
        names = file.split('_')
        chrom_names=names[0]
        #element one will print the start portion of the title
        start_seq=names[1]
        #print(start_seq)
        #element 2 will print the end position of the title
        end_seq=file.replace('_',' ').replace('.',' ').split()[2]

        #iterates through the list containing all the start sites
        for i in start_list:
            #iterates through the list containing all the top sites
            for j in stop_list:
                #if the start in the title matches the value of the start in the start_list, does the same for the end
                if(int(i)==int(start_seq) and int(j)==int(end_seq)):
                    for l, line in enumerate(open_file):
                        #replaces the spaces into tabs and reports the details for the scaffold, start, end, probe, and temp for each probe in the loop
                        replace_space = line.replace(' ', '\t')
                        bed_contents = replace_space.split()
                        scaffold=bed_contents[0]
                        start_base=bed_contents[1]
                        end_base=bed_contents[2]
                        probe_seq=bed_contents[3]
                        temp=bed_contents[4]
                        #writes all of these details to a large list that contains the TRF output information and then the concatenated probe file information
                        #the concatenated probe file contains information that a .bed output would contain from blockParse, but with the probes from all bedfiles into one file
                        bed_file_contents.append(str(chrom_names)+ '\t' + str(i) + '\t' + str(j) + '\t'+ str(scaffold) +
                                    '\t' + str(start_base) + '\t' +str(end_base) +
                                    '\t' + str(probe_seq) + '\t' + str(temp) + '\n')
                        concat_probe_file.append(str(scaffold) + '\t' + str(start_base) + '\t' +
                                                    str(end_base) + '\t' + str(probe_seq) + '\t' + str(temp) + '\n')
                            #write the total candidate probe information output with the TRF information
                    with open('candidate_probe_kmerFilter_output.txt', 'w') as iter:
                        for b in bed_file_contents:
                            iter.write("%s" % b)
                    iter.close()
                    #write the concatenated probes into one big text file
                    with open('concatenated_probes.txt','w') as probes:
                        for p in concat_probe_file:
                            probes.write("%s" % p)
                    probes.close()

#contains all probe information
kmer_candidates="candidate_probe_kmerFilter_output.txt"
#contains all concatenated probes
candidate_probes="concatenated_probes.txt"
probe_list=[]

#open the kmer_candidates file
with open(kmer_candidates,'r') as kfile:
    #open the TRF output file
    with open(TRF_output, 'r') as tfile:
        file1 = pd.read_csv(kmer_candidates, delimiter="\t")
        #make file1 into pandas dataframe containing all probe related relevant information
        file1.columns = ['chr', 'start', 'end', 'scaffold', 'start_base', 'end_base', 'probe_seq', 'temp']
        file2=pd.read_csv(TRF_output, delimiter='\t')
        #make file2 into pandas dataframe containing all TRF related information for the probes
        file2.columns = ['chr', 'start', 'end', 'm_length', 'monomer', 'copy_num', 'r_length', 'probes']

        #merge the two dataframes on the start values since we know that this is unique
        file3=pd.merge(file1,file2,on='start')

        #remove the excessive columns containing the chr and end values
        file3 = file3.drop(["chr_y", "end_y"], axis=1)

        #round up the copy number by making a new column that goes to cieling
        file3['ciel_copy_num']=file3['copy_num'].apply(np.ceil)

        #multiplies this copy number value by the K_ADJUST argument specified above
        #creates a new column called new_K, this value is to be used as the K argument for kmerFilter
        file3['new_k']=(file3['ciel_copy_num']*K_ADJUST).apply(np.ceil)

        #write this intermediate text file as a test that contains all this information above from the data frame
        file3.to_csv('intermediate_test.txt', header=True, index=None, sep='\t')

        #with the candidate probe file above, you will want to find a way to run this with kmerFilter.py
        with open(candidate_probes, 'r') as pfile:
            for i, line in enumerate(pfile):
                replace_space = line.replace(' ', '\t')
                elements = replace_space.split()
                scaffold = elements[0]
                start_base = elements[1]
                stop_base = elements[2]
                probe=elements[3]

                #this is a list that appends all the probe sequences
                probe_list.append(probe)

        k_val_list=[]
        all_probe_information="intermediate_test.txt"
        probe_info = pd.read_csv(all_probe_information, delimiter="\t")
        probe_info.columns=['chr','start','end',
                            'scaffold','start_base',
                            'end_base','probe_seq','temp',
                            'm_length','monomer','copy_num','r_length',
                            'probes','ciel_copy_num','new_k']

        probe_info['Grouper'] = probe_info['start'] + probe_info['end'] + probe_info['new_k']
        probe_info = probe_info.drop_duplicates('Grouper')

        subset_probe_info = probe_info[['chr', 'start','end', 'new_k']]
        #print(subset_probe_info)
        subset_probe_info.to_csv('probe_coordinates_kval.txt', header=False, index=None, sep='\t')

        probe_coord_k="probe_coordinates_kval.txt"
        with open(probe_coord_k,'r') as p_kval:
            for i, line in enumerate(p_kval):
                columns=line.split()
                for file in file_list:
                    file_name=file.replace("_",'\t').replace(".",'\t')
                    split_file_name=file_name.split()

                    if(columns[0]==split_file_name[0] and columns[1]==split_file_name[1]):
                        k_val_list.append(file + '\t' + columns[3] + '\n')

                    with open('file_name_k_val.txt','w') as file_info:
                        for f in k_val_list:
                            file_info.write("%s" % f)
                    file_info.close()

        file_k="file_name_k_val.txt"
        with open(file_k, 'r') as f_k:
            for i, line in enumerate(f_k):
                instance=line.split()
                kmerFilter.runFilter(instance[0], instance[0], 18, 'hg38_alt_patch.jf', float(instance[1]), None, False, False, True, 0)

        meta_file_list = []
        chr_split=[]
        start_split=[]
        stop_split=[]
        meta = glob.glob('*_meta.txt')
        for i in meta:
            meta_file_list.append(i)
            i = i.replace("_", '\t').replace(".", '\t')
            file_name_split = i.split()
            chr_split.append(file_name_split[0])
            start_split.append(file_name_split[1])
            stop_split.append(file_name_split[2])

        to_merge_meta_contents=[]
        store_file_names=[]


        for file in meta_file_list:
            with open(file, 'r') as o_file:
                for i, line in enumerate(o_file):
                    split_meta_column = line.split()
                    store_file_names.append(split_meta_column[0])
                    split_on=split_meta_column[0].replace("_",'\t').replace(".","\t")
                    file_split=split_on.split()
                    to_merge_meta_contents.append(file_split[0] + '\t' +
                          file_split[1] + '\t' +
                          file_split[2] + '\t' +
                          split_meta_column[5] + '\t' +
                          split_meta_column[6] + '\n')

            with open('meta_contents_merged.txt', 'w') as meta:
                for m in to_merge_meta_contents:
                    meta.write("%s" % m)
                meta.close()

        meta_contents="meta_contents_merged.txt"
        with open(meta_contents,'r') as m:
            meta_info = pd.read_csv(m, delimiter="\t")
            meta_info.columns=['chr','start','end','found','total']

            collapse_dataframes=(pd.merge(file3, meta_info, on='start'))

            collapse_dataframes['Kmer_Enrich'] = (collapse_dataframes['found']/collapse_dataframes['total'])*100

            collapse_dataframes = collapse_dataframes.drop(["scaffold","chr", "end"], axis=1)

            collapse_dataframes.columns=['chr','start','end','start_base',
                                         'end_base','probe_seq','temp','m_length',
                                         'monomer','copy_num','r_length','probes',
                                         'ciel_copy_num','new_k','found','total',
                                         'Kmer_Enrich']

            collapse_dataframes.to_csv('final_meta_contents.txt', header=True, index=None, sep='\t')
