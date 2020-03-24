###################
###5/15/2019
###Robin Aguilar
###Beliveau lab
###Input: the value of the kmerFilterCopyNum.py script
###Output: a new dataframe with the maximum contents of the kmer count in the repeat v the entire genome
###Dependencies: requires Jellyfish
###Purpose: Used to conduct kmer count analysis on different probes and provide a dataframe to be used for further specificity analysis 
##################

import time
import io
import sys
import re
import csv
from collections import defaultdict
from itertools import islice
from operator import itemgetter, attrgetter
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
import subprocess
import numpy as np
import pandas as pd
import glob
import os
import timeit
import itertools


#open the file that contains all probe information

all_candidate_probe_info="final_meta_contents_1.5.txt"
chr_l=[]
start_l=[]
end_l=[]
probe_l=[]
start_base_l=[]
fasta_list=[]
with open(all_candidate_probe_info,'r') as c_probe:
    for i, line in enumerate(c_probe):
        if i > 0:
            split_line=line.split()
            chr=split_line[0]
            start=split_line[1]
            end=split_line[2]
            start_base=split_line[3]
            probe=split_line[5]
            fasta_list.append(">"+chr + "_" + start + "_" + end + "_" + probe + '\n' + probe + '\n')
    with open('large_probe_fasta.txt', 'w') as iter:
        for f in fasta_list:
            iter.write("%s" % f)
    iter.close()

    print("Making large fasta file ...")

#this will write each of the individual probes into it's respective .fa file, note that this should be deleted in the end!
    probe_fasta='large_probe_fasta.txt'
    f=SeqIO.parse(probe_fasta,"fasta")
    for rec in f:
        of=open("%s_.fa" % (rec.id), "w")
        SeqIO.write(rec,of,"fasta")
        of.close()
    f.close()

    print("Splitting into temp.fa ... ")

jellyfish_files = glob.glob('*.jf')
list_jf_file=[]
for i in jellyfish_files:
    list_jf_file.append(i)
#print(list_jf_file)

probe_fastas = glob.glob('*_.fa')
probe_fasta_l=[]
for i in probe_fastas:
    probe_fasta_l.append(i)

merge_jf_probe=[]

for i in list_jf_file:
    replace_jf = i.replace("%3A", '\t').replace(":",'\t').replace("-", '\t').replace(".",'\t')
    jf_name_split = replace_jf.split()

    for j in probe_fasta_l:
        replace_probe=j.replace("_", '\t').replace(".",'\t')
        probe_name_split=replace_probe.split()

        if(jf_name_split[1]==probe_name_split[1]):
            merge_jf_probe.append(j + '\t' + i + '\n')

    with open('probes_and_jf_map.txt', 'w') as p_jf:
        for n in merge_jf_probe:
            p_jf.write("%s" % n)
    p_jf.close()

    print("Mapping temp.fa probes to jellyfish file ... ")

    # Call jellyfish and count kmers in the probe sequences.
    probe_map = "probes_and_jf_map.txt"
    with open(probe_map, 'r') as p_m:
        for i, line in enumerate(p_m):
            instance = line.split()
            subprocess.call(['jellyfish', 'query', instance[1], '-s',
                            instance[0], '-o', instance[0]+'_temp.txt'],
                            stderr=None, shell=False)

temp = glob.glob('*_temp.txt')

list_max_coord=[]
for item in temp:
    with open(item, 'r') as jfData:
        jf_read = [line.strip() for line in jfData]

        max_list=[]
        for i in range(0, len(jf_read), 1):
            if ([int(jf_read[i].split(' ')[1])] ==0):
                pass
            else:
                max_list.append([int(jf_read[i].split(' ')[1])])
    list_max_coord.append(str(item) + '\t' + str(max(max_list)) + '\n')
    jfData.close()

for file in temp:
    os.remove(file)

with open('probe_and_max.txt', 'w') as max:
    for m in list_max_coord:
        max.write("%s" % m)
max.close()

print("Finding the max value of mer in each probe ... ")

probe_and_max="probe_and_max.txt"

probe_max_merge_l=[]
with open(probe_and_max,'r') as probe_max:
    for i,line in enumerate(probe_max):
        replace_space = line.replace(' ', '\t')
        file_name = replace_space.replace("_", '\t').replace(".", '\t').replace("[","").replace("]","")
        elements=file_name.split()
        probe_max_merge_l.append(elements[0] + '\t' + elements[1] + '\t' +
                                 elements[2] + '\t' + elements[3] + '\t' +
                                 elements[7] + '\n')

    with open('probe_max_merge.txt', 'w') as probe_merge:
        for merge in probe_max_merge_l:
            probe_merge.write("%s" % merge)
    probe_merge.close()

print("Merging probe seq with concatenated probe information ... ")

to_merge_probe_max = "probe_max_merge.txt"
to_merge_all_probe_info="final_meta_contents_1.5.txt"

with open(to_merge_probe_max, 'r') as probe_merge:
    probe_merge_info = pd.read_csv(probe_merge, delimiter="\t")
    probe_merge_info.columns=["chr","start","end","probe_seq","max_count_repeat"]

    with open(to_merge_all_probe_info,'r') as probe_all:
        probe_all_info=pd.read_csv(probe_all, delimiter="\t")

        collapse_dataframes = (pd.merge(probe_all_info, probe_merge_info, on=['chr','start','end','probe_seq']))

        #collapse_dataframes = collapse_dataframes.drop(["chr_y", "start_y", "end_y"], axis=1)

        collapse_dataframes.to_csv('meta_contents_max_merge.txt', header=True, index=None, sep='\t')

        collapse_dataframes.columns = ['chr', 'start', 'end','start_base',
                                       'end_base', 'probe_seq', 'temp', 'm_length',
                                       'monomer', 'copy_num', 'r_length', 'probes',
                                       'ciel_copy_num', 'new_k', 'found', 'total',
                                       'Kmer_enrich','max_count_repeat']

    print("Merge with into total probe list and making column ... ")
#now for the human genome iteration

hg_patch="hg38_alt_patch_36.jf"
hg_patch_list=[]
for i in probe_fasta_l:
    hg_patch_list.append(i + '\t' + hg_patch + '\n')

    with open('probe_fas_hg_jf.txt', 'w') as hg_jf:
        for merge in hg_patch_list:
            hg_jf.write("%s" % merge)
    hg_jf.close()

    print("Mapping probe to jellyfish file ... ")

probe_hg_jf="probe_fas_hg_jf.txt"
with open(probe_hg_jf, 'r') as p_hg:
    for i,line in enumerate(p_hg):
        instance = line.split()
        print(instance)

        subprocess.call(['jellyfish', 'query', instance[1], '-s',
                    instance[0], '-o', instance[0] + '_temp.txt'],
                    stderr=None, shell=False)

temp_f = glob.glob('*_temp.txt')
hg_list_max_coord = []

for item in temp_f:
    with open(item, 'r') as jfData_h:
        jf_read = [line.strip() for line in jfData_h]

        hg_max_list = []
        for i in range(0, len(jf_read), 1):
            hg_max_list.append([int(jf_read[i].split(' ')[1])])

    hg_list_max_coord.append(str(item) + '\t' + str(max(hg_max_list)) + '\n')
    jfData_h.close()

for file in temp_f:
    os.remove(file)

with open('hg_probe_and_max.txt', 'w') as max_l:
    for m in hg_list_max_coord:
        max_l.write("%s" % m)
max_l.close()

print("Finding max for each probe sequence ... ")

probe_and_max = "hg_probe_and_max.txt"

probe_max_merge_l = []
with open(probe_and_max, 'r') as probe_max:
    for i, line in enumerate(probe_max):
        replace_space = line.replace(' ', '\t')
        file_name = replace_space.replace("_", '\t').replace(".", '\t').replace("[", "").replace("]", "")
        elements = file_name.split()
        probe_max_merge_l.append(elements[0] + '\t' + elements[1] + '\t' +
                                    elements[2] + '\t' + elements[3] + '\t' +
                                    elements[7] + '\n')

    with open('hg_probe_max_merge.txt', 'w') as probe_merge:
        for merge in probe_max_merge_l:
            probe_merge.write("%s" % merge)
    probe_merge.close()

    print("Merging max value with the rest of probe information ... ")

to_merge_probe_max = "hg_probe_max_merge.txt"
to_merge_all_probe_info = "meta_contents_max_merge.txt"

with open(to_merge_probe_max, 'r') as probe_merge:
    probe_merge_info = pd.read_csv(probe_merge, delimiter="\t")
    probe_merge_info.columns = ["chr", "start", "end", "probe_seq", "max_count_hg"]

    with open(to_merge_all_probe_info, 'r') as probe_all:
        probe_all_info = pd.read_csv(probe_all, delimiter="\t")
        collapse_dataframes = (pd.merge(probe_all_info, probe_merge_info, on=['chr','start','end','probe_seq']))
        print(collapse_dataframes)
        #collapse_dataframes = collapse_dataframes.drop(["chr_y", "start_y", "end_y"], axis=1)

        collapse_dataframes.to_csv('meta_contents_max_merge_all.txt', header=True, index=None, sep='\t')

        collapse_dataframes.columns = ['chr', 'start', 'end', 'start_base',
                                        'end_base', 'probe_seq', 'temp', 'm_length',
                                        'monomer', 'copy_num', 'r_length', 'probes',
                                        'ciel_copy_num', 'new_k', 'found', 'total',
                                        'Kmer_enrich', 'max_count_repeat','max_count_hg']

        print("Merging all information with total dataframe ... ")

temp_fasta = glob.glob('*_.fa')
for file in temp_fasta:
    os.remove(file)
