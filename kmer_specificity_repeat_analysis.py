#import libraries needed 

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
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein
import time
from collections import OrderedDict
import collections

#declare a timer
start_time=time.time()

#need to write a function that's going to know to open the chr_regions.fasta file
#need to strip that string in the last column of the dataframe and then use it to match against the fasta file

#probe_test="hg38_probes_df_2019_08_09_400.bed"

probe_test=sys.argv[1]
hg38_jf_index=sys.argv[2]

#clean the file remove first line and then first column 
#write a function that takes a small chunk of that dataframe with all the probe info
#generate a test file that at least has different probes and different, index k-mer positions

def parse_probe(probe_file,hg38_jf):
    probe_data=pd.read_csv(probe_file, delimiter = '\t',names=["chrom","p_start","p_end","probe_seq","Tm","scaffold"])
    probe_data[['chr','r_start','r_end']] = probe_data.scaffold.str.replace(":","-").str.split("-",expand=True)
    probe_data=probe_data.drop(columns=['chrom'])
    #print(probe_data)
    
    for i in probe_data['chr']:
        with open(i + "_regions.fa") as fasta:
            for rec in SeqIO.parse(fasta, "fasta"):
                id = rec.id
                seq = rec.seq
                id_file = open(id + "_split.fa", "w")
                id_file.write(">"+str(id)+"\n"+str(seq))
                id_file.close()
            fasta.close()
    
    names=[]
    split_fasta=glob.glob('*_split.fa')
    for i in split_fasta:
        names.append(i.split("_")[0])

    #while column is the same entry then generate a fasta file from the probe_seqs
    probes_list=[]
    split_probe_data = dict(tuple(probe_data.groupby('scaffold')))
    for key,value in split_probe_data.items():
        probe_id_start=value["p_start"]
        start=probe_id_start.values.tolist()
        probe_id_stop=value["p_end"]
        stop=probe_id_stop.values.tolist()
        seq=value["probe_seq"]
        sequence=seq.values.tolist()
        scaffold_id=value["scaffold"]
        scaffold=scaffold_id.values.tolist()
        chrom=value["chr"]
        chromosomes=chrom.values.tolist()
        
        for (a,b,c,d,e) in zip(chromosomes,start,stop,scaffold,sequence):
            split_probe_file=open(key+"_sep_probes.fa","a")
            split_probe_file.write(">"+str(a)+":"+str(b)+"-"+str(c)+"_"+str(key)+"\n"+str(e)+"\n")
            split_probe_file.close()
    
    #then if the scaffold matches a file name, do a fasta count for 18mers
        repeat_fasta=glob.glob(str(key)+"_split.fa")
        repeat_jf=str(key)+".jf"
        for item in repeat_fasta:
             subprocess.call(['jellyfish','count','-s', '3300M', '-m', '18',
                         '-o', repeat_jf, str(item)],stderr=None,shell=False)
        
    #now do a jf query on hg38 and each of those repeat regions
        repeat_jf_list=glob.glob(str(key)+".jf")
        split_probes=glob.glob(str(key)+"_sep_probes.fa")
        repeat_query_output=str(key)+"_repeat_query.txt"
        print(repeat_jf_list)
        for (a,b) in zip(split_probes,repeat_jf_list):
            print(a,b)
            subprocess.call(['jellyfish', 'query', str(b), '-s', 
                str(a), '-o', repeat_query_output], stderr=None, shell=False)
        
        hg_query_output=str(key)+"_hg_query.txt"
        for a in split_probes:
            subprocess.call(['jellyfish','query',hg38_jf,'-s',
                            str(a), '-o', hg_query_output])
                
        subprocess.call(['rm',str(key)+'_split.fa'],stderr=None,shell=False)
        subprocess.call(['rm',str(key)+'.jf'],stderr=None,shell=False)
        subprocess.call(['rm',str(key)+'_sep_probes.fa'],stderr=None,shell=False)

    repeat_queries=glob.glob('*_repeat_query.txt')

    hg_queries=glob.glob('*_hg_query.txt')

    probe_indices=glob.glob('*_probe_indices.txt')

    repeats=[]
    for item in repeat_queries:
        repeats.append(item.split("_")[0])

    hg_repeat=[]
    for item in hg_queries:
        hg_repeat.append(item.split("_")[0])

    probe_regions=[]
    for item in probe_indices:
        probe_regions.append(item.split("_")[0])

    repeats.sort()
    hg_repeat.sort()
    probe_regions.sort()

    repeat_contents=[]
    hg_contents=[]
    index_contents=[]

    for (a,b,c) in zip(repeats,hg_repeat,probe_regions):
        with open(str(a)+"_repeat_query.txt", 'r') as repeat_jf, open(str(b)+"_hg_query.txt", 'r') as hg_jf, open(str(c)+"_probe_indices.txt",'r') as probes:
            repeat_info=pd.read_csv(repeat_jf, delimiter = ' ',names=["kmer","count1"])
            hg_info=pd.read_csv(hg_jf, delimiter = ' ',names=["kmer","count2"])
            probe_info=pd.read_csv(probes, delimiter = ' ',names=["index"])

        df_col = pd.concat([repeat_info,hg_info,probe_info], axis=1)
        df_col['region']=a
        grouped=df_col.groupby(['index','region'])

        max_region=grouped['count1'].agg(np.max)
        max_hg=grouped['count2'].agg(np.max)

        merged_max=pd.concat([max_region, max_hg], axis=1)
        merged_max['region']=a
    
        merged_max.to_csv('indices_loc.txt', mode='a', header=False, index=False, sep="\t")
        
        subprocess.call(['rm','*_repeat_query.txt'],stderr=None,shell=False)
        subprocess.call(['rm','*_hg_query.txt'],stderr=None,shell=False)
        subprocess.call(['rm','_probe_indices.txt'],stderr=None,shell=False)

        index_file='indices_loc.txt'
    
        #probe_data=pd.read_csv(probe_test, delimiter = '\t',names=["chrom","p_start","p_end","probe_seq","Tm","region"])
        #probe_data[['chr','r_start','r_end']] = probe_data.region.str.replace(":","-").str.split("-",expand=True)
        #probe_data=probe_data.drop(columns=['chrom'])
        index_data=pd.read_csv(index_file,delimiter="\t",names=['repeat_count','hg_count','region'])

    merged_all=probe_data.merge(index_data, on="region", left_index=True, right_index=True)
    subprocess.call(['rm','indices_loc.txt'],stderr=None,shell=False)
    merged_all["k-mer enrich"]=merged_all['repeat_count']/(merged_all['hg_count']+1)
    print(merged_all)
    merged_all.to_csv('probes_final_score.txt', header=False, index=False, sep="\t")

parse_probe(probe_test,hg38_jf_index)
