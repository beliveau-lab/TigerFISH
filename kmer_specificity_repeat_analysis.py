"""
Author: Robin Aguilar
Create date: 07/17/2019, Python3.7
Date last modified: 08/25/2019
Beliveau and Noble Labs
Dependencies: Jellyfish
Purpose: this script implements a kmer count approach for identifying regions of interest containing repeat elements
This script should be able to keep track of indices from a fasta file and provide them for kmers identified from .jf
"""
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
from collections import OrderedDict
import collections

#import modules needed for blockParse refactor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from itertools import groupby
from Bio.Alphabet import generic_dna, generic_protein

#declare a timer
start_time=time.time()

########################################################################################################################
"""
The purpose of this function is to parse the probe file, perform a jf query for probes dedigned against a k-mer enrich region
on that ID'd region. This script then calculates the k-mer enrichment score, and filters probes above thersholds.
"""

def parse_probe(probe_file,hg38_jf,probe_fa):
    probe_data=pd.read_csv(probe_file, delimiter = '\t',names=["index","p_start","p_end","probe_seq","Tm","region","chr","r_start","r_end"])
    probe_data[['chr','r_start','r_end']] = probe_data.region.str.replace(":","-").str.split("-",expand=True)
    probe_data=probe_data.drop(columns=['index'])
    
    with open(probe_fa) as fasta:
        for rec in SeqIO.parse(fasta, "fasta"):
            id = rec.id
            seq = rec.seq
            id_file,fasta_id_path=tempfile.mkstemp(suffix="_split.fa", prefix=id+"_", dir=".")
            id_file= os.fdopen(id_file,"w")
            id_file.write(">"+ str(id) + "\n" + str(seq))
            id_file.close()

########################################################################################################################

    #in dataframe, group by region
    #write a temp file for each k-mer enriched region that contains all probes designed against that region
    
    probes_list=[]
    split_probe_data = dict(tuple(probe_data.groupby('region')))
    for key,value in split_probe_data.items():
        seq=value["probe_seq"]
        sequence=seq.values.tolist()
        seqList=[]
        temp_probe_file, probe_file_path=tempfile.mkstemp(suffix="_sep_probes.fa", prefix=key+"_", dir=".")
        with os.fdopen(temp_probe_file, 'a') as tmp:
            for i in sequence:
                tmp.write('>' + str(key) + "\n" + str(i) + "\n")
                seqList.append(i)
            tmp.close()
     
    #make a file, that stores the indices of each probe when broken into 18mers for that region
    #this becomes useful for appending the correct k-mer enrichment score to a particular probe in a given region
    
        indexList = []
        for i in range(0, len(seqList), 1):
            probeLength = len(seqList[i])
            merWindows = probeLength - merLengthVal + 1
            for j in range(0, merWindows, 1):
                indexList.append('%d' % i)
        index_file, index_file_path=tempfile.mkstemp(suffix="_probe_indices.txt", prefix=key+"_", dir=".")
        index_file= os.fdopen(index_file,"a")
        for item in indexList:
            index_file.write("%s\n" % item)
        index_file.close()
 
########################################################################################################################

    #generate a list the k-mer enriched fasta regions
    repeat_probes=glob.glob("*_sep_probes.fa")
    repeat_fasta=glob.glob("*_split.fa")

    #generate a tmp directory to remove all results at end
    dirpath = tempfile.mkdtemp(dir=".")

########################################################################################################################

    #generate names for the .jf files you will make from jellyfish
    jf_name=[]
    for item in repeat_fasta:
        jf_name.append(dirpath+"/"+str(item.split("_")[0]) + ".jf")

    #sort the names of the k-mer enriched fastas, and the jf files. So files with same header can be iterated together
    jf_name.sort()
    repeat_fasta.sort()
    repeat_probes.sort()
    
    #zips the two lists above, then creates .jf index files for each k-mer enrich. region, puts in tmp dir
    for (jf, fasta) in zip(jf_name,repeat_fasta):
        print(jf,fasta)
        subprocess.call(['jellyfish','count','-s', '3300M', '-m', '18','-o', jf, '-L', '2', str(fasta)],stderr=None,shell=False)

########################################################################################################################

    #generate names for the probe_repeat_query performed with jf query
    repeat_name=[]
    for item in repeat_fasta:
       repeat_name.append(dirpath+"/"+str(item.split("_")[0]) + "_repeat_query.txt")
    
    #sort the items in this list, so it matches the sort of jf_name and repeat_fasta
    repeat_name.sort()

    #zip all three lists together to run a jf query - decomposes probes into 18mers, returns count of 18mer in k-mer enr. region
    for (probes,jf,name) in zip(repeat_probes,jf_name,repeat_name):
        print(probes,jf,name)
        subprocess.call(['jellyfish', 'query', str(jf), '-s', str(probes), '-o', str(name)], stderr=None, shell=False)
        
########################################################################################################################
    
    #generate names for the genome_wide query performed with jf query
    hg_name=[]
    for item in repeat_fasta:
        hg_name.append(dirpath+"/"+str(item.split("_")[0])+"_hg_query.txt")

    #sort these names so the ordering matches those in repeat_fasta
    hg_name.sort()

    #decomposes probes into 18mers, returns count of 18mer in hg38 query
    for (probes,hg) in zip(repeat_probes,hg_name):
        print(probes,hg)
        subprocess.call(['jellyfish','query',hg38_jf,'-s',str(probes), '-o', str(hg)])
        
########################################################################################################################
    
    #removes files no longer needed
    for item in glob.glob('*_split.fa'):
         os.remove(item)

    for item in glob.glob('dirpath/*.jf'):
        os.remove(item)

    for item in glob.glob('*_sep_probes.fa'):
        os.remove(item)

    #takes the probe indices generated earlier
    probe_indices=glob.glob('*_probe_indices.txt')

    #sorts repeat queries, hg queries, and probe_indices
    repeat_name.sort()
    hg_name.sort()
    probe_indices.sort()
    
    #zips list together, to open each and read into individual dataframes
    for (a,b,c) in zip(repeat_name,hg_name,probe_indices):
        with open(str(a), 'r') as repeat_jf, open(str(b), 'r') as hg_jf, open(str(c),'r') as probes:
            repeat_info=pd.read_csv(repeat_jf, delimiter = ' ',names=["kmer","r_count"])
            hg_info=pd.read_csv(hg_jf, delimiter = ' ',names=["kmer","h_count"])
            probe_info=pd.read_csv(probes,names=["index"])

            #keep track of id_regions in probe_indices
            id_regions=[]
            for item in probe_indices:
                id_regions.append(item.split("_")[0])

        #concat into one dataframe
        df_col = pd.concat([repeat_info,hg_info,probe_info], axis=1)
        
        for name in id_regions:
            df_col['region']=name
        grouped=df_col.groupby(['index','region'])

        #take the max k-mer in each group for the comparison (aggregate by max)
        max_region=grouped['r_count'].agg(np.max)
        max_hg=grouped['h_count'].agg(np.max)

        #merge indices with which to perform the k-mer enrich computation
        merged_max=pd.concat([max_region, max_hg], axis=1)
        for name in id_regions:
            merged_max['region']=name
        merged_max.to_csv('indices_loc.txt', mode='a', header=False, index=False, sep="\t")

        index_file='indices_loc.txt'
        index_data=pd.read_csv(index_file,delimiter="\t",names=['repeat_count','hg_count','region'])
        
########################################################################################################################

    #using the index file, now calculate k-mer enrich score
    merged_all=probe_data.merge(index_data, on="region", left_index=True, right_index=True)
    
    #calculate k-mer enrich score, need to add one because of how hg_file was made, does not count itself
    merged_all["k-mer enrich"]=merged_all['repeat_count']/(merged_all['hg_count']+1)
    merged_all=merged_all[(merged_all['k-mer enrich'] >= 0.50) & (merged_all['repeat_count'] >= 10)]

    #clear all redundant/temp files
    text_files=glob.glob('*.txt')
    for files in text_files:
        subprocess.call(['rm',files],stderr=None,shell=False)
    shutil.rmtree(dirpath)
    
    #keep final file of probes
    merged_all.to_csv('probes_final_score.txt', header=False, index=False, sep="\t")
    
########################################################################################################################
# MAIN
########################################################################################################################

def main:
    global USAGE:
        
        probe_test=sys.argv[1]
        hg38_jf_index=sys.argv[2]
        probe_fasta=sys.argv[3]

        merLengthVal=18
        
        parse_probe(probe_test,hg38_jf_index,probe_fasta)
        print("---%s seconds ---"%(time.time()-start_time))
        
