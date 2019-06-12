#!/usr/bin/env probeMining

##################
###4/24/2019
###Robin Aguilar
###Beliveau lab
###repeatFinder.py
###Input: The files from the trf shell script output
###Output: A dataframe containing the coordinates of each repeat and it's corresponding designed probes as .fastq files (the output of blockparse.py)
###Dependencies: the blockparse.py script. It calls this as a subprocess. 
###Purpose:This script is used to organize and run blockparse on the repeats identified from the output of tandem repeat finder. The format of the dataframe is the following:
###chr,start,stop,r_length,copy_num,monomer,m_length
##################


import time
import io
import sys
import re
import csv
from collections import defaultdict
from itertools import islice
from operator import itemgetter, attrgetter
import blockParse
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
import numpy as np
import pandas as pd
import glob
import os

#you will want to make this an argument so it can read in any text file that's passed into it from snakemake
TRF_output= sys.argv[1]
#TRF_output="chrY.fa.2.7.7.80.10.50.500.dat"
chrom_names=TRF_output.split('.')[0]
sequence_num = 0
sequence_dict = {}
all_fa_sequences=[]
repeat_information=[]
HEADER=15
SPAN_L=3000
MONOMER_L=80
MIN_L=80
cat_files=[]

with open(TRF_output,'r') as file:
    for i, line in enumerate(file):
        #may need to remove this if statement if you use snakemake to remove the first 15 lines anyway
        if i >= HEADER:
            #replaces the spaces in the formatted text with a tab
            replace_space=line.replace(' ','\t')
            elements=replace_space.split()
            start=elements[0]
            stop=elements[1]
            copy_num=elements[3]
            span=int(elements[1])-int(elements[0])
            #if the difference between the start and the end is greater than 3000 bases
            if int(span) >= SPAN_L:
                #if the length of the monomer is >= 80 then, add that information from the repeat to the list
                if len(elements[13])>=MONOMER_L:
                    repeat_information.append(chrom_names+'\t'+str(elements[0]+'\t'+str(elements[1])+'\t'+str(span)) +'\t'+str(copy_num)+'\t'+str(elements[13])+'\t'+str(len(elements[13]))+'\n')
                    all_fa_sequences.append(elements[13])
                #if the length of the monomer is < 80 then you need to make the monomer at least 80 bases
                if len(elements[13])<MONOMER_L:
                    min_length=MIN_L
                    #this is the equation to make the monomer greater than 80, it just tacks on the bases from the start until 80 bases
                    edited_elements=(elements[13]*(min_length//len(elements[13])+1))
                    all_fa_sequences.append(edited_elements)
                    repeat_information.append(chrom_names+'\t'+str(elements[0]+'\t'+str(elements[1])+'\t'+str(span))+'\t'+str(copy_num)+'\t'+str(edited_elements)+'\t'+str(len(edited_elements))+'\n')
                    #write all elements in the list to a text file that contains start, stop, length of span, the monomer, and length of monomer
                    with open(chrom_names+'_out_sequences.txt','w') as f:
                        for i in repeat_information:
                            f.write("%s" % i)

#take the chromosome text file to pass into function
TRF_parsed=chrom_names+'_out_sequences.txt'
all_monomers=[]
fasta_list=[]
test_list=[]
with open(TRF_parsed,'r') as file:
    df_toparse=pd.read_csv(TRF_parsed,delimiter="\t")
    df_toparse.columns=['chr','start','end','r_length','copy_num','monomer','m_length']
    final_df=df_toparse.groupby((df_toparse.end.shift() - df_toparse.start).lt(0).cumsum()).agg(
        {'chr': 'first', 'start': 'first', 'end': 'last','r_length':'sum','copy_num':'sum','monomer':'sum','m_length':'sum'})
    print(final_df)
    final_df.to_csv(chrom_names+'_parsed_repeats.txt', header=None, index=None, sep='\t')

f_repeats=chrom_names+'_parsed_repeats.txt'
#read the text file as a pandas dataframe
df=pd.read_csv(f_repeats,delimiter="\t")

#defines a function that will iterate through blockParse as a module
def fasta_blockParse(fasta_line, output_name):
    #make a temp fasta file
	fasta = open("temp.fa", "w")
    #write the contents in a fasta format
	fasta.write(">temp\n%s" % (fasta_line))
	fasta.close()
	# run blockParse on temp fasta
	blockParse.runSequenceCrawler("temp.fa", 36, 41, 20, 80, mt.DNA_NN3, 42, 47, 'AAAAA,TTTTT,CCCCC,GGGGG', 390, 50, 0, 25, 25, None, True , False, False, False, False, False, output_name)
    #will remove all the temp files
    	subprocess.call(['rm','temp.fa'],stderr=None,shell=False)


    #will iterate through the length of the dataframe
for i in range(0, len(df), 1):
    #will take the chr, start, stop as the name of the output
	outname=str(df.iloc[i,3])+"_"+str(df.iloc[i,2])+"_"+str(df.iloc[i,6])
	#will take the output above and write the sequence to run through the blockParse sequenceCrawler
	fasta_blockParse(df.iloc[i,1],outname)

'''
#reads in all files in the current directory that end in *.fastq
read_files=glob.glob("*[0-9].fastq")

#will take all the current *.fastq in the directory and concatenate them together as one outfile
with open(chrom_names+"_result.fastq","wb") as outfile:
	for f in read_files:
		with open(f,"rb") as infile:
			outfile.write(infile.read()+"\n")

#finds any *.fastq file that matches the guidelines of having end position numbers (these are the files you want to remove)
filelist=glob.glob("*[0-9].fastq")
for file in filelist: 
	os.remove(file)

#now you will be left with only one chromosome output file that contains the concat output together from blockParse
'''

