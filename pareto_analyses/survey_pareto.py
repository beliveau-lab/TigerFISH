"""
Created by Robin Aguilar
Beliveau and Noble Labs
Date Created: 09.22.2020
Purpose: To analyze all possibilities for the Pareto frontier (modifiable parameters in first two steps),
to see how these modify the total number of satellite DNAs identified by Tigerfish
"""

#load libraries
import time
import io
import sys
import re
import csv
import os
import numpy as np
import pandas as pd
import glob 
import itertools
from os import listdir
from os.path import isfile, join
import subprocess
from subprocess import check_output
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import combinations

##############################################################################################

def parse_file_names(step_2_file,step_4_file,step_2_f_name,step_4_f_name):

    global step_2_name
    global step_4_name

    base_2 = os.path.basename(step_2_file)
    step_2_name = os.path.splitext(base_2)[0]

    base_4 = os.path.basename(step_4_file)
    step_4_name = os.path.splitext(base_4)[0]

    #store the probe names in the lists to make the final proportion dataframe
    step_2_f_name.append(step_2_name)
    step_4_f_name.append(step_4_name)

    return step_2_name, step_4_name
    
##############################################################################################

def write_repeat_masker_to_bed(repeat_masker_file):
    
    colnames=["score","perc_div","perc_del","perc_ins",'chrom',"start","end",
              "left","strand", "sat_type","repeat_class","position_begin",
              "position_end","repeat_left","ID"]
    
    #open the ground probe file and read
    #trick - clean up by removing the headers using awk, then upload using whitespace delim
    rm_df = pd.read_csv(repeat_masker_file, names=colnames, header=None, delim_whitespace=True)

    rm_df=rm_df.drop(['score', 'perc_div','perc_del','perc_ins','left','strand',
                 'position_begin','position_end','repeat_left','ID'], axis=1)
    
    #total_number of satellites
    satellite_df = rm_df[rm_df['repeat_class'].str.contains("Satellite")]

    #total number of everything but satellite
    remaining_repeat_df = rm_df[~rm_df['repeat_class'].str.contains("Satellite")]

    #remove rows with alt scaffolds
    satellite_df = satellite_df[~satellite_df.chrom.str.contains("_")]

    #keep a satellite_df with the repeat_labels
    satellite_repeat_labels = satellite_df

    #drop the satellite column
    satellite_df = satellite_df.drop(['repeat_class'], axis=1)

    satellite_df.to_csv('satellite_regions.bed',header=None, index=None, sep='\t')

    #contains two dfs - the satellites only, all other repeats
    return satellite_df, remaining_repeat_df, satellite_repeat_labels 

##############################################################################################

def load_probes_designed_file(design_probe_file,step_2_t_probes_list,step_2_name):
        
    colnames=["chrom","start","end","probe","Tm","region"]
    
    #open the ground probe file and read
    designed_probes = pd.read_csv(design_probe_file, names=colnames, header=None, sep = '\t')

    #this will append the total number of probes at this step
    designed_probes=designed_probes.drop_duplicates(subset=['probe'])

    step_2_t_probes_list.append(len(designed_probes))

    designed_probes=designed_probes.drop(["chrom","start","end","probe","Tm"], axis=1)

    #remove duplicate regions and parse columns into bed file ~48k regions
    designed_probes=designed_probes.drop_duplicates(subset=['region'])

    designed_probes[['chr','start_end']] = designed_probes['region'].astype(str).str.split(':', expand=True)
    designed_probes[['start', 'end']] = designed_probes['start_end'].astype(str).str.split('-', expand=True)

    designed_probe_regions=designed_probes.drop(['region','start_end'], axis=1)

    probes_2_name = "step_2_regions/" +  str(step_2_name) + ".bed"
    
    #make a bed file here
    designed_probe_regions.to_csv(probes_2_name,header=None, index=None, sep='\t')

    return designed_probe_regions,probes_2_name

##############################################################################################

def load_probes_post_lda_file(post_lda_file,step_4_t_probes_list,step_4_name):
    
    colnames=["region","probe","Tm","start_index","end_index","h_count","r_count","pb_score",
              "bt2_index"]
    
    #open the ground probe file and read
    post_lda_probes = pd.read_csv(post_lda_file, names=colnames, header=None, sep = '\t')
    
    post_lda_probes=post_lda_probes.drop(["probe","Tm","start_index","end_index","h_count",
                                          "r_count","pb_score","bt2_index"], axis=1)

    #remove duplicate regions and parse columns into bed file ~48k regions
    post_lda_probes=post_lda_probes.drop_duplicates(subset=['region'])

    #this will append the total number of probes at this step
    step_4_t_probes_list.append(len(post_lda_probes))
    
    post_lda_probes[['chr','start_end']] = post_lda_probes['region'].astype(str).str.split(':', expand=True)
    post_lda_probes[['start', 'end']] = post_lda_probes['start_end'].astype(str).str.split('-', expand=True)

    post_lda_probe_regions=post_lda_probes.drop(['region','start_end'], axis=1)

    probes_4_name = "step_4_regions/" + str(step_4_name) + ".bed"

    #make a bed file here
    post_lda_probe_regions.to_csv(probes_4_name,header=None, index=None, sep='\t')

    return post_lda_probe_regions, probes_4_name

##############################################################################################

def parse_subprocess(sub_obj):

    #used to parse out the subprocess output, might want to make it, it's own function
    if sys.version_info[0] < 3:
        from StringIO import StringIO
    else:
        from io import StringIO

    subp_string = StringIO(sub_obj.communicate()[0].decode('utf-8'))

    return subp_string

##############################################################################################

def run_bedtools(satellite_DNA_file, probes_2_name , probes_4_name , satellite_repeat_labels, step_2_name, step_4_name):

    step_2_bedtools = subprocess.Popen(['bedtools','intersect','-wa','-wb','-a',satellite_DNA_file,'-b',probes_2_name],stdout=subprocess.PIPE)

    step_2_to_df = parse_subprocess(step_2_bedtools)

    #make column names for the satellites found dataframe
    colnames=["chrom","start","end",'sat_type',"chr_r","r_start","r_end"]

    #generate the dataframe
    satellites_step_2 = pd.read_csv(step_2_to_df, names = colnames, header = None, sep="\t")

    #now we want to say that if the sat_start and sat_end match with the values in the start and end as repeat
    satellites_step_2 = satellites_step_2.merge(satellite_repeat_labels, how='left')

    step_4_bedtools = subprocess.Popen(['bedtools','intersect','-wa','-wb','-a',satellite_DNA_file,'-b',probes_4_name],stdout=subprocess.PIPE)

    step_4_to_df = parse_subprocess(step_4_bedtools)

    #generate the dataframe
    satellites_step_4 = pd.read_csv(step_4_to_df, names = colnames, header = None, sep="\t")

    #now we want to say that if the sat_start and sat_end match with the values in the start and end as repeat
    satellites_step_4 = satellites_step_4.merge(satellite_repeat_labels, how='left')

    #write the satellites you're going to keep
    probes_2_name = "kept_sats_2/" + str(step_2_name) + "_kept_sats.bed"
    satellites_step_2.to_csv(probes_2_name,header=None, index=None, sep='\t')

    probes_4_name = "kept_sats_4/" + str(step_4_name) + "_kept_sats.bed"
    satellites_step_4.to_csv(probes_4_name,header=None, index=None, sep='\t')

    return satellites_step_2, satellites_step_4

##############################################################################################

def remove_non_pf_values(param_df):

    param_df.reset_index(inplace = True) 
 
    #prop_list
    prop_list = param_df['proportion'].tolist()
    print("Entire list")
    print(prop_list)

    non_pareto_vals = []

    n = len(prop_list)

    brr = [0]*n; l = 1;  
  
    brr[0] = prop_list[0];  
    for i in range(1,n): 
        if (brr[l - 1] <= prop_list[i]): 
            brr[l] = prop_list[i];  
            l += 1  
  
    # Print the sorted array  
    for i in range(l) : 
        non_pareto_vals.append(brr[i]); 

    param_df = param_df[param_df['proportion'].isin(non_pareto_vals)]
 
    return param_df

##############################################################################################

def analyze_repeat_population(step_2_f_name,step_4_f_name,step_2_sat_p,step_4_sat_p,step_2_t_probes_list,step_4_t_probes_list,step_2_alphas,step_4_alphas):

    #lets make a dataframe with the file name, probe proportion, probe total
    all_satellite_2_df = pd.DataFrame(list(zip(step_2_f_name, step_2_sat_p, 
                                         step_2_t_probes_list)), columns =['params', 'proportion', 'probes'])
    alpha_satellite_2_df = pd.DataFrame(list(zip(step_2_f_name, step_2_alphas,
                                         step_2_t_probes_list)), columns =['params', 'proportion', 'probes'])

    all_satellite_4_df = pd.DataFrame(list(zip(step_4_f_name, step_4_sat_p, 
                                         step_4_t_probes_list)), columns =['params', 'proportion', 'probes'])
    alpha_satellite_4_df = pd.DataFrame(list(zip(step_4_f_name, step_4_alphas,  
                                         step_4_t_probes_list)), columns =['params', 'proportion', 'probes'])

    #make the satellite_p and alpha_p columns in ascending order
    #rank in ascending order of probes
    all_satellite_2_df = all_satellite_2_df.sort_values(by=['probes'],ascending=True)
    all_satellite_4_df = all_satellite_4_df.sort_values(by=['probes'],ascending=True)

    #alpha sort
    alpha_satellite_2_df = alpha_satellite_2_df.sort_values(by=['probes'],ascending=True)
    alpha_satellite_4_df = alpha_satellite_4_df.sort_values(by=['probes'],ascending=True)

    #run the pareto function
    all_satellite_2_df = remove_non_pf_values(all_satellite_2_df)
    all_satellite_4_df = remove_non_pf_values(all_satellite_4_df)

    alpha_satellite_2_df = remove_non_pf_values(alpha_satellite_2_df)
    alpha_satellite_4_df = remove_non_pf_values(alpha_satellite_4_df)

    #label new column to describe data set
    all_satellite_2_df['type'] = "all_satellite"
    all_satellite_4_df['type'] = "all_satellite"

    alpha_satellite_2_df['type'] = "alpha_satellite"
    alpha_satellite_4_df['type'] = "alpha_satellite"

    #combine the dataframes together
    satellite_step_2_df = all_satellite_2_df.append(alpha_satellite_2_df, ignore_index=True)
    satellite_step_4_df = all_satellite_4_df.append(alpha_satellite_4_df, ignore_index=True)

    step_2_file = "step_2_stats_summary.html"
    #write the satellites you're going to keep
    probes_2_name = "stats_2/" + str(step_2_file)
    satellite_step_2_df.to_html(probes_2_name)

    step_4_file = "step_4_stats_summary.html"
    probes_4_name = "stats_4/" + str(step_4_file)
    satellite_step_4_df.to_html(probes_4_name)

    return satellite_step_2_df, satellite_step_4_df

##############################################################################################

def survey_missing_satellites(satellites_step_2,satellites_step_4,satellite_repeat_labels,step_2_sat_p,step_4_sat_p, step_2_name, step_4_name):

    #drop the columns that are not needed to do an anti-intersect to find the missing satellites
    captured_satellites_2 = satellites_step_2.drop(['chr_r', 'r_start','r_end'], axis=1)

    #this merges the satellites you have with the total labeled set to find those that were not id'd
    merged = satellite_repeat_labels.merge(captured_satellites_2, how='left', indicator=True)
    merged = merged[merged['_merge']=='left_only']

    #drop the merged column and rename dataframe    
    missing_satellites_2 = merged.drop(['_merge'], axis=1)
    
    #you need to get the proportions here
    step_2_sat_p.append(((len(satellite_repeat_labels)-len(missing_satellites_2))/len(satellite_repeat_labels)))

    #drop the columns that are not needed to do an anti-intersect to find the missing satellites
    captured_satellites_4 = satellites_step_4.drop(['chr_r', 'r_start','r_end'], axis=1)

    #this merges the satellites you have with the total labeled set to find those that were not id'd
    merged = satellite_repeat_labels.merge(captured_satellites_4, how='left', indicator=True)
    merged = merged[merged['_merge']=='left_only']

    #drop the merged column and rename dataframe    
    missing_satellites_4 = merged.drop(['_merge'], axis=1)

    #you need to get the proportions here
    step_4_sat_p.append(((len(satellite_repeat_labels)-len(missing_satellites_4))/len(satellite_repeat_labels)))

    #generate a file that returns the missing satellites as a file
    probes_2_name = "missing_sats_2/" + str(step_2_name) + "_missing_sats.bed"
    missing_satellites_2.to_csv(probes_2_name,header=None, index=None, sep='\t')

    probes_4_name = "missing_sats_4/" + str(step_4_name) + "_missing_sats.bed"
    missing_satellites_4.to_csv(probes_4_name,header=None, index=None, sep='\t')  


##############################################################################################

def plot_pareto_2(satellite_step_2_df):

    fig, ax = plt.subplots()

    for key, grp in satellite_step_2_df.groupby(['type']):
        ax = grp.plot(ax=ax, x='probes',
                      drawstyle='steps-pre', y='proportion', label=key)

    plt.legend(loc='best')

    plt.show()
    plt.savefig('pareto_plots/sat_2_pareto.png')
    plt.clf()

##############################################################################################

def plot_pareto_4(satellite_step_4_df):

    fig, ax = plt.subplots()

    for key, grp in satellite_step_4_df.groupby(['type']):
        ax = grp.plot(ax=ax, x='probes',
                      drawstyle='steps-pre', y='proportion', label=key)

    plt.legend(loc='best')
    plt.show()
    plt.savefig('pareto_plots/sat_4_pareto.png')
    plt.clf()

##############################################################################################

def identify_alphas(step_2_alphas,step_4_alphas,satellite_repeat_labels,satellite_step_2,satellite_step_4):

    total_alphas = satellite_repeat_labels[satellite_repeat_labels['sat_type'].str.contains("Alpha")]
    total_alphas = total_alphas.drop_duplicates(subset=['chrom','start','end'])

    #in satellite_repeat_labels find the total number of unique alphas total
    alphas_2 = satellite_step_2[satellite_step_2['sat_type'].str.contains("Alpha")]
    alphas_4 = satellite_step_4[satellite_step_4['sat_type'].str.contains("Alpha")]

    #then find number of unique alphas in step 2 and step 4
    alphas_2 = alphas_2.drop_duplicates(subset=['chrom','start','end'])
    alphas_4 = alphas_4.drop_duplicates(subset=['chrom','start','end'])

    #then appen this value to the list
    step_2_alphas.append((len(alphas_2))/len(total_alphas))
    step_4_alphas.append((len(alphas_4))/len(total_alphas))

    return step_2_alphas, step_4_alphas

##############################################################################################

def main():
    
    #start your timer
    start_time=time.time()

    repeat_masker_file="../../data/2020_09_30_hg38_default_probes/hg38_repeat_masker_nh.bed"
    satellite_DNA_file = "satellite_regions.bed"

    #whole_dataset
    step_2_files = [file for file in glob.glob("../../data/2020_10_05_survey_pareto/step_2_probes/*.bed")]
    step_4_files = [file for file in glob.glob("../../data/2020_10_05_survey_pareto/step_4_probes/*.bed")]

    #test_files
    #step_2_files = [file for file in glob.glob("../../data/2020_09_30_hg38_default_probes/design_probes/*.bed")]
    #step_4_files = [file for file in glob.glob("../../data/2020_09_30_hg38_default_probes/post_lda/*.bed")]

    #you need to store the file names in the order that they come in to make a dataframe
    step_2_f_name = []
    step_4_f_name = []

    #you need to store the total number of probes in each step
    step_2_t_probes_list = []
    step_4_t_probes_list = []

    #you need to store the proportion of the satellites identified
    step_2_sat_p = []
    step_4_sat_p = []

    #you need to store the proortion of the alphas identified
    step_2_alphas = []
    step_4_alphas = []

    #while the stems of the two files are the same
    for f1,f2 in zip(step_2_files,step_4_files):

        step_2_name,step_4_name = parse_file_names(f1,f2,step_2_f_name,step_4_f_name)

        satellite_df, remaining_repeat_df, satellite_repeat_labels = write_repeat_masker_to_bed(repeat_masker_file)

        designed_probe_regions,probes_2_name = load_probes_designed_file(f1,step_2_t_probes_list,step_2_name)

        post_lda_probe_regions,probes_4_name = load_probes_post_lda_file(f2,step_4_t_probes_list,step_4_name)

        satellites_step_2, satellites_step_4 = run_bedtools(satellite_DNA_file, probes_2_name , probes_4_name , satellite_repeat_labels, step_2_name,step_4_name)

        survey_missing_satellites(satellites_step_2,satellites_step_4,satellite_repeat_labels, step_2_sat_p, step_4_sat_p, step_2_name, step_4_name)

        step_2_alphas, step_4_alphas = identify_alphas(step_2_alphas,step_4_alphas,satellite_repeat_labels,satellites_step_2,satellites_step_4)

    print("---%s seconds ---"%(time.time()-start_time))

    satellite_step_2_df, satellite_step_4_df = analyze_repeat_population(step_2_f_name,step_4_f_name,step_2_sat_p,step_4_sat_p,step_2_t_probes_list,step_4_t_probes_list,step_2_alphas,step_4_alphas)

    print("---%s seconds ---"%(time.time()-start_time))
        
    plot_pareto_2(satellite_step_2_df)

    print("---%s seconds ---"%(time.time()-start_time))

    plot_pareto_4(satellite_step_4_df)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
