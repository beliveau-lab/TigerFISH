#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
"""
Created on Thu Jul  1 19:09:52 2021
@author: Robin Aguilar
Beliveau and Noble Labs
University of Washington | Department of Genome Sciences
"""
##############################################################################

#import libraries
import pandas as pd
import matplotlib.pyplot as plt
import re
from sklearn import preprocessing
import time
import argparse

##############################################################################

def generate_dfs(chr_track_file,chr_overlap_bed,pDups_scores):
    """
    Function takes the genome bin file, the bedtools intersect file, and
    the pairwise file containing pdups computes to process them as 
    dataframes.
    
    Parameters
    ----------
    chr_track_file : file
        file containing the 1Mb bins of the genome
    chr_overlap_bed : file
        bedtools intersect output
    pDups_scores : file
        file from alignment output

    Returns
    -------
    chr_track : dataframe
        bed file containing the coordinates of 1Mb bins
    chr_overlap : dataframe
        contains the overlaps of bedtools intersect
    pairs_pdups : dataframe
        contains the pairwise output

    """

    #label columns for ground file with the bed tracks
    colnames = ["chrom","start","stop"]
    chr_track = pd.read_csv(chr_track_file, names=colnames,header=None,
                            sep='\t')
    
    chr_track=chr_track[~chr_track.chrom.str.contains("M")]

    #label columns for matching overlaps
    colnames = ["chrom", "start", "stop", "chrom_b","start_b","stop_b"]
    chr_overlap = pd.read_csv(chr_overlap_bed, names=colnames,header=None,
                          sep='\t')

    #label columns for data where you have the nupack data
    colnames = ["probe_type","parent","derived",
                "derived_coords","start","pDups"]

    pairs_pdups = pd.read_csv(pDups_scores, names=colnames,header=None,
                          sep='\t')

    return chr_track,chr_overlap,pairs_pdups

##############################################################################

def map_columns(pairs_pdups,chr_overlap):
    """
    Function maps the bedtools intersect overlaps with the pdups scorese

    Parameters
    ----------
    pairs_pdups : dataframe
        contains the pairwise probes with pdups computes
    chr_overlap : dataframe
        contains the genomic overlaps from bedtools intersect

    Returns
    -------
    bin_sums : dataframe
        contains the derived chrom, start, stop, and pdups score

    """

    pairs_pdups = pairs_pdups.drop(['parent', 'derived'], axis = 1)
    
    #overlap pairs_pdups to chr_overlap on o_chrom, o_start, o_stop
        
    pairs_pdups = pd.merge(pairs_pdups,chr_overlap, on='start')

    bin_sums = pairs_pdups.groupby(['chrom_b','start_b',
                           'stop_b'])['pDups'].sum().reset_index()
    
    #sort columns in ascending order
    bin_sums = bin_sums.sort_values(by=['chrom_b','start_b','stop_b'])
    
    
    bin_sums.columns = ['chrom', 'start','stop','pDups']
    
    return bin_sums

##############################################################################

def intersect_chr_tracks(bin_sums,chr_track):
    """
    Function takes the where probes map to genomic bins and their pdups
    scores and merges this to all other genomic bins. All other bins that
    do not have pdups scores mapping to them receive a 0.

    Parameters
    ----------
    bin_sums : dataframe
        contains chromosome tracks and corresponding pdups scores
    chr_track : dataframe
        contains the genomic chr, start, stop for 1Mb bins

    Returns
    -------
    merged : dataframe
        contains all genomic bins with corresponding mapped scores

    """
    
    #drop pdups for now 
    bins_w_vals = bin_sums.drop(['pDups'], axis = 1)

    #subset rows not in bins with vals
    common = bins_w_vals.merge(chr_track,on=['chrom', 'start','stop'])
    
    not_overlapping = pd.merge(chr_track, common, how='left', indicator=True) \
           .query("_merge == 'left_only'") \
           .drop('_merge',1)
           
    not_overlapping['pDups'] = 0.0
    
    #now merge with nupack
    frames = [bin_sums,not_overlapping]
    
    merged = pd.concat(frames)
    
    merged = merged.sort_values(by=['chrom','start','stop'])

    return merged

##############################################################################

def generate_plot(merged,out_plot):
    """
    Function takes the tracks of genomic regions and maps them
    accordingly based on scaffold.

    Parameters
    ----------
    merged : dataframe
        chrom, start, stop, and pdups sum

    Returns
    -------
    .png of in silico predicted binding of oligo signal
    """
    
    #sort merged
    chrom_list = merged['chrom'].tolist()
    
    #implement sort nicely
    def tryint(s):
        try:
            return int(s)
        except ValueError:
            return s

    def alphanum_key(s):
        return [tryint(c) for c in re.split('([0-9]+)', s)]

    def sort_nicely(l):
        return sorted(l, key=alphanum_key)
    
    chrom_list = sort_nicely(chrom_list)
    
    #keep the first instanc of chrom
    condensed_list = list(dict.fromkeys(chrom_list))
        
    chrom_nums = []
    for i in range(0,len(condensed_list)):
        chrom_nums.append(int(i))
        
    #make a dict that binds the chrom to the num
    chrom_dict = dict(zip(condensed_list, chrom_nums))
    
    merged['chrom_num'] = merged['chrom'].map(chrom_dict)
        
    merged = merged.sort_values(by=['chrom_num','start','stop'])
        
    #scales the values of pdups between 0 - 255 using normalization
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 255))
    y_transformed = min_max_scaler.fit_transform(merged[['pDups']])
    y_transformed_list = y_transformed.ravel()
    merged['nupack_trans'] = list(y_transformed_list)
    
    merged = merged.drop(['pDups'], axis=1)
    
    #generate linked lists for the chrom bins
    bins_ll=[]
    #values for transformed pdups
    y_vals_trans = []
    #chrom name
    standard_chrom = []
    
    #iterate by chromosome to append y transform values to scaffold
    grouped_probes=merged.groupby('chrom',sort=False)
    for name,group in grouped_probes:
        
        y_list = group['nupack_trans'].tolist()
        y_vals_trans.append(y_list)
        
        bins_list = group['start'].tolist()
        bins_ll.append(bins_list)
        standard_chrom.append(name)

    #plot figure size
    fig = plt.figure(figsize=(24, 8))
    
    #generate figure width
    fig_width_list = []
    for ll in bins_ll:
        fig_width_list.append((ll[-1]))
    
    #scale figure height
    height_list = []
    for i in range(len(fig_width_list)):
        height_list.append(24)
    
    #zip the bins and y vals
    for i,j,chrom,hw in zip(list(range(0,len(bins_ll))),
                  list(range(1,(len(y_vals_trans))+1)),
                  standard_chrom,fig_width_list):
        
        #make each chrom a subplot
        ax = fig.add_subplot(6,4,j)

        #for each figure make plotting params
        params = {"ytick.color" : "w",
          "xtick.color" : "w",
          "axes.labelcolor" : "w",
          "axes.edgecolor" : "w"}
        
        #update plot params
        plt.rcParams.update(params)
        
        #title of each plot is chrom name
        plt.title(str(chrom),color='w',fontsize=14,fontweight='bold')
        
        #set cmap colors
        cmap = plt.cm.Blues
        cmap.set_bad(color='w')

        #remove axes for easy viewing
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.spines['bottom'].set_linewidth(1.0)
        
        #plot as pileup of signal (histogram)
        plt.hist(bins_ll[i], weights= y_vals_trans[i],color = "magenta")
        ax.set_ylim([0,300])
        ax.set_xlim(0,hw)
                
    #organize subplots
    fig.subplots_adjust(bottom=0.10, right=0.6,top=0.9,hspace = 2.0,
                        wspace = 0.8)
    
    plt.savefig(out_plot,transparent=True,dpi=300)
    
    plt.show()
    
    
##############################################################################

def main():
    
    start_time=time.time()

    userInput = argparse.ArgumentParser(description=\
        'Requires the genome bin file, derived sequence bedtools intersect'
        'and the file containing pariwise sequences')
        
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-c_t', '--chrom_track', action='store', 
                               required=True, help='The genome file'
                               'containing chrom sizes')
    requiredNamed.add_argument('-c_o', '--chrom_overlaps', action='store',
                               required=True, help='bedtools intersect'
                               'output')
    requiredNamed.add_argument('-p', '--pairwise_pdups', action='store',
                               required=True, help='file with pairwise'
                               'pdups vals')
    requiredNamed.add_argument('-pl', '--out_plot', action='store',
                               required=True, help='output genome plot')

    args = userInput.parse_args()
    chr_track_file = args.chrom_track
    chr_overlap_bed = args.chrom_overlaps
    pdups_scores = args.pairwise_pdups
    out_plot = args.out_plot
    
    chr_track,chr_overlap,pairs_pdups = generate_dfs(chr_track_file,
                                                     chr_overlap_bed,
                                                     pdups_scores)

    print("---%s seconds ---"%(time.time()-start_time))

    
    bin_sums = map_columns(pairs_pdups,chr_overlap)

    print("---%s seconds ---"%(time.time()-start_time))
    
    merged = intersect_chr_tracks(bin_sums,chr_track)

    print("---%s seconds ---"%(time.time()-start_time))

    generate_plot(merged,out_plot)

    print("---%s seconds ---"%(time.time()-start_time))
    
    print("done")
    
if __name__== "__main__":
    main()
