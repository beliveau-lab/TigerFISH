"""
Robin Aguilar
Beliveau and Noble labs
Created: 05.09.2021
Purpose: In this script the purpose of it is to consolidate and run bt2 by RR
after this, sam2 pairwise is implemented to return alignments,
these alignments are then checked to see if all alignments returned are < k
if so, proceed with storing region as dataframe with alignments
consolidate seqs to run nupack and filter as necessary, then store df
The returned file should summarize the probes you want to keep and their 
corresponding alignments.
The script after this should just map the returned probes to repeat regions
from the chm13 genome and a summary file that reports "done".
"""

import argparse
import time
import pandas as pd
import nupack
from nupack import *
import tempfile
from Bio.Seq import reverse_complement as rev_comp
import subprocess
from subprocess import Popen
import itertools
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import GC
from Bio.Seq import Seq

###################################################################################

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--probe_file', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')
requiredNamed.add_argument('-o', '--out_file', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')
args = userInput.parse_args()
p_file = args.probe_file
o_file = args.out_file

bowtie_string = "--local -N 1 -L 20 -R 3 -D 20 -i C,4 --score-min G,1,4"
bt2_idx = "/net/beliveau/vol1/home/eaguil/tigerfish/bt2_indices/t2t/t2t"

# configure nupack
NUPACK_MODEL = nupack.Model(
    material = 'dna', 
    celsius = 74.5, 
    sodium = 0.39, 
    magnesium = 0.0,
    ensemble = 'stacking')

strand_conc_a = 1e-6
strand_conc_b = 1e-12

region_threshold = 5000

#default will be 0.95
pdups_prop = 0.95

start_time=time.time()

###################################################################################

def read_probe_filter(p_file):
    
    #read in the file with the appropriate column names
    colnames = ["chrom","p_start","p_end","probe","Tm","region","r_count","h_count","max_mer","k_score"]
    
    #read as a dataframe
    probe_df = pd.read_csv(p_file, delimiter = '\t', names = colnames)

    probe_df[['chrom','region']] = probe_df['region'].str.split(':',
                expand=True)

    probe_df[['r_start','r_end']] = probe_df['region'].str.split('-',
                expand=True)

    #cast r columns as ints
    probe_df[['r_start','r_end']] = probe_df[['r_start',
                'r_end']].astype(int)

    #make a column for chrom:start-stop that matches true genomic coords
    probe_df['start'] = probe_df['p_start'] + probe_df['r_start']
    probe_df['end'] = probe_df['p_end'] + probe_df['r_start']

    probe_coords_list = []
    #make lists to make a new region
    chrom_list = probe_df['chrom'].tolist()
    start_list = probe_df['start'].tolist()
    end_list = probe_df['end'].tolist()

    #make the probe coords column
    for chrom,start,end in zip(chrom_list,start_list,end_list):
        probe_coords_list.append('%s:%s-%s' % (chrom,start,end))

    probe_df['probe_coords'] = probe_coords_list

    #now do the same for region again
    r_start_list = probe_df['r_start'].tolist()
    r_end_list = probe_df['r_end'].tolist()

    r_coords_list = []
    for chrom,r_start,r_end in zip(chrom_list,r_start_list,r_end_list):
        r_coords_list.append('%s:%s-%s' % (chrom,r_start,r_end))

    probe_df['region'] = r_coords_list

    #cull columns
    columns = ['chrom','p_start','p_end','r_start','r_end','start','end','max_mer']
    probe_df = probe_df.drop(columns,axis=1)

    #rearrange cols
    probe_df = probe_df[['probe_coords','region','probe','Tm','r_count','h_count','k_score']]

    return probe_df

###################################################################################

def filter_by_region(probe_df):

    #need to make lists to store the probe names, on target, off target, prop
    keep_probe_names_list = []
    keep_on_target_list = []
    keep_off_target_list = []
    keep_probe_prop_list = []
    keep_probe_seq = []

    #initiate groupby
    grouped_region=probe_df.groupby('region',sort=False)

    for name,group in grouped_region:

        #this takes into account single instances of probe in repeat
        probe_list = group['probe'].tolist()
        probe_coords_list = group['probe_coords'].tolist()

        #make a list of the repeat regions
        probe_regions_list = group['region'].tolist()

        probe_region_dict = dict(zip(probe_coords_list,probe_regions_list))

        threshold_count = 0
         
        while threshold_count <= region_threshold and len(probe_list) >= 1:

            #make the call for the top probe to get the pairwise df
            top_probe_pairwise = generate_pairwise_df(probe_list[0],probe_coords_list[0])

            #compute the on target sum for the top probe
            total_prop_dict,total_on_target_dict,total_off_target_dict = target_nupack_sum(top_probe_pairwise,probe_region_dict)

            for key,val in total_prop_dict.items():
                if val >= pdups_prop and len(keep_probe_names_list) == 0:

                    keep_probe_names_list.append(key)
                    keep_on_target_list.append(total_on_target_dict[key])
                    keep_off_target_list.append(total_off_target_dict[key])
                    keep_probe_prop_list.append(total_prop_dict[key])
                    keep_probe_seq.append(probe_list[0])
                    
                    threshold_count += total_on_target_dict[key]

                if val >= pdups_prop and len(keep_probe_names_list) >= 1:
                    pdups_track_list = []
                    for probe in keep_probe_seq:
                        pdups_compare = pdups(probe_list[0],probe)
                        pdups_track_list.append(pdups_compare)

                    if max(pdups_track_list) <= 0.25:

                        keep_probe_names_list.append(key)
                        keep_on_target_list.append(total_on_target_dict[key])
                        keep_off_target_list.append(total_off_target_dict[key])
                        keep_probe_prop_list.append(total_prop_dict[key])
                        keep_probe_seq.append(probe_list[0])
                        threshold_count += total_on_target_dict[key]

            probe_list.pop(0)
            probe_coords_list.pop(0)

    #make dictionaries of the on target, off target, prop
    on_target_dict = dict(zip(keep_probe_names_list,keep_on_target_list))
    off_target_dict = dict(zip(keep_probe_names_list,keep_off_target_list))
    prop_target_dict = dict(zip(keep_probe_names_list,keep_probe_prop_list))

    return on_target_dict,off_target_dict,prop_target_dict

###################################################################################

def generate_final_df(probe_df,on_target_dict,off_target_dict,prop_target_dict):

    #make the three dicts into a dataframe
    on_target_df = pd.DataFrame(on_target_dict.items(), columns=['probe_coords', 'on_target_sum'])
    off_target_df = pd.DataFrame(off_target_dict.items(), columns=['probe_coords', 'off_target_sum'])
    prop_target_df = pd.DataFrame(prop_target_dict.items(), columns=['probe_coords', 'on_target_prop'])

    #then merge the dicts into the main probe df (lengths should match)
    on_off_merged = pd.merge(on_target_df, off_target_df, on="probe_coords")
    pdups_binding_df = pd.merge(on_off_merged, prop_target_df, on="probe_coords")

    #write final df
    keep_probes_df = pd.merge(probe_df, pdups_binding_df, on="probe_coords")

    keep_probes_df.to_csv(o_file, header=False, index=False, sep="\t")

###################################################################################

def generate_pairwise_df(probe_seq,probe_coords):

        #begin looping over each RR to write a tmp file
        with tempfile.TemporaryDirectory() as tmpdir:

            #make temporary file name
            fastq_filename = tmpdir + '/region.fastq'

            #reads the file to start writing
            fastq = open(fastq_filename,'w')

            #writes the file as a temp fastq
            quals = '~' * len(probe_seq)
            fastq.write('@'+ str(probe_coords) + '\n%s\n+\n%s' % (probe_seq, quals))
            fastq.close()

            #now you want to take the fastq and I think run bt2
            sam_file = tmpdir + '/derived.sam'
            sam_file = bt2_call(fastq_filename,sam_file,200000)

            bam_file = tmpdir + '/derived.bam'
            bam_file = samtools_call(sam_file,bam_file)

            #now make the pairwise file
            pairwise_file = tmpdir + '/derived.out'
            s2p = (sam2pairwise_call(bam_file).decode('utf-8'))

            with open(pairwise_file,'w') as pf:
                for line in s2p:
                    pf.write(line)
                pf.close()

            #write a function to clean the pairwise sequences
            pairwise_df,unique_df = process_pairwise(pairwise_file)

            #for any probes after it, nupack to see if pdups < 0.50
            parent_seq = unique_df["parent"].tolist()
            derived_seq = unique_df["derived"].tolist()

            pdups_vals_list = []
            probe_names_list = []
            for ps,ds in zip(parent_seq,derived_seq):

                #add names and compute pdups to list
                probe_names_list.append(ps + "_" + ds)
                #compute the pdups values
                pdups_val = pdups(ps,ds)
                #add the values into a list
                pdups_vals_list.append(pdups_val)

            #add names and pdups into zipped dict
            probe_pdups_dict = dict(zip(probe_names_list,pdups_vals_list))

            #make a list for all pdups scores
            all_pdups_list = []
            #now do this for the larger df to append vals
            all_parent_seq = pairwise_df["parent"].tolist()
            all_derived_seq = pairwise_df["derived"].tolist()

            #zip the two lists together and append value to key in dict
            for ps,ds in zip(all_parent_seq,all_derived_seq):
                all_dict_key = ps + "_" + ds
                all_pdups_list.append(probe_pdups_dict[all_dict_key])

            pairwise_df['pdups'] = all_pdups_list

            return pairwise_df

###################################################################################

def bt2_call(fq_file,sam_file,k_val):

    subprocess.call(['bowtie2', '-x', bt2_idx,
         '-U', str(fq_file),
         '-k',str(k_val),
         bowtie_string,
         '-S', str(sam_file)])

    return sam_file

###################################################################################

def samtools_call(sam_file,bam_file):

    subprocess.call(["samtools", "view", "-bS", sam_file],
                    stdout=open(bam_file,'w'))

    return bam_file 

###################################################################################

def sam2pairwise_call(bam_file):

    sv = subprocess.Popen(['samtools', 'view', bam_file], \
    stdout=subprocess.PIPE)

    s2p = (subprocess.check_output('sam2pairwise', stdin=sv.stdout))

    return s2p

###################################################################################

def process_pairwise(pairwise_file):

    #create list to add dictionary of rows
    rows_list = []

    with open(pairwise_file) as f:
        for probe_ID, parent, aligns, reference in itertools.zip_longest(*[f]*4):
            #create dict to store the row
            row_dict = {}
            row_dict.update(probe_ID = probe_ID.split('\t')[0],
                            parent = str(parent.replace('-','').strip('N\n')),
                            derived = str(Seq(reference.replace('-','').strip('N\n'))),
                            align_chr = probe_ID.split('\t')[2],
                            align_start = probe_ID.split('\t')[3])

            rows_list.append(row_dict)
        
        pairwise_df = pd.DataFrame(rows_list, columns = ['probe_ID',
                                        'parent',
                                        'derived',
                                        'align_chr',
                                        'align_start'])
    
    #collapse duplicates
    unique_probes_df = pairwise_df.drop_duplicates(subset=['parent',
                                                  'derived'],
                                          keep='first') 

    return pairwise_df, unique_probes_df

###################################################################################

def pdups(seq1,seq2):
    
    a = nupack.Strand(seq1, name='a')
    b = nupack.Strand(rev_comp(seq2), name='b')
    t = nupack.Tube({a: strand_conc_a, b: strand_conc_b}, complexes=nupack.SetSpec(max_size=2), name='t')

    # calculate partition function and resulting concentrations
    tube_result = nupack.tube_analysis(tubes=[t], model=NUPACK_MODEL)
    
    # calculate duplex probability
    conc = tube_result.tubes[t].complex_concentrations[nupack.Complex([a, b])]
    pdups_score = conc / strand_conc_b

    return pdups_score

###################################################################################

def target_nupack_sum(pairwise_df,probe_region_dict):

    #make a col that includes the derived end
    dp_lengths=[]
    derived_probe_list = pairwise_df['derived'].tolist()
    for probe in derived_probe_list:
        dp_lengths.append(len(probe))

    derived_end_list =[]
    pairwise_df.align_start = pairwise_df.align_start.astype(int)

    derived_starts_list = pairwise_df['align_start'].tolist()

    for dpl,ds in zip(dp_lengths,derived_starts_list):
        derived_end_list.append(dpl+ds)

    pairwise_df['align_end'] = derived_end_list

    repeat_coords_list = []
    probe_coords_list = pairwise_df['probe_ID'].tolist()

    for pc in probe_coords_list:
        if pc in probe_region_dict:
            repeat_coords_list.append(probe_region_dict[pc])

    #append repeat coords as a col in df
    pairwise_df['repeat_coords'] = repeat_coords_list

    #format columns with the repeat_coords
    pairwise_df[['parent_chrom','region']] = pairwise_df['repeat_coords'].str.split(':',
                expand=True)
        
    pairwise_df[['r_start','r_end']] = pairwise_df['region'].str.split('-',
                expand=True)
    
    #cast r columns as ints
    pairwise_df[['r_start','r_end']] = pairwise_df[['r_start',
                'r_end']].astype(int)

    #make dicts to store the on_target_sum, off_target, total

    total_prop_dict = {}
    total_on_target_dict = {}
    total_off_target_dict = {}

    #perform a group by intersection on the coords
    grouped_region=pairwise_df.groupby('probe_ID',sort=False)

    for key,item in grouped_region:

        on_total = 0
        off_total = 0

        probe_chrom = key.split(":")[0]

        #probe list
        derived_chrom_list = item['align_chr'].tolist()
        parent_chrom_list = [probe_chrom for i in range(len(derived_chrom_list))]

        #repeat coords list
        repeat_start_list = item['r_start'].tolist()
        repeat_end_list = item['r_end'].tolist()

        #derived probe starts and ends
        derived_start_list = item['align_start'].tolist()
        derived_end_list = item['align_end'].tolist()

        #pdups list
        pdups_list = item['pdups'].tolist()

        on_target_pdups = []
        off_target_pdups = []

        #modify so instead it provides a ratio
        for pc,dc,ds,de,rs,re,pdups in zip(parent_chrom_list,derived_chrom_list,
                               derived_start_list,derived_end_list,
                               repeat_start_list,repeat_end_list,pdups_list):

            if pc == dc and int(ds) >= int(rs) and int(de) <= int(re):
                on_total += pdups
            else:
                off_total += pdups

        #instead of a list this should be a dictionary
        total_prop_dict[key] = (on_total)/(on_total+off_total)
        total_on_target_dict[key] = on_total
        total_off_target_dict[key] = off_total

    return total_prop_dict,total_on_target_dict,total_off_target_dict

###################################################################################

def main():

    probe_df = read_probe_filter(p_file)

    print("---%s seconds ---"%(time.time()-start_time))

    on_target_dict,off_target_dict,prop_target_dict = filter_by_region(probe_df)

    print("---%s seconds ---"%(time.time()-start_time))

    generate_final_df(probe_df,on_target_dict,off_target_dict,prop_target_dict)

if __name__== "__main__":
    main()

