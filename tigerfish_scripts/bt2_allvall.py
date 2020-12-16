"""
Robin Aguilar
Beliveau and Noble Labs
Created: 07/29/2020
Purpose: To perform the last filtering step that implements machine learning (LDA)
"""

#import all functions/modules needed
import time
import io
import sys
import re
import csv
import argparse
import subprocess
import numpy as np
import pandas as pd
import tempfile
import fasta2matrix_function as f_matrix
from itertools import product
from itertools import combinations
from itertools import combinations_with_replacement 
import os 

from sklearn.model_selection import StratifiedKFold
from sklearn import svm, datasets
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import classification_report
import pDups as pDuplex

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

###################################################################################
"""
Argparse: The purpose of this is to have all of the arguments documented here
"""
userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--probe_file', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')

requiredNamed.add_argument('-o', '--out_file', action='store', required=True,
                              help='The desired name of the output file')

userInput.add_argument('-t_l', '--threshold_local', action='store', default = 0, type = int,
                               help='parameter to filter in clf predict probas, default=3')

userInput.add_argument('-t_g', '--threshold_global', action='store', default = 0, type = int,
                               help='parameter to filter in clf predict probas, default=3')

###################################################################################
"""
Variable names: These are arguments and variable names that are used throughout the 
script.
"""
args = userInput.parse_args()
p_file = args.probe_file
o_file=args.out_file
threshold_l=args.threshold_local
threshold_g=args.threshold_global

start_time=time.time()

#changed to test if job gets killed
bowtie_string= "-k 1000 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4"

###################################################################################
"""
Data structures. These are structures used to organize sam scores, and probe sets.
"""

sam_scores_list=[]
sequence_list=[]

max_scores_list=[]
seq_names_list=[]


###################################################################################
"""
This function reads the probe file and parses it into a dataframe.
"""

def read_filtered_probes(p_file):

    #label columns
    colnames = ['chr','p_start','p_end','probe','Tm','region','repeat_num','hg_num','score']

    #open the ground probe file and read
    probe_filter = pd.read_csv(p_file, names=colnames, header=None, sep='\t')

    return probe_filter

###################################################################################
"""
This function takes the dataframe, orders it by repeat number, and appends the probes
into a list and a list of their indices. These two lists will be important for
the appending the lookup table aspect of this.
"""

def generate_ordered_df(probe_filter):

    #I wonder if you want to use
    #let's sort in descending order based on the repeat number 
    probe_filter=probe_filter.sort_values(by='repeat_num', ascending=False)

    #remove probe duplicates
    probe_filter = probe_filter.drop_duplicates(subset=['probe'], keep='first')

    #make the probes into a list
    probes_list=probe_filter['probe'].tolist()
    probe_filter['index_col'] = probe_filter.index

    return probe_filter,probes_list

###################################################################################
"""
This function takes a sequence as a string and will generate the RC of that seq.
This is used when feeding in pairs of scripts into the LDA clf. 
"""

def make_seq_RC(sequence):

    probe_seq = Seq(sequence,generic_dna)
    seqRC=str(probe_seq.reverse_complement())

    return seqRC

###################################################################################

"""
This function runs bowtie2 while also generating a lookup table of the best
possible alignment score for each entry, if there is no entry from the sam
file output, then that key gets a value of a 0. 
"""

def all_v_all_lookup_table(scaffold_list,read_id_list,probe_idx_list):

    with tempfile.TemporaryDirectory() as tmpdir:
        #make the list of fasta seqs to build an alignment index
        fasta_file = tmpdir + '/all_v_all.fa'
        with open(fasta_file, 'w') as f:
            for idx,probe in zip(probe_idx_list,scaffold_list):
                f.write('>' + str(idx) + "\n" + str(probe) + "\n")
            f.close()

        #put the sequence to align as a fastq in a tmpdir
        fastq_file = tmpdir + '/all_v_all.fastq'
        quals_list = ['~' * len(probe) for probe in read_id_list]

        fastq_seqs_list=["@"+str(idx)+"\n"+str(probe)+"\n"+"+"+"\n"+str(qual)+"\n" for idx,probe,qual in zip(probe_idx_list,read_id_list,quals_list)]
        with open(fastq_file,'w') as f:
            for item in fastq_seqs_list:
                f.write("%s" % item)

        #build the alignment index
        index_stem = tmpdir + '/derived'
        subprocess.call(['bowtie2-build', fasta_file, index_stem, '-q'],stderr=None, shell=False)

        #align sequences with the index
        sam_file = tmpdir + '/all_v_all.sam'

        subprocess.call(['bowtie2', '-x', index_stem, '-U', fastq_file, bowtie_string, '-S', sam_file, '--no-hd', '--quiet'])

        cols = ['read_id', 'col2', 'scaffold_id', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10', 
            'col11', 'AS', 'col13', 'col14', 'col15', 'col16', 'col17', 'col18', 'col19','col20']

        sam_chunk = pd.read_csv(sam_file, header = None, sep = '\t', names=cols ,chunksize=1000000)

        combos_dict = {}
        for chunk in sam_chunk:
            chunk = chunk.drop(columns = ['col2','col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10','col11','col13','col14', 'col15', 'col16', 'col17', 'col18', 'col19','col20'])
            read_id_list = chunk['read_id'].tolist()
            scaffold_id_list = chunk['scaffold_id'].tolist()
            sam_scores_list = chunk['AS'].tolist()
            sam_scores_list = [(int(string.split(':')[2])) for string in sam_scores_list]
            read_scaffold_key = []
            sam_score_val = []
            for read_id,scaffold_id,sam_scores in zip(read_id_list,scaffold_id_list,sam_scores_list):
                if read_id != scaffold_id:
                    read_scaffold_key.append(str(read_id) + "_" + str(scaffold_id))
                    sam_score_val.append(sam_scores)
            for key,val in zip(read_scaffold_key,sam_score_val):
                if key not in combos_dict:
                    combos_dict[key]=val

        return combos_dict

###################################################################################
"""
This function loads the clf.
"""

def load_clf():

    clf_vals=[[-6.68168328e+00,  2.17578108e+00,  3.15727413e+00, -5.29375005e+00,
   4.32334947e+00, -7.98422848e+00, -5.85814615e+00, -8.48674966e+00,
   6.13222104e+00,  3.53446861e+00, -6.87893257e+00,  1.84971554e+00,
   7.35756937e+00,  1.68303182e+00, -6.58832891e-01, -7.90274741e-01,
  -1.90870942e+00, -4.48892817e+00,  2.02860715e+00,  6.30674636e+00,
   2.35769736e-01,  7.49763000e+00,  4.24633674e+00, -2.23214467e+00,
   4.88854888e+00, -2.19090169e+00,  9.37368877e+00, -6.99511092e+00,
   1.66653139e+00,  1.09171275e+00,  3.59848235e+00, -1.60945783e+01,
   1.26021857e+01,  1.18789071e+01,  7.44128392e+00,  1.18586132e+01,
  -1.22534359e+00,  1.14276004e+00, -1.55429412e+00, -5.52848033e+00,
  -6.69948185e+00, -7.92510100e+00, -2.35094772e-01, -4.53439932e+00,
  -2.57966800e-01, -2.26677043e+00,  4.87012474e-01,  8.44936883e+00,
   4.91941573e+00,  2.37927313e+00, -3.29474246e+00, -5.40560434e+00,
  -5.03575590e+00, -4.73740025e+00, -4.34427594e+00, -1.01251945e+01,
  -5.30328510e+00, -8.16154486e+00, -4.80858452e+00, -4.97403054e+00,
  -2.03617762e+00, -6.45708444e+00, -5.43257952e+00, -3.18516883e+00,
   9.01199755e+00,  3.60309688e+00, -9.63275217e-01, -6.65560079e+00,
  -1.76675059e+01,  7.43902757e+00,  2.61895115e+00, -2.23878931e+01,
   8.16155889e+00,  4.75678855e+00,  2.13174201e+01,  6.59311597e-01,
   6.07077797e+00,  7.75152196e+00,  1.04978174e+01, -6.83486380e+00,
   7.73574320e+00, -5.95086239e+00,  3.00479217e-01,  7.30205284e+00,
  -1.05755265e+00, -3.62505706e+00, -2.63345722e+00, -3.53787174e+00,
  -5.91659769e+00, -2.79449413e+00, -6.36546431e+00, -6.48780248e+00,
  -4.48652605e+00, -3.80988585e+00, -4.24337729e+00, -5.34102501e+00,
   2.44046584e+00,  2.31804019e+00,  1.44142353e+00,  2.27805081e+00,
  -6.05864944e-01,  1.09392138e+00,  2.37745318e+00,  8.31072154e-01,
   4.17912982e+00,  5.98858128e+00,  6.46676615e+00,  4.89531966e+00,
   1.06936566e+00, -7.45143557e-01, -1.39229532e+00,  8.58445591e-01,
   1.11085284e+01,  1.20724843e+01,  1.23464193e+01,  1.18185656e+01,
   3.33868194e+00,  3.69178652e+00,  4.92686565e+00,  2.35221126e+00,
   2.51380985e+00,  1.78592619e+00,  6.18056234e+00,  3.05313041e+00,
   3.47711071e-01,  1.66634647e+00,  1.34165371e+00,  1.35726201e+00,
   1.40968318e+01,  1.88927124e+01,  1.93596798e+01,  2.05189852e+01,
   1.43247804e+01, -3.04261281e+00, -5.06434155e-01,  1.09361064e+01,
  -1.01930018e+01, -7.44163137e+00, -1.47613940e+01, -1.21095733e+00,
  -7.84692086e+00, -9.98461057e+00, -1.07403460e+01,  3.23166546e+00,
  -7.23687353e+00,  2.34592416e+00, -6.38232027e+00, -7.90189555e+00,
   4.76884304e-01, -2.12618661e+00, -9.66981522e-01, -1.39607384e+00,
  -5.45523615e+00, -2.70671165e+00, -6.58750684e+00, -5.92465317e+00,
  -4.40561850e+00, -3.83635020e+00, -5.03767972e+00, -5.27635053e+00,
   1.26258708e+00,  9.74717141e-01,  8.87886213e-01,  1.17350586e+00,
  -1.13663777e-01,  1.89877846e+00,  3.59548636e+00,  1.48968553e+00,
   1.96873590e+00,  2.70914456e+00,  2.07308461e+00,  3.56673361e+00,
  -2.22049920e+00, -2.25808925e+00, -3.30828724e+00, -2.63804668e+00,
  -5.00522763e-01, -5.60494585e-01, -3.88084752e-01, -1.07888136e+00,
   2.02556931e+00,  5.14197741e+00,  4.51660597e+00,  1.47637151e+00,
   3.97367111e+00,  2.45646845e+00,  6.99672412e+00,  3.91494617e+00,
   9.89578621e-01,  3.01580510e+00,  3.11067776e+00,  3.33593289e+00,
  -5.55071713e+00, -5.12259254e+00, -4.55075507e+00, -4.03817762e+00,
   1.23275051e+01, -5.75357479e+00, -1.22740335e+00,  6.97701377e+00,
  -9.86521577e+00, -7.14051044e+00, -1.27258365e+01, -7.35047042e-01,
  -7.39078572e+00, -1.06625531e+01, -9.66450814e+00,  2.74845723e+00,
  -2.95371594e-01,  9.16885823e+00, -7.49243873e-01, -1.49917755e+00,
   6.29733087e-01, -3.28029935e+00, -9.38063482e-01, -2.42048645e+00,
  -2.09537243e+00, -2.72675748e-01, -3.68093321e+00, -3.73718880e+00,
  -3.92371446e+00, -2.06822451e+00, -3.15281457e+00, -4.81411795e+00,
   2.60362841e+00,  4.36999587e+00,  3.65681812e+00,  3.35238971e+00,
  -1.30999127e-01,  1.41794706e+00,  2.56554165e+00,  1.02419461e+00,
   4.32238434e+00,  5.64454230e+00,  4.84064480e+00,  5.47280834e+00,
  -1.37006809e+00,  1.33369002e-01, -1.27567320e+00, -1.03923802e+00,
   7.75979511e-01,  4.37434410e-01,  9.32806776e-01,  6.94006017e-01,
   6.81176699e-01,  2.31794080e+00,  3.27230796e+00,  8.58465370e-01,
   3.98465639e+00,  2.57136611e+00,  7.63257445e+00,  4.70361024e+00,
   1.63175717e+00,  2.94263615e+00,  2.88151801e+00,  3.52941851e+00,
  -2.87842047e+00,  1.28286557e+00,  1.35485556e+00,  3.87788783e+00,
   1.26668421e+01, -4.47664421e+00, -3.94750505e-01,  7.58887663e+00,
  -9.93097552e+00, -7.75428673e+00, -1.46264476e+01, -7.98796651e-01,
  -9.15520232e+00, -1.14977090e+01, -1.09460421e+01,  1.67693209e+00,
  -2.80457996e-01,  8.81273597e+00,  8.08417062e-01, -2.76084945e+00,
  -3.98043822e-01, -3.43484085e+00, -2.99881284e-01, -2.80045689e+00,
  -4.86296221e+00, -7.55345806e-01, -6.41202629e+00, -7.04184611e+00,
  -4.92614424e+00, -2.63682501e+00, -4.95322022e+00, -5.04215742e+00,
   8.12042127e+00,  8.66969910e+00,  7.49509999e+00,  5.79792498e+00,
  -2.38847183e-01,  1.95128761e+00,  4.16327478e+00,  2.19199611e+00,
   2.82040454e+00,  3.29697138e+00,  1.93055586e+00,  2.87976362e+00,
  -2.99718243e+00, -1.75036211e+00, -5.00310302e+00, -2.60106714e+00,
   1.31685173e-02, -1.16622651e+00,  6.50601597e-02, -1.42832094e+00,
   2.04455962e+00,  3.66761268e+00,  4.20157676e+00,  1.67144421e+00,
   3.46318049e+00,  3.28047228e+00,  7.37695026e+00,  5.07921611e+00,
   4.51220202e-01,  1.39139790e+00,  1.95910887e+00,  1.31551340e+00,
  -9.91661138e-01,  2.54115734e+00,  2.94402566e+00,  4.11130749e+00,
   1.33468486e+01, -4.63846505e+00, -2.20065482e+00,  7.65047453e+00,
  -9.94945912e+00, -6.72137015e+00, -1.39258227e+01, -2.30355577e-01,
  -7.89810756e+00, -1.00514773e+01, -1.07544237e+01,  3.22180383e+00,
  -7.64044663e+00,  1.76108487e+00, -7.08231368e+00, -7.83292484e+00,
   3.45441238e+00, -5.73915089e+00, -4.85108570e+00,  2.67519300e+00,
  -9.68266216e-02, -1.21710622e+00, -1.30613248e+00, -3.45054751e+00,
  -4.76642898e-01,  7.69861388e-01,  1.93194527e+00, -2.02330021e+00,
   1.70784123e+00, -3.48534053e-01, -1.51822149e+00, -2.91794675e+00,
   8.85593865e+00,  4.50671885e+00,  1.96176135e+00,  1.28800781e-01,
  -3.43592685e+00,  1.10713113e+00, -3.96851621e-01, -6.82871881e-02,
  -2.45922962e-02,  2.17384933e+00,  1.66208121e+00,  2.18016553e+00,
  -1.57433215e+00,  2.97584059e+00,  2.61019590e+00,  2.39461570e+00,
  -7.97354870e+00, -1.39267023e+00, -1.67234504e-01,  2.68197032e-01,
  -8.02304344e-02,  3.38361107e+00,  2.29703177e+00,  2.31332536e+00,
   1.05588227e+00,  4.40306013e+00,  2.76322474e+00,  2.83554956e+00,
   8.35413227e-02,  5.02833250e+00,  4.94369948e+00,  3.16729116e+00,
  -5.25564253e+00,  1.05031155e+00,  1.61134925e+00,  1.47769087e+00,
  -1.74622579e+00,  1.84008149e+00,  7.58273912e-01,  3.41302233e-01,
   2.22334347e+00,  4.74592081e+00,  5.36449916e+00,  3.44124759e+00,
   1.02596232e+00,  5.87414501e+00,  5.40984700e+00,  4.65472943e+00,
  -5.24175019e+00,  1.18298593e+00,  1.22230272e+00,  2.25083334e+00,
  -8.64935566e+00, -5.77320770e+00, -7.31377877e+00, -6.88650199e+00,
  -3.08263675e+00, -1.62629837e+00, -3.16098028e+00, -1.99806810e+00,
  -3.22356022e+00,  6.74233523e-01,  7.72659201e-02, -6.69080300e-01,
  -8.90870045e+00, -3.01912899e+00, -2.50261758e+00, -2.74290469e+00,
   8.55497035e-01,  2.06779665e+00,  2.31409264e+00,  1.55088006e+00,
  -1.58927264e+00, -3.57273562e-01, -1.41151742e-02, -1.63469043e+00,
  -7.11240289e-01,  4.76862599e-01,  1.31533995e+00, -2.47600010e-01,
  -7.24460499e-01,  1.42523918e-01,  5.93447838e-02, -1.14966811e+00,
  -1.51027731e+00,  1.00634323e+00,  1.04953893e+00,  2.85056392e-02,
  -6.55710854e-01, -4.25891512e-01, -2.38589870e-01, -4.11703174e-01,
  -1.23656576e+00,  1.45748693e+00,  6.49022509e-01, -3.58529296e-02,
  -2.95311357e+00, -1.25422744e+00, -3.19362750e-01, -1.56072542e+00,
  -2.15053274e-02,  1.54233133e+00,  2.23028917e+00,  1.26129395e+00,
  -1.90845615e+00, -5.62485111e-01, -2.29876361e+00, -8.74335550e-01,
  -2.01244780e+00,  1.15734881e-01, -9.16170049e-02, -1.07807268e+00,
  -3.55324385e+00, -4.04547131e-01, -2.64999328e-01, -2.13896024e+00,
   7.17772162e+00,  8.73758851e+00,  9.07167618e+00,  6.84910157e+00,
   1.77238224e+00,  3.25817933e+00,  3.90682423e+00,  2.53623302e+00,
   6.04774290e-01,  2.78905373e+00,  2.29636542e+00,  1.80357261e+00,
  -4.28162306e-01,  1.82832984e+00,  2.15769235e+00,  8.02882648e-01,
   2.41041526e+00,  2.75360026e+00,  3.98690718e+00,  2.35491001e+00,
   1.40933154e-03,  1.50814654e+00,  2.07546329e+00,  3.06093741e-01,
   5.13268271e-01,  2.68321772e+00,  2.34978886e+00,  9.65712221e-01,
  -3.97219352e-01,  1.60627945e+00,  2.06459611e+00,  4.87766764e-02,
   5.55511238e-01,  1.21875962e+00,  1.84097757e+00,  1.05917458e+00,
  -5.29755802e-01,  1.04254361e+00,  8.73748026e-02,  7.39986447e-02,
  -4.84929367e-02,  6.20192964e-01,  1.07095392e+00,  2.04624940e+00,
  -1.70742437e+00,  4.64905000e-01,  1.04074606e+00, -6.33604046e-01,
  -2.17318521e-01,  2.24450779e+00,  1.88334187e+00,  1.11097568e+00,
  -2.00684296e+00, -1.49814973e+00, -8.67867704e-01, -5.77894721e-01,
  -1.48908126e+00,  2.78059096e-01,  6.25957824e-02, -2.02944706e+00,
  -3.17988134e+00, -5.15314737e-01,  9.75115349e-01, -1.28491403e+00,
   7.28456300e+00,  8.37595737e+00,  9.51828797e+00,  7.77902438e+00,
   2.32644149e+00,  3.91076965e+00,  3.70597367e+00,  2.92268072e+00,
   1.22159705e+00,  3.98483369e+00,  4.18660563e+00,  2.77823623e+00,
   2.24096502e-01,  2.70485912e+00,  3.58432866e+00,  2.10322350e+00,
   2.52623556e+00,  2.56262026e+00,  3.16300497e+00,  2.32306980e+00,
  -1.32128840e-01,  9.82419672e-01,  7.70491672e-01,  4.78780028e-01,
   2.36803339e-01,  1.91129374e+00,  2.32800262e+00,  1.00622396e-01,
  -1.00493557e+00,  2.07821275e+00,  2.20231870e+00, -1.08147107e-02,
  -1.77431256e-01,  1.88821015e+00,  2.23224370e+00,  7.37039990e-01,
   4.08398684e-01,  1.02863671e+00,  1.14897738e+00,  8.06065159e-02,
  -1.95860538e+00,  2.89118681e-01,  6.00709871e-01,  1.15715242e+00,
  -1.88862972e+00,  1.96142642e-01,  1.28106450e+00, -1.75679780e-01,
   1.38907548e+00,  2.94089100e+00,  3.24979018e+00,  1.95575241e+00,
  -5.26973962e-01,  5.20860702e-01, -1.09008219e+00,  1.37067677e-01,
  -9.44781298e-01,  2.58477666e-01,  1.07186564e+00, -5.75390725e-01,
  -2.50015991e+00,  3.59779273e-01,  1.60501544e+00, -8.74375440e-01,
   6.89366676e+00,  8.89274791e+00,  9.91894571e+00,  8.08845313e+00,
   1.89016322e+00,  4.09037523e+00,  2.97015009e+00,  2.88159584e+00,
   1.60789961e+00,  3.82856661e+00,  4.61878828e+00,  3.36449199e+00,
   6.58331157e-01,  2.92818280e+00,  2.85782503e+00,  1.49102433e+00,
  -4.10324452e-01, -4.49756142e-01,  8.27962868e-01, -1.37727810e-02,
  -2.17614497e+00, -7.17692557e-01, -1.85260442e+00, -2.50276858e+00,
  -2.14863436e+00, -8.83294567e-01, -6.86990572e-01, -1.37406618e+00,
  -2.26717583e+00,  1.42317493e-02, -5.93583103e-01, -1.83609058e+00,
  -2.06878687e+00, -7.55058505e-01, -8.77824492e-01, -1.25303069e+00,
  -1.89250824e+00, -1.00445350e+00, -2.20038755e-01, -9.17436003e-01,
  -1.42028636e+00,  7.66247744e-01,  3.90367219e-01, -7.06191195e-01,
  -3.46936442e+00, -9.15477062e-01, -1.19713489e+00, -2.33625757e+00,
  -5.18007335e-01,  6.10090624e-01,  1.26570498e+00,  4.69246991e-01,
  -2.39962571e+00, -1.00801453e+00, -2.62149901e+00, -1.51556496e+00,
  -2.50050920e+00, -7.88488989e-01, -2.46897670e-01, -1.36617013e+00,
  -3.40124514e+00, -9.98990539e-01,  2.28833392e-01, -2.32857553e+00,
   5.31896098e+00,  6.83549362e+00,  7.63098493e+00,  6.10937360e+00,
   1.04857204e+00,  1.93931898e+00,  2.32158621e+00,  9.37240036e-01,
   1.85849476e-01,  1.66141804e+00,  2.24777049e+00,  7.23567115e-01,
  -6.54038605e-02,  1.44023133e+00,  1.94489494e+00,  3.82818364e-01,
   2.27494108e-01]]

    clf_intercept = [-8.6334261]
    class_list=[-1,1]

    coef_array=np.asarray(clf_vals)
    inter_array=np.asarray(clf_intercept)
    class_array=np.asarray(class_list)

    clf=LinearDiscriminantAnalysis()

    clf.coef_=coef_array
    clf.intercept_=inter_array
    clf.classes_=class_array

    return clf

###################################################################################

def repeat_region_filter(probe_filter,clf):

    probes_all = probe_filter['probe'].tolist()
    index_all = probe_filter['index_col'].tolist()    
    scores_dict = all_v_all_lookup_table(probes_all,probes_all,index_all)

    grouped_probes=probe_filter.groupby('region')

    idx_df_dict = {}

    rm_set = set()

    #iterate over the groups
    for name,group in grouped_probes:

        group_probe_list=group['probe'].tolist()
        group_idx_list=group['index_col'].tolist()

        #scores_dict = all_v_all_lookup_table(group_probe_list,group_probe_list,group_idx_list)
           
        while len(group_probe_list) and len(group_idx_list) > 1:

           idx_key = []
           seq_RC_list = []
           seq_AS = []

           #set up head_probe_list
           head_probe = group_probe_list[0]
           head_probe_l = [head_probe for i in range(1,len(group_probe_list))]
            
           #set up head probe idx
           head_probe_idx = group_idx_list[0]
           head_probe_idx_l = [head_probe_idx for idx in range(1,len(group_idx_list))]
            
           #set up compare probe list
           compare_probe_l = group_probe_list[1:]

           #return seq_RC for each probe to compare
           seq_RC_list = [make_seq_RC(cp) for cp in compare_probe_l]

           #set up the idx for compare idx
           compare_probe_idx_l = group_idx_list[1:]
            
           #make the seq_AS list
           idx_key = [str(hpi) + "_" + str(cpi) for hpi,cpi in zip(head_probe_idx_l, compare_probe_idx_l)]

           #add the seq_AS score
           for idx in idx_key:
               if idx in scores_dict:
                   seq_AS.append(scores_dict[idx])
               else:
                   seq_AS.append(0)

           #make this loop go through so that anything with a 0 can be removed
           for hp,cp,cp_rc,cpi,idxs,as_score in zip(head_probe_l,
                                                            compare_probe_l,
                                                            seq_RC_list,
                                                            compare_probe_idx_l,
                                                            idx_key,seq_AS):
            
               if as_score != 0:
                
                   df_hp_list = []
                   df_cp_list = []
                   df_cp_rc_list = []
                   df_idx_key_list = []
                   df_cpi_list = []
                   df_seq_AS = []
                
                   df_hp_list.append(hp)
                   df_cp_list.append(cp)
                   df_cp_rc_list.append(cp_rc)
                   df_cpi_list.append(cpi)
                   df_seq_AS.append(as_score)       
                   df_idx_key_list.append(idxs)

                   seq_AS_l = [float(i) for i in df_seq_AS]
                   seq_as_arr = np.asarray(seq_AS_l)
                   seq_as_arr = seq_as_arr.reshape(seq_as_arr.shape[0],-1)

                   head_probe_arr = f_matrix.fasta2matrix(4,df_hp_list)
                   compare_probes_arr = f_matrix.fasta2matrix(4,df_cp_rc_list)
            
                   head_pair_concat = np.hstack((head_probe_arr,compare_probes_arr))
                   head_pair_np = np.append(head_pair_concat,seq_as_arr,1)
                   df_result = clf.decision_function(head_pair_np)
                   df_result_list = df_result.tolist()

                   idx_key_dict = dict(zip(df_idx_key_list,df_result_list))
                   idx_df_dict.update(idx_key_dict)
                       
                   #probe_df = pd.DataFrame(list(zip(df_hpi_list,
                   #                df_cpi_list, df_result_list)), 
                   #columns =['hpi','cpi','df']) 

                   for cpi,cp,df in zip(df_cpi_list,df_cp_list,
                                            df_result_list):
                       if df >= float(threshold_l):
                           rm_set.add(cpi)
                           group_idx_list.remove(cpi)
                           group_probe_list.remove(cp)
                       

           group_probe_list.pop(0)
           group_idx_list.pop(0)
               
    #conver sets to list
    rm_list = list(rm_set)

    idx_df_dict = {key:val for key, val in idx_df_dict.items() if val < threshold_l}

    #in probe filter remove the probes to remove
    probe_filter = probe_filter.loc[~probe_filter['index_col'].isin(rm_list)]

    all_probes = probe_filter['probe'].tolist()
    all_idx = probe_filter['index_col'].tolist()
    
    #print(len(probe_filter))

    #print(probe_filter)

    return probe_filter,all_probes,all_idx,idx_df_dict,scores_dict

###################################################################################

def filter_all(probe_filter,clf,all_probes,all_idx,idx_df_dict,scores_dict):

    rm_set = set()

    #scores_dict = all_v_all_lookup_table(all_probes,all_probes,all_idx)
    
    while len(all_probes) and len(all_idx) > 1:

        idx_key = []
        seq_RC_list = []
        seq_AS = []
        
        #set up head_probe_list
        head_probe = all_probes[0]
        head_probe_l = [head_probe for i in range(1,len(all_probes))]
            
        #set up head probe idx
        head_probe_idx = all_idx[0]
        head_probe_idx_l = [head_probe_idx for idx in range(1,len(all_idx))]
            
        #set up compare probe list
        compare_probe_l = all_probes[1:]
        seq_RC_list = [make_seq_RC(cp) for cp in compare_probe_l]
        
        compare_probe_idx_l = all_idx[1:]
            
        #make the seq_AS list
        idx_key = [str(hpi) + "_" + str(cpi) for hpi,cpi in zip(head_probe_idx_l,
                   compare_probe_idx_l)]
           
        for idx in idx_key:
            if idx in scores_dict:
                seq_AS.append(scores_dict[idx])
            else:
                seq_AS.append(0)
                
        #make this loop go through so that anything with a 0 can be removed
        for hp,cp,cp_rc,cpi,idxs,as_score in zip(head_probe_l,compare_probe_l,seq_RC_list,
                                         compare_probe_idx_l,idx_key,
                                         seq_AS):
            
            if as_score != 0 and idxs not in idx_df_dict:
                
                df_hp_list = []
                df_cp_list = []
                df_cp_rc_list = []
                df_cpi_list = []
                df_seq_AS = []
                
                df_hp_list.append(hp)
                df_cp_list.append(cp)
                df_cp_rc_list.append(cp_rc)
                df_cpi_list.append(cpi)
                df_seq_AS.append(as_score)
                
                
                #make lists to use of each item
                seq_AS_l = [float(i) for i in df_seq_AS]
                seq_as_arr = np.asarray(seq_AS_l)
                seq_as_arr = seq_as_arr.reshape(seq_as_arr.shape[0],-1)

                head_probe_arr = f_matrix.fasta2matrix(4,df_hp_list)
                compare_probes_arr = f_matrix.fasta2matrix(4,df_cp_rc_list)
            
                head_pair_concat = np.hstack((head_probe_arr,compare_probes_arr))
                head_pair_np = np.append(head_pair_concat,seq_as_arr,1)
            
                df_result = clf.decision_function(head_pair_np)
                df_result_list = df_result.tolist()           

                #probe_df = pd.DataFrame(list(zip(df_hpi_list,
                #      df_cpi_list, df_result_list)), 
                #       columns =['hpi','cpi','df']) 

                for cp, cpi, df in zip(df_cp_list,df_cpi_list,
                                    df_result_list):
                    if df >= float(threshold_g):
                        rm_set.add(cpi)
                        all_idx.remove(cpi)
                        all_probes.remove(cp)
            
        all_probes.pop(0)
        all_idx.pop(0)
    
    #conver sets to list
    rm_list = list(rm_set)
    #print(rm_list)
    #in probe filter remove the probes to remove
    probe_filter = probe_filter.loc[~probe_filter['index_col'].isin(rm_list)]

    #print(probe_filter)

    return probe_filter

###################################################################################

def write_probe_sets(probe_filter):

    probe_filter.to_csv(str(o_file), sep='\t',index=False,header=False)

###################################################################################

def main():

    probe_filter = read_filtered_probes(p_file)

    print("---%s seconds ---"%(time.time()-start_time))

    probe_filter,probes_list = generate_ordered_df(probe_filter)

    print("---%s seconds ---"%(time.time()-start_time))

    clf = load_clf()

    print("---%s seconds ---"%(time.time()-start_time))

    probe_filter,all_probes,all_idx,idx_df_dict,scores_dict = repeat_region_filter(probe_filter,clf)

    print("---%s seconds ---"%(time.time()-start_time))
    
    probe_filter = filter_all(probe_filter,clf,all_probes,all_idx,idx_df_dict,scores_dict)

    print("---%s seconds ---"%(time.time()-start_time))

    write_probe_sets(probe_filter)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
               
