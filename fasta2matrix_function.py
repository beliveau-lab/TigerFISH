#!/usr/bin/env python
# FILE: fasta2matrix.py
# AUTHOR: William Stafford Noble
# CREATE DATE: 27 September 2005
# MODIFIED BY: Robin Aguilar
# EDIT DATE: 17 December 2019

import sys
import math
import numpy as np
from numpy import array
import pandas as pd
import argparse
import bt2_allvall as pdups_model
#import model_performance as pdups_model

#############################################################################
def make_kmer_list (k, alphabet):

    # Base case.
    if (k == 1):
        return(alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return([])

    # Error case.
    if (k < 1):
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in make_kmer_list(k-1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])

    return(return_value)

##############################################################################
def make_upto_kmer_list (k_values,
                         alphabet):

    # Compute the k-mer for each value of k.
    return_value = []
    for k in k_values:
        return_value.extend(make_kmer_list(k, alphabet))

    #print("Values of k_values: " + " " + str(k_values))

    return(return_value)

##############################################################################
def normalize_vector (normalize_method,
                      k_values,
                      vector,
                      kmer_list):

    # Do nothing if there's no normalization.
    if (normalize_method == "none"):
        return(vector)

    # Initialize all vector lengths to zeroes.
    vector_lengths = {}
    for k in k_values:
        vector_lengths[k] = 0

    # Compute sum or sum-of-squares separately for each k.
    num_kmers = len(kmer_list)
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        if (normalize_method == "frequency"):
            vector_lengths[kmer_length] += count
        elif (normalize_method == "unitsphere"):
            vector_lengths[kmer_length] += count * count

    # Square root each sum if we're doing 2-norm.
    if (normalize_method == "unitsphere"):
        for k in k_values:
            vector_lengths[k] = math.sqrt(vector_lengths[k])

    # Divide through by each sum.
    return_value = []
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        vector_length = vector_lengths[kmer_length]
        if (vector_length == 0):
            return_value.append(0)
        else:
            return_value.append(float(count) / float(vector_length))

    return(return_value)

##############################################################################
def make_sequence_vector (sequence,
                          revcomp,
                          revcomp_dictionary,
                          normalize_method,
                          k_values,
                          mismatch,
                          alphabet,
                          kmer_list,
                          pseudocount):

    # Make an empty counts vector.
    kmer_counts = {}

    # Iterate along the sequence.
    for k in k_values:
        seq_length = len(sequence) - k + 1
        for i_seq in range(0, seq_length):
            kmer=sequence[i_seq : i_seq + k]

            # If we're doing reverse complement, store the count in the
            # the version that starts with A or C.
            if (revcomp == True):
                rev_kmer = find_revcomp(kmer, revcomp_dictionary)
                if (cmp(kmer, rev_kmer) > 0):
                    kmer = rev_kmer

            # Increment the count.
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1

            # Should we also do the mismatch?
            if (mismatch != False):

                # Loop through all possible mismatches.
                for i_kmer in range(0, k):
                    for letter in alphabet:

                        # Don't count yourself as a mismatch.
                        if (kmer[i_kmer:i_kmer+1] != letter):

                            # Find the neighboring sequence.
                            neighbor = substitute(i_kmer, letter, kmer)

                            # If we're doing reverse complement, store the
                            # count in the version that starts with A or C.
                            if (revcomp == 1):
                                rev_kmer = find_revcomp(kmer,
                                                        revcomp_dictionary)
                                if (cmp(kmer, rev_kmer) > 0):
                                    kmer = rev_kmer


                           # Increment the neighboring sequence.
                            if (kmer_counts.has_key(neighbor)):
                                kmer_counts[neighbor] += mismatch
                            else:
                                kmer_counts[neighbor] = mismatch

    # Build the sequence vector.
    sequence_vector = []
    for kmer in kmer_list:
        if kmer in kmer_counts:
            sequence_vector.append(kmer_counts[kmer] + pseudocount)
        else:
            sequence_vector.append(pseudocount)

    # Normalize it
    return_value = normalize_vector(normalize_method,
                                    k_values,
                                    sequence_vector,
                                    kmer_list)

    return(return_value)

##############################################################################

def fasta2matrix(k,my_sequences):

    k_val=k
    sequences=my_sequences
    upto = k_val
    revcomp = False
    normalize_method = "unitsphere"
    alphabet = "ACGT"
    mismatch=False
    pseudocount=0

    # Check for reverse complementing non-DNA alphabets.
    if ((revcomp == True) and (alphabet != "ACGT")):
      sys.stderr.write("Attempted to reverse complement ")
      sys.stderr.write("a non-DNA alphabet (%s)\n" % alphabet)

    # Make a list of all values of k.
    k_values=list(range(k_val,-1,-1))
    k_values.sort() 
    # Make a list of all k-mers.
    kmer_list = make_upto_kmer_list(k_values, alphabet)

    # Set up a dictionary to cache reverse complements.
    revcomp_dictionary = {}

    # Use lexicographically first version of {kmer, revcomp(kmer)}.
    if (revcomp == True):
      new_kmer_list = []
      for kmer in kmer_list:
          rev_kmer = find_revcomp(kmer, revcomp_dictionary)
          if (cmp(kmer, rev_kmer) <= 0):
              new_kmer_list.append(kmer)
      kmer_list = new_kmer_list;

    vector_list=[]
    #Loop through each sequence in list
    for sequence in my_sequences:
        #generate a vector for each sequence
        vector = make_sequence_vector(sequence,
                                revcomp,
                                revcomp_dictionary,
                                normalize_method,
                                k_values,
                                mismatch,
                                alphabet,
                                kmer_list,
                                pseudocount)
        #append vector item into list
        vector_list.append(vector)

    #append total list into dataframe
    vector_dataframe=pd.DataFrame(vector_list)

    #column headers are the corresponding kmers in lexico. order
    vector_dataframe.columns=kmer_list

    return(vector_dataframe)


##############################################################################
def find_revcomp (sequence,
                  revcomp_dictionary):

    # Save time by storing reverse complements in a hash.
    if sequence in revcomp_dictionary:

        return(revcomp_dictionary[sequence])

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if (letter == "A"):
            return_value = return_value + "T"
        elif (letter == "C"):
            return_value = return_value + "G"
        elif (letter == "G"):
            return_value = return_value + "C"
        elif (letter == "T"):
            return_value = return_value + "A"
        elif (letter == "N"):
            return_value = return_value + "N"
        else:
            sys.stderr.write("Unknown DNA character (%s)\n" % letter)
            sys.exit(1)

    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value

    return(return_value)


##############################################################################
# MAIN
##############################################################################

"""
Usage: This script is implemented as a function in pDups_model_clean.py
The primary function - fasta2matrix, will take a list of DNA sequences
and will populate a matrix with kmer count frequency. Users can specify the
value of k which is an argument in the pDups_model_clean.py script. Here, the
script is called as a pDups object.The output is a dataframe of the matrix.

The default parameters include unitsphere normalization and reverse complement
consideration.
"""

def main():
    
    #Define variables and default options
    k_val=pdups_model.k_val
    sequences=pdups_model.seq_list
    upto = k_values

    #Generate vector of specified k and list of sequences
    vector_dataframe=fasta2matrix(k_val,sequences)

if __name__== "__main__":
    main()

