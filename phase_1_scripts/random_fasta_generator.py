"""
Robin Aguilar
Beliveau and Noble Labs
08/15/2019
Purpose: To write a script that will generate a fasta file containing a series of random repetitive elements
of different types. This algorithm should store the locations of the occurence of these instances and should 
containing the following types of repeats: tandem repeats, inverted repeats, interspersed repeats
"""

#import modules
import numpy as np

BASES = ('A', 'C', 'T', 'G')
P = (0.2, 0.2, 0.3,0.3)
repeat_indices_start=[]
repeat_indices_end=[]
rand_indices_start=[]
rand_indices_end=[]
n_indices_start=[]
n_indices_end=[]

def generate_random_fasta(random_length,repeat_span,short_copy_num,long_copy_num):
    
    #make first random seq of given length
    sequence1=''.join(np.random.choice(BASES, p=P) for i in range(random_length))
    
    #make a second random seq
    sequence2=''.join(np.random.choice(BASES, p=P) for i in range(random_length))

    #make a short tandem repeat
    short_tandem_repeat=(''.join((np.random.choice(BASES,p=P)) for i in range(repeat_span)))*short_copy_num
    
    #make a long tandem repeat
    long_tandem_repeat=(''.join((np.random.choice(BASES,p=P)) for i in range(repeat_span)))*long_copy_num
    
    #make an unassembled region
    string_of_n = "N"*50
    
    #cat string together to make fasta
    total_sequence=sequence1+long_tandem_repeat+sequence2+string_of_n+short_tandem_repeat
        
    #append the indices of coords to report where all these start and stop regions happen for each item above
    repeat_indices_start.append(total_sequence.lower().find(str(short_tandem_repeat).lower()))
    repeat_indices_start.append(total_sequence.lower().find(str(long_tandem_repeat).lower()))

    repeat_indices_end.append(total_sequence.lower().find(str(short_tandem_repeat).lower())+len(short_tandem_repeat))
    repeat_indices_end.append(total_sequence.lower().find(str(long_tandem_repeat).lower())+len(long_tandem_repeat))

    rand_indices_start.append(total_sequence.lower().find(str(sequence1).lower()))
    rand_indices_end.append(total_sequence.lower().find(str(sequence1).lower())+len(sequence1))
    
    rand_indices_start.append(total_sequence.lower().find(str(sequence2).lower()))
    rand_indices_end.append(total_sequence.lower().find(str(sequence2).lower())+len(sequence2))
    
    n_indices_start.append(total_sequence.lower().find(str(string_of_n).lower()))
    n_indices_end.append(total_sequence.lower().find(str(string_of_n).lower())+len(string_of_n))


    zip_repeats = zip(repeat_indices_start, repeat_indices_end)
    dict_repeats = dict(zip_repeats)
    repeats=("Indices of tandem repeats: " + str(dict_repeats))
    
    zip_rand=zip(rand_indices_start,rand_indices_end)
    dict_rand=dict(zip_rand)
    random=("Indices of the random dna: " + str(dict_rand))
    
    zip_n=zip(n_indices_start,n_indices_end)
    dict_n=dict(zip_n)
    n=("Indices of n bases: " + str(dict_n))
    
    #writes two files, the fasta and the key files that tells you the indices of all the sequence events
    with open("fasta_reference_15_25_5_10.txt","w") as ref:
        ref.write(repeats + "\n" + random + "\n" + n)
        ref.close()
    
    with open("fasta_seq_15_25_5_10.fa","w") as fasta:
        fasta.write(">randseq"+"\n"+ total_sequence)
        fasta.close()
        
generate_random_fasta(15,25,5,10)
