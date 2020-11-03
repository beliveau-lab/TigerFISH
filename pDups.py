import os
import argparse
import random
import csv
import timeit
import time
import tempfile

import subprocess
from subprocess import Popen, PIPE

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def calculateDup(seq, seqRC, correctedTemp, sodiumConc):
    # Generate random number for parallelization
    randomInt = random.randrange(100000000)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write to .in file
        in_File = open('%s/strand_%d_temp.in' % (tmpdir, randomInt), 'w')
        in_File.write('2\n%s\n%s\n2\n' % (seq, seqRC))
        in_File.close()
    
        # Write to .con file
        con_File = open('%s/strand_%d_temp.con' % (tmpdir, randomInt), 'w')
        con_File.write('1E-6\n1E-12\n')
        con_File.close()

        # Run "complexes"
        subprocess.call(['complexes', '-T', correctedTemp, '-sodium', sodiumConc,
                        '-material', 'dna1998', '-quiet', '%s/strand_%d_temp'
                        % (tmpdir, randomInt)])

        # Run "concentrations"
        subprocess.call(['concentrations', '-cutoff', '1E-20', '-sort', '0', '-quiet', '%s/strand_%d_temp'
                        % (tmpdir, randomInt)])
    
        # Extract concentration value
        eq_File = open('%s/strand_%d_temp.eq' % (tmpdir, randomInt), 'r')
        for line in eq_File:
            components = line[:-1].split('\t')
            if len(components) > 4 and components[1] == '1' and components[2] == '1':
                conc_val = components[4]

        # Return concentration value
        if conc_val is not None:
            return conc_val

def main():
    startTime = timeit.default_timer()

    userInput = argparse.ArgumentParser()
    userInput.add_argument('-inputFile', '--inputFile', action='store', type=str,
                            help='file of randomly generated sequences, txt format')
    userInput.add_argument('-T', '--T', action='store', type=float,
                            help='valid temperature, without formamide correction')
    userInput.add_argument('-F', '--F', action='store', type=float,
                            help='float between 0 and 100')
    userInput.add_argument('-S', '--S', action='store', type=float,
                            help='valid sodium concentrtion (molar)')
    userInput.add_argument('-outputFile', '--outputFile', action='store', type=str,
                            help='file suffix')

    args = userInput.parse_args()
    inputFile = args.inputFile
    temperature = args.T
    percentFormamide = args.F
    sodiumConc = args.S
    outputFile = args.outputFile

    # corrected temperature given formamide
    correctedTemp = temperature + (.65 * percentFormamide)

    # Create "temp" directory
    # Cluster path: '/net/beliveau/vol1/home/chernr/biotemp/'
    path = '.'
    if not os.path.exists(path):
        os.mkdir(path)

    seqList = []
    pDuplexList = []

    num = 0
    with open(inputFile, 'r') as f:
        for line in f:
                # Convert str to seq using BioSeq
                seq = Seq(line.strip(), generic_dna)
                seqRC = seq.reverse_complement()

                # Calculate pDuplex
                conc_val = calculateDup(str(seq), str(seqRC), str(correctedTemp), str(sodiumConc))
                pDuplex = ((float(conc_val))/(float('1.0e-12')))

                # Add result
                seqList.append(str(seq))
                pDuplexList.append(str(pDuplex))
    
    # Creates file if it does not already exist
    if not os.path.exists('%s' % outputFile):
        subprocess.call(['touch', '%s' % outputFile])
        header = open('%s' % outputFile, 'w')
        header.write('seq\t')
        header.write('pDuplex\t')
        header.write('conditions\n')
        header.close()
    
    # Write to file
    out_File = open('%s' % outputFile, 'a')

    conditions = 'T%s_F%s_S%s' % (temperature, percentFormamide, sodiumConc)

    for i in range(len(seqList)):
        out_File.write(seqList[i])
        out_File.write('\t')
        out_File.write(pDuplexList[i])
        out_File.write('\t')
        out_File.write(conditions)
        out_File.write('\n')

    out_File.close()

    # Print start --> end time
    endTime = timeit.default_timer()
    print('Program took %f seconds' % (endTime - startTime))

if __name__ == '__main__':
    main()

