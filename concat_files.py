"""
Robin Aguilar
Beliveau and Noble labs
Date Modified: 03/01/2020
Modification: To make some functions more efficient and commented with user 
statements 
"""
#import libraries needed
from os import listdir
from os.path import isfile, join
import time

###################################################################################

userInput = argparse.ArgumentParser(description=\
        '%Requires a FASTA file as input. And Jellyfish Index file of genome '
        'single-entry or multi-lined FASTA files are supported.  Returns a dataframe which '
        'can be used for further parsing with probe information.')

requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--file_path', action='store', required=True,
                               help='The filtered probe file that contains all sequences and regions')

requiredNamed.add_argument('-o', '--out_path', action='store', required=True,
                              help='The desired name of the output file')

args = userInput.parse_args()
file_path = args.file_path
out_path=args.out_path

###################################################################################

def list_files(file_path):
    
    #take files in directory and make a list
    files_list = [f for f in listdir(file_path) if isfile(join(file_path, f))]

    return files_list
    
###################################################################################

def concat_files(files_list):
    
    with open(out_path, 'w') as outfile:
        for fname in files_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
                    
###################################################################################
    
def main():
    
    files_list = list_files(file_path)

    print("---%s seconds ---"%(time.time()-start_time))

    concat_files(files_list)

    print("---%s seconds ---"%(time.time()-start_time))

    print("Done")

if __name__== "__main__":
    main()
