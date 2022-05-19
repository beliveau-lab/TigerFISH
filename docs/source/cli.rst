Command Line Interface
######################

`Tigerfish` was written by Robin Aguilar while in the Department of Genome Sciences at the University of Washington.

Official code website:

More documentation and examples: 

The `Tigerfish` workflow is managed by Snakemake which is used to scale and automatically deploy pipeline jobs in parallel as described in our getting started page.

Essentially, each script in the Tigerfish workflow has a set of parameters that may be modified by users. Each script in the main workflow implements these parameters which are called from the pipeline's config.yml file. As a user, one would only need to modify arguments from the config.yml file in order to get Tigerfish working. But to be thorough, named arguments for each script are defined in detail below. 

Main workflow
#############

Here is a simplified DAG of how Snakemake implements scripts in the Tigerfish workflow:

.. image:: imgs/snakemake.svg
   :width: 400
   :alt: A picture of a DAG of the Tigerfish workflow


The scripts below are presented in the order that they are executed by the `Tigerfish` workflow via Snakemake. Parameters are defined so that they are easy to understand how they can be modified in the config.yml file that executes the pipeline.

generate_jf_count
-----------------

Purpose: Generates a Jellyfish index file that is used for counting k-mers downstream.
Input: Genome reference FASTA file
Output: A genome query hash file containing counts of all k-mers in the form of a .jf file. 

.. code-block:: bash

   jellyfish count -s 3300M -m {params.mer_val} -o {output} {input.fasta_file}


**config.yml parameters**

fasta_file: File path. The genomic reference file used for probe design. Should includes all scaffolds of interest in the provided genome in proper FASTA format. 

mer_val: Integer value. The k-mer size of interest when designing Jellyfish query files, hash tables, and during comparative sequence composition analyses. Default = 18. 


generate_bt2_index
------------------

Purpose: Generates genome Bowtie2 indices which is used for aligning probes to the entire genome of interest.

Input: Genome reference FASTA file

Output: Collection of Bowtie2 indices placed in a Bowtie2 directory of your choosing.

.. code-block:: bash

   bowtie2-build --threads 4 {input} {BOWTIE2_DIR}/{ASSEMBLY}

**config.yml parameters**

ASSEMBLY = String. The name of the genome assembly being used. Example) ASSEMBLY = "CHM13"

BOWTIE2_DIR = File path. The location of where generated Bowtie2 indices will live. Default: 'pipeline_output/01_reference_files/02_bt2_idx/{ASSEMBLY}'

generate_jf_index
-----------------

Purpose: To generate k-mer count index files using the derived jellyfish hash table from the `generate_jf_count` step. Generates independent k-mer count index files for each scaffold. 

Inputs: Genome reference FASTA file (FASTA_FILE). Output Jellyfish hash table generated from the `generate_jf_count` step (JF_INDEXFILE)

Outputs: An output jellyfish k-mer count file containing all k-mers within a selected scaffold and it's corresponding counts (JF_OUT). A file that is used to reference the base position of where the start of each k-mer in the output count file occurs (J_INDEX_OUT). A seperated FASTA file of each selected scaffold (SCAFFOLD_FA_OUT).

.. code-block:: bash

    usage: generate_jf_idx.py [-h] -f FASTA_FILE -j JF_INDEXFILE -c CHR_NAME -f_o
                          SCAFFOLD_FA_OUT -j_o JF_OUT -i J_INDEX_OUT -m
                          MER_VAL


**config.yml parameters**

CHR_NAME = String. Described as sample in config.yml file. Each sample can be one or more scaffolds present in a given genome. Scaffold names should match FASTA file headers. 

Example format in config.yml:

sample:
    - "chr1"

Example format in FASTA file:

>chr1
dnasequence

repeat_ID
---------

Purpose: 
Input:
Output:
 
.. code-block:: bash

    usage: repeat_ID.py [-h] -j JF_COUNT -i INDEX_FILE -chr CHR_NAME -st START
                    [-w WINDOW_LENGTH] [-t THRESHOLD] [-c COMPOSITION_SCORE]
                    -o_b BED_FILE -m MER_LENGTH


