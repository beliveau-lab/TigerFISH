# Pipeline Steps

Description of workflow for the Tigerfish pipeline.

## Overview

The pipeline consists of six primary stages that are described in further detail below.

1. First, reference files and indices are created from input files.

2. Next, if repeat identification is invoked, scaffolds are searched for regions of elevated k-mer counts, where repetitive regions are most abundant. This step is optional, and if users provide a .bed file with a particular genomic region of interest, this step will be skipped.

3. When provided with a .bed file, probe design is run on the provided region of interest, making use of reference files.

4. Here, candidate probes are refined by selecting probes that are most likely to bind to their target repeat regions using a k-mer similarity search. As redundant probes with sequence similarity and those with low on-target binding are removed, probes are then rank ordered based on target binding abundance.

5. Probe candidates then undergo alignment based filtering to determine in-silico predicted target binding using the tool NUPACK. If the repeat idenfication step is implemented, all repeat regions surveyed are processed in parallel using snakemake for efficient runtime.

6. Finally, a set of final output files are generated with a summary of probes found and predicted on target binding within each repeat region. Probes of interest may be selected for downstream analyses and ideogram generation. 

<div align="center">
  <a href="#overview_of_pipeline"><img src="./img/tigerfish_steps.png" width="800"></a>
<div>

## 1. Generating reference files

The first step of the workflow is to parse the assembly fasta file and generate `bowtie2` and `jellyfish` indices. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is a NGS aligner with very sensitive parameters and is used downstream to check alignment of candidate probe sequences to the reference genome. To survey k-mer frequency analysis, [jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) acts as a k-mer counter to facilitate with these steps.  User provided names for scaffolds of interest to query are taken into account for proceeding steps and for file name generation.

## 2. Repeat identification (optional)

Tigerfish implements logic that identifies repetitive genomic sequences that are highly enriched for k-mers with high occurrence counts. These counts come from a genome-wide index that is pre-computed by Jellyfish from the provided genome assembly as described in the figure below. Here, Tigerfish identifies regions of elevated k-mer counts based on three user defined parameters. First, the algorithm binarizes the integer counts produced by the Jellyfish query by applying the threshold T (default 10), where counts must be >= T. Second, the algorithm applies a sliding window, W (default = 3000), along the vector to sum the number of bases where count is > T. Third, these summed counts are divided by the window width W and compared to a user defined composition score, C (default = 0.5). When the proportion of counts in W that exceed the threshold T surpasses the defined value of C, the start and end indices where C is true are recorded. These defined coordinates highlight regions of k-mer enriched genomic sequences which are then used for probe design.

<div align="center">
  <a href="#repeat_identification_step"><img src="./img/tigerfish_1.png" width="400"></a>
<div>

## 3. Probe design


<div align="center">
  <a href="#probe_design_step"><img src="./img/tigerfish_2.png" width="400"></a>
<div>

## 4. K-mer based specificity checking

<div align="center">
  <a href="#kmer_specificity_check"><img src="./img/tigerfish_3.png" width="400"></a>
<div>

## 5. Alignment based filtering

<div align="center">
  <a href="#alignment_filtering"><img src="./img/tigerfish_4.png" width="800"></a>
<div>

## 6. Output and postprocessing


