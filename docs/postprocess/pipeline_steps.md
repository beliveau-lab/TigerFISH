# Pipeline Steps

Description of workflow for the Tigerfish postprocess pipeline

## Overview

The pipeline consists of five primary stages that are described in further detail below. 

1. First, reference files of the binned queries genome are created from input files.

2. Next, bed files are generated from probes of interest what were candidates from the output of the Tigerfish main pipeline. Each probe is processed independently and in parallel.

3. Each probe undergoes alignment with the reference genome and thermodynamic probabilities of probe to target alignment are computed.

4. Bedtools intersect is implemented to identify where predicted pileup of predicted binding events occur in the binned reference genome.

5. Predicted binding events are aggregated, normalized, and plotted for easy visual reference. Repetitive DNA targets are also visualized on scaffolds using the chromoMap R package.


## 1. Binning the queried genome

The first step in the pipeline is to parse a provided genome assembly fasta file to bin the genome into 1 Mb bins using [bedtools](https://bedtools.readthedocs.io/en/latest/). This binned genome file is used downstream for creating visualization of aggregate pileups of predicted thermodynamic binding sites in the genome.

## 2. Generate approriate probe bed files

For selected probes of interest from the output `.tsv` file of the main Tigerfish pipeline, these probe coordinates are organized in a `.bed` format, where each probe is split into a seperate file. Each probe is split seperately to generate alignments and visualizations independently.

## 3. Align probes to the reference genome individually

Here, bowtie2 is implemented to generate alignments of each target probe to a provided query genome. Here, users have the option to specify the maximum number of returned alignments for each probe target. For each of these alignments, we implement NUPACK to predict the thermodynamic likelihood that the probe of interest will likely bind to the derived alignment sequence. This computed value is stored for all alignments corresponding to a given probe target.

## 4. Intersect alignments with genome bins

Then, all returned alignment targets in the genome are processed into a `.bed` format, where bedtools is implemented to intersect alignment targets with the genomic bins generated in part 1 of the workflow.

## 5. Generating in-silico predicted target maps and chromoMaps of repeat targets

By mapping to where each of these alignments correspond to genomic bins, this allows for aggregation of predicted thermodynamic binding sites to be done for each 1 Mb genomic bin. Using this framework, these aggregated sums are normalized using [sklearn.preprocessing](https://scikit-learn.org/stable/modules/preprocessing.html) to standardize values across genomic bins before being plotted using matplotlib to generate barplots representing binding enrichment over each scaffold. For each target repeat, a chromoMap is generated to demonstrate visually where the probe target should be found within a given scaffold. Bedgraphs generated from these analyses are also provided and may be used for further analyses on a genome browser beyond Tigerfish use.

