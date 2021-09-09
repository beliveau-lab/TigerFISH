# Pipeline Steps

Description of workflow for the Tigerfish pipeline.

## Overview

The pipeline consists of six primary stages that are described in further detail below.

1. First, reference files are created including input Jellyfish query files and relevant alignment indices.

2. Next, if repeat identification is invoked, scaffolds are searched for regions of elevated k-mer counts, where repetitive regions are most abundant. This step is optional, and if users provide a .bed file with a particular genomic region of interest, this step will be skipped.

3. When provided with a .bed file, probe design is run on the provided region of interest, making use of reference files.

4. Here, candidate probes are refined by selecting probes that are most likely to bind to their target repeat regions using a k-mer similarity search. As redundant probes with sequence similarity and those with low on-target binding are removed, probes are then rank ordered based on target binding abundance.

5. Probe candidates then undergo alignment based filtering to determine in-silico predicted target binding using the tool NUPACK. If the repeat idenfication step is implemented, all repeat regions surveyed are processed in parallel using snakemake for efficient runtime.

6. Finally, a set of final output files are generated with a summary of probes found and predicted on target binding within each repeat region. Probes of interest may be selected for downstream analyses and ideogram generation. 

<div align="center">
  <a href="#overview_of_pipeline"><img src="./img/tigerfish_steps.png" width="800"></a>
<div>

## 1. Generating reference files

## 2. Repeat identification (optional)

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


