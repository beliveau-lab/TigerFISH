# Pipeline Steps

Description of workflow for the Tigerfish postprocess pipeline

## Overview

The pipeline consists of five primary stages that are described in further detail below. 

1. First, reference files of the binned queries genome are created from input files.

2. Next, bed files are generated from probes of interest what were candidates from the output of the Tigerfish main pipeline. Each probe is processed independently and in parallel.

3. Each probe undergoes alignment with the reference genome and thermodynamic probabilities of probe to target alignment are computed.

4. Bedtools intersect is implemented to identify where predicted pileup of predicted binding events occur in the binned reference genome.

5. Predicted binding events are aggregated, normalized, and plotted for easy visual reference.


## 1. Binning the queried genome

## 2. Generate approriate probe bed files

## 3. Align probes to the reference genome individually

## 4. Intersect alignments with genome bins

## 5. Generating in-silico predicted target maps and chromoMaps of repeat targets

