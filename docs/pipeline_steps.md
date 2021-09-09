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
</div>

## 1. Generating reference files

The first step of the workflow is to parse the assembly fasta file and generate `bowtie2` and `jellyfish` indices. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is a NGS aligner with very sensitive parameters and is used downstream to check alignment of candidate probe sequences to the reference genome. To survey k-mer frequency analysis, [jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) acts as a k-mer counter to facilitate with these steps.  User provided names for scaffolds of interest to query are taken into account for proceeding steps and for file name generation.

## 2. Repeat identification (optional)

Tigerfish implements logic that identifies repetitive genomic sequences that are highly enriched for k-mers with high occurrence counts. These counts come from a genome-wide index that is pre-computed by Jellyfish from the provided genome assembly as described in the figure below. Here, Tigerfish identifies regions of elevated k-mer counts based on three user defined parameters. First, the algorithm binarizes the integer counts produced by the Jellyfish query by applying the threshold T (default 10), where counts must be >= T. Second, the algorithm applies a sliding window, W (default = 3000), along the vector to sum the number of bases where count is > T. Third, these summed counts are divided by the window width W and compared to a user defined composition score, C (default = 0.5). When the proportion of counts in W that exceed the threshold T surpasses the defined value of C, the start and end indices where C is true are recorded. These defined coordinates highlight regions of k-mer enriched genomic sequences which are then used for probe design.

<div align="center">
  <a href="#repeat_identification_step"><img src="./img/tigerfish_1.png" width="400"></a>
</div>

## 3. Probe design

Candidate probe sequences are mined from the genome using [OligoMiner](https://github.com/beliveau-lab/OligoMiner). Here, a multi-entry fasta is produced to then identify specificity of each candidate probe to the region it was designed against or derived from.

<div align="center">
  <a href="#probe_design_step"><img src="./img/tigerfish_2.png" width="400"></a>
</div>

## 4. K-mer based specificity checking

After probes are mined, each candidate then undergoes specificity filtering by checking k-mer similarity of each probe to that of the target sequence. These steps are briefly described in the image below, but are also discussed to greater depth in the Tigerfish preprint.

Briefly, the occurence of all k-mers within a given candidate probe are counted within a given target sequence over that of the whole queried genome. Probes are then rank ordered in descending order using a normalization value using two constants c1 and c2 which were applied to the sum of on target 18-mers within each repeat region and the proportion of k-mer binding for each probe, respectively. By rank ordering probes by these two values, this prioritizes probes that are likely to demonstrate high on target repeat binding as shown when we implement Bowtie2 to assess off target potential.

<div align="center">
  <a href="#kmer_specificity_check"><img src="./img/tigerfish_3.png" width="400"></a>
</div>

## 5. Alignment based filtering

Remaining probes then proceed with a filtering approach that independently implements Bowtie2 on each probe against sequences in the human genome. Here, alignments are returned as a BAM file (citation) for each probe which is then processed from the resulting SAM file using [SAMtools](https://samtools.github.io). Using this SAM file, [sam2pairwise](https://github.com/mlafave/sam2pairwise) is called to return derived alignment sequences that are returned for each oligo probe. With these provided pairs of probe sequence and derived alignment sequence, we implement [NUPACK 4.0](http://www.nupack.org) to compute the predicted thermodynamic likelihood that each alignment pair will form duplexes under FISH conditions. Here, we compute a predicted on target alignment score which is determined by taking the sum of all predicted duplexing scores for derived alignments that are found within the target repeat region. Off target alignment scores are computed by taking the aggregate sum of all predicted duplexing scores from derived alignments that are found outside the target repeat. The predicted in silico on target binding for each oligo is then computed as a proportion of on-target binding over all described predicted duplexing scores over all alignments. After the first probe is identified as a valid repetitive specific probe that satisfies a value below the maximum predicted off target binding sum, all probes that are sequentially surveyed for valid on target binding, are compared to one another using NUPACK to survey the likelihood that probes in the final set will result in secondary structure under FISH conditions.


<div align="center">
  <a href="#alignment_filtering"><img src="./img/tigerfish_4.png" width="800"></a>
</div>

## 6. Output and postprocessing

Lastly, Tigerfish provides a file of probes found within each repeat region surveyed in a .tsv format, along with aggregate on target and off target binding sums based on thermodynamic prediction for all probes found within each repeat region. If repeat identification is implemented over a scaffold, a final output file contains all candidate probes that have satisfied specificity filtering. Candidates from this file may then be taken to downstream postprocessing.

In the postprocessing pipeline, users may select particular probes of interest and generate further information regarding in-silico binding predictions as well as chromomaps where anticipated target probes are likely to be imaged using metaphase spreads when conducting FISH. Further documentation about the postprocessing pipeline may be found here.

