# Pipeline Input

Tigerfish postprocess pipeline input specification

### Overview

The Tigerfish postprocess pipeline is intended for analysis of specific oligo probes of interest after Tigerfish has been successfully run. Here, users may take selected probes directly from the final Tigerfish probe output file and generate plots of predicted thermodynamic binding sites for each scaffold. Maps of repeat location on each target scaffold are also generated using [chromoMap](https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html). Output bedgraphs of normalized alignment pileup over 1Mb bins may be useful for other genomic analyses beyond Tigerfish use. 

### Config file

The paths for reference files is specified in the [config.yml](../../example_run/config.yml). Samples to run the Tigerfish postprocessing pipeline are included in the provided [example directory](../../example_run/data/).

**Required parameters for implementation**

* `genome_windows` is a file containing the lengths of each scaffold provided in a given genome assembly that is referred to as a `chrom.sizes` file. Here, for this example, the `chrom.sizes` file for the chm13 genome can be found here and is provided in this example exercise here.

* `probe_file` is a file containing a subset of the selected probes of interest after the Tigerfish main pipeline has been successfully run and completed. The exact file format and columns of the .tsv contain the following columns. For this exercise, subset probes from the probe file generated after the Tigerfish main run example is used. For reference to the full original file, that may be found here. For further descriptions about the output of Tigerfish and run modes, please refer to the guides provided in the Tigerfish main output and Tigerfish steps markdowns.

The columns in DNA probe files are: 

| # | Column | Description |
|---|--------|-------------|
| 0 | probe coords | the chrom:start-stop of each probe |
| 1 | repeat coords | the chrom:start-stop of each repeat |
| 2 | probe | DNA sequence of the probe |
| 3 | Tm | Melting temperature of the probe DNA sequence |
| 4 | k-mer count in repeat target | sum of all probe k-mer counts within target repeat |
| 5 | k-mer count in whole genome | sum of all probe k-mer counts within queried genome |
| 6 | on-targe k-mer binding proportion | the values of column 4 over column 5 |
| 7 | normalized rank ordering value | used to rank order probes within a target repeat, considers cols 4 and 6 |
| 8 | total on-target predicted duplexing | aggregate sum of predicted duplexing values within a target repeat |
| 9 | total off-target predicted duplexing | aggregate sum of predicted duplexing values outside a target repeat |
| 10 | proportion of on-target predicted duplexing | the values of column 8 over column 9 |


* `bt2_alignments`, when generating alignments for candidate probes to the entire queried genome, users can specify the maximum number of alignments that bowtie2 can return (default=300000).

###


