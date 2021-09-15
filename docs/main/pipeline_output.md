# Tigerfish Output

Tigerfish pipeline output specification.

### Overview

| Folder        | Description                                                       |
|---------------|-------------------------------------------------------------------|
| [01_reference_files/](../../example_run/tigerfish_main/expected_pipeline_output/01_reference_files)   | files that may be of use in other pipelines or analyses |
| [02_intermediate_files/](../../example_run/tigerfish_main/expected_pipeline_output/02_intermediate_files)  | large intermediate files, useful when debugging, but otherwise disposable |
| [03_output_files/](../../example_run/tigerfish_main/expected_pipeline_output/03_output_files) | all DNA/RNA probe sets as .tsv files |

All unwanted files may be safely deleted once pipeline is run and completed.

**NOTE:** To minimize disk usage, it may be recommended to keep only files found after probes are designed under the intermediate directories. Some of these Jellyfish count files can be quite large if running at the scale of hg38 and chm13 human genomes.

### Probes

#### File locations

| Item        | Location                                                       |
|---------------|-------------------------------------------------------------------|
| DNA FISH probes in .tsv format | [03_output_files/01_dna_probes/](../../example_run/expected_pipeline_output/03_output_files/01_dna_probes) | 

#### Repetitive DNA probe sets

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

#### Summary of probes within each target repeat identified

The columns in this file are:

| # | Column | Description |
|---|--------|-------------|
| 0 | repeat coords | the repeat coordinates from the DNA probe file |
| 1 | probe count | the total number of probes found in each repeat region |
| 2 | on target sum | the total sum of all values for a repeat from col 9 |
| 3 | off target sum | the total sum of all values for a repeat from col 10 |


### Reporting

An HTML report with diagnostics and detailed pipeline information can by generated with the following command:

```
$ snakemake --snakefile path/to/Snakefile --configfile path/to/config.yml --report pipeline_output/report.html
```

An example report is available in the example output tutorial. For a visualization of the pipeline DAG structure, see: [pipeline.pdf](../../example_run/expected_pipeline_output/pipeline.pdf) or [pipeline.svg](../../example_run/expected_pipeline_output/pipeline.svg)
