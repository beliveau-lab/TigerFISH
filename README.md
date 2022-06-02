
<p align="center">
  <img width="300" height="300" src="https://user-images.githubusercontent.com/46250421/132409143-8072bc08-c21f-407d-bc58-9019701a0c7f.png">
</p>

# TigerFISH Pipeline

## Overview

Tigerfish is a technology and software tool that enables users to design oligonucleotide FISH experiments specific for satellite DNA and repetitive DNA families at the scale of genomes. Tigerfish is comprised of the following resources for repetitive probe design:

1. A pipeline that allows for repetitive probe discovery and design. Users may provide a particular region of interest for oligo probe design, or provide a genomics scaffold for repeat discovery and probe design.

2. A post processing pipeline, where users may perform analysis on probes of interest to determine on and off target probe specificity in silico and generate karyoplots using ChromoMaps.

3. An interactive web application for probe design that will include diverse model organism genomes containing information about relevant repetitive DNA probes for FISH experiments (FISHtank, forthcoming Q1 2022).  

Our documents are publicly found and hosted on Read the Docs.

## Installation

1. Install [conda](https://docs.conda.io/en/latest/miniconda.html) as needed for your system.

2. Proceed with installing Mamba to assist with snakemake installation, as recommended in the [Snakemake installation tutorial](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

```
$ conda install -n base -c conda-forge mamba
```

3. Clone this repo, and create and activate the provide environment for implementing Snakemake for Tigerfish:

```
$ git clone https://github.com/beliveau-lab/TigerFISH.git \
    && cd TigerFISH/ \
    && mamba env create -f snakemake_env.yml \
    && conda activate snakemake_env
```

## Running the pipeline

A full tutorial for running Tigerfish using both repeat discovery and probe design from user provided bed coordinates in included to test the pipeline installation. To run the pipeline using both cases on sample files:

```
$ cd example_run/main/main_pipeline/ && . run_pipeline.sh
```

Repeat_ID mode is set at the default in the [config.yml](example_run/main/main_pipeline/config.yml) file, but the defined_coords mode can also be set to True to enable this run mode. Please see our config.yml file and Read the Docs documentation to learn more about parameters and Tigerfish run modes.


When these examples are run, expected outputs may be compared [here](example_run/main/main_pipeline/expected_pipeline_output/defined_coords_output/) when coordinates are provided and [here](example_run/main/main_pipeline/expected_pipeline_output/repeat_ID_output/) for repeat discovery mode.

To survey how postprocessing is done, a seperate analysis pipeline may be invoked using the command below once both items in the main implementation tutorial have been completed. 

```
$ cd example_run/postprocess/ && . run_pipeline.sh
```
To deploy Tigerfish on your own data, update the file paths in config.yml with the paths to your genome assembly. For more information on input and output files, please read further documentation here.

This pipeline is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), and distributed according to [best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). If you are interested in learning more about Snakemake, please follow their tutorials to learn more about their resources for getting started. 

## Documentation

Please read our comprehensive [Read the Docs page and tutorials to learn more about Tigerfish](https://beliveau-lab-tigerfish.readthedocs-hosted.com/en/latest/).

## Questions

If you have questions or issues, please open an issue on GitHub

## Citation

For usage of the pipeline, please cite according to citation.bib

## License

We provide this open source software without any warranty under the MIT license.

