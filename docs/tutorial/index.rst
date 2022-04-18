Tigerfish tutorial
=======================================

In this tutorial you will learn how to install and run Tigerfish with the example files provided. Further documentation will describe in detail the features and parameters used in the workflow.

This tutorial is aimed for all who are interested in learning how to use Tigerfish for their own research to design oligonucleotide FISH probes specific at imaging repetitive DNA in current and emerging genomes. 

The only things you will need to follow are a web browswer, an Internet connection, a GitHub account, and access to remote computing resources.


Getting Started
=======================================

Downloading Tigerfish
=======================================

1. Install conda as needed for your system.


2. Proceed with installing Mamba to assist with snakemake installation, as recommended in the Snakemake installation tutorial.


3. Clone this repo, and create and activate the provide environment for implementing Snakemake for Tigerfish:

Running Tigerfish
=======================================

A full tutorial for running Tigerfish using both repeat discovery and probe design from user provided bed coordinates in included to test the pipeline installation. To run the pipeline using both cases on sample files:

When repeat coordinates are provided:


When repeat discovery is initiated:


When these examples are run, expected outputs may be compared here when coordinates are provided and here for repeat discovery mode.

To survey how postprocessing is done, a seperate analysis pipeline may be invoked using the command below once both items in the main implementation tutorial have been completed.



To deploy Tigerfish on your own data, update the file paths in config.yml with the paths to your genome assembly. For more information on input and output files, please read further documentation here.

This pipeline is implemented using Snakemake, and distributed according to best practices. If you are interested in learning more about Snakemake, please follow their tutorials to learn more about their resources for getting started.

