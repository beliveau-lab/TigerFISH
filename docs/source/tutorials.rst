
Tigerfish User Tutorials
################################

Overview
--------
These tutorials for **Tigerfish** are intended to provide users with an overview of the tool's features and functionality. Here, users will learn how to run **Tigerfish** using repeat identification mode and probe design mode. Additionally, all default parameters, config file usage, and parameters implemented in the workflow will be documented for further exploration. 

Similarly, the post process pipeline used to provide in silico predictions of probe binding in FISH experiments is documented for user implementation. Following this tutorial, users will have the resources to run **Tigerfish** for their own research applications.

These tutorials were also written in mind for those who are newer to running these scripts in a command line environment! Here, you will follow tutorials on how to make sure that you have installed tool correctly and this will walk you through test cases and how to get Tigerfish started on generating real data in the CHM13 v2.0 genome. 

So let's get started! |:microscope:| |:tiger:|

Main Pipeline
=============

The following tutorial is used to cover **Tigerfish** functionality on a toy genome composed of the chr4 D4Z4 and chrX DXZ4 repeats in the latest version (ver. 2) of the `CHM13 human <https://github.com/marbl/CHM13>`_ genome. This reference file can be found within the project's Github repo `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/data/example.fa>`_. 

.. image:: imgs/tutorials_summary.png
     :width: 500
     :alt: Tigerfish tutorials overview

Config file
===========

**Tigerfish** is a command line workflow that is implemented using Snakemake. This `config.yml <https://github.com/beliveau-lab/TigerFISH/blob/master/example_run/main/main_pipeline/config.yml>`_ file summarizes parameters that users are able to modify for their own research. However, we provide default parameters summarized in our paper for recommended use. These parameters for probe design are also described within this documentation along with definitions on our command line definitions page. 

Repeat Identification Mode
--------------------------

Repeat Identification mode can simply be activated when `repeat_ID` = TRUE and `probe_design` = FALSE in the config.yml file provided in the working directory to run Tigerfish. An example of what this looks like in the header of the config.yml file is shown below:

.. code-block:: bash

    #path to genome fasta
    fasta_file: "data/example.fa"

    #path to file containing primary chromosome sizes
    chrom_sizes_file: "data/chm13.chrom.sizes"

    #if coordinates are provided for probe design, file goes here
    #bed_file: "data/dxz4_synthetic.bed"

    #option for probe design that directs pipeline implementation
    defined_coords: "FALSE"
    repeat_discovery: "TRUE"

Here, the goal of this pipeline is to identify valid repeat regions within the provided toy genomic FASTA to design valid imaging probes against them. 

**Pipeline input**

Here, users provide a toy genomic FASTA file containing hte DXZ4 and D4Z4 repeats along with a chrom.sizes file summarizing the lengths of these test scaffolds.

**Pipeline output**

All expected output files can be found within `this directory <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/expected_pipeline_output/repeat_ID_output>`_. Here, a collection of probes for both repeats found on each scaffold are provided in independent directories.

Probe Design Mode
-----------------

Here, a test file approximating the location of the DXZ4 within this test FASTA file is provided so that the `repeat_discovery` mode is skipped so probes are directly mined and mapped within the defined coordinates provided. The specifications for the config.yml file should be modified as such:


.. code-block:: bash

    #path to genome fasta
    fasta_file: "data/example.fa"

    #path to file containing primary chromosome sizes
    chrom_sizes_file: "data/chm13.chrom.sizes"

    #if coordinates are provided for probe design, file goes here
    bed_file: "data/dxz4_synthetic.bed"

    #option for probe design that directs pipeline implementation
    defined_coords: "TRUE"
    repeat_discovery: "FALSE"


**Pipeline input**

In addition to the genomic FASTA and chrom.sizes file, users also specify that a BED file containing the coordinates of the repeat(s) of interest are provided.

**Pipeline output**

Similar to that of `repeat_discovery` mode, an independent directory contains the probes of interest that map to the repeat region provided in the input BED file. This output directory can be found `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/expected_pipeline_output/repeat_ID_output>`_. 



Postprocess Pipeline
====================

The Tigerfish postprocess pipeline is intended for analysis of specific oligo probes of interest after Tigerfish has been successfully run. Here, users may take selected probes directly from the final Tigerfish probe output file and generate plots of predicted thermodynamic binding sites for each scaffold. Maps of repeat location on each target scaffold are also generated using `chromoMap <https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html>`_. Output bedgraphs of normalized alignment pileup over 1Mb bins may be useful for other genomic analyses beyond Tigerfish use. Here, collections or individual designed probes are validated to check each probe(s) predicted binding behavior.

Config file
===========

This `config.yml <https://github.com/beliveau-lab/TigerFISH/blob/master/example_run/postprocess/config.yml>`_ file summarizes parameters that users are able to modify for their own research. This workflow is also implemented in Snakemake and provides example outputs that users may compare.

Implementing the workflow
=========================

**Pipeline input**

To implement the post process workflow, users must provide a probe file that was derived as the output from the main workflow. An example probe that is used for testing is one generated for DXZ4. Here, users may provide collections of probes that map to the same repeat, or those that map to different repeats on different scaffolds of interest. 

Users must also provide the directory for where Bowtie2 indices were generated from the main pipeline and a chrom.sizes file. These test files may be viewed within the provided paths shown within the config.yml provided.

**Pipeline output**

Here, users will receive a directory containing genome wide binding maps of aggregate binding for each chromosome repeat target, a summary of which genome bins map to binding signal reported by thermodynamic data, as well as a chromoMap to demonstrate where binding is anticipated to occur during a FISH experiment. These expected outputs can be found `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/postprocess/expected_pipeline_output>`_, for user comparison.

