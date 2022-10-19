
Tigerfish User Tutorials
########################

Overview
--------
These tutorials for **Tigerfish** are intended to provide users with an overview of the tool's features and functionality. Here, users will learn how to run **Tigerfish** using *"Repeat Identification Mode"* and *"Probe Design Mode"* using example cases and on real data in the CHM13 V2.0 genome. Following implementation of the main pipeline, users will also walk through example and real-world cases of the *post-process pipeline*. Following this tutorial, users will have all the resources to run **Tigerfish** for their own research applications. 

This tutorial guide was written in mind for those who are newer to running these scripts in a command line environment. These tutorials will help users ensure that they have installed and are running **Tigerfish** properly on their own command line system. 

So let's get started! |:microscope:| |:tiger:|

Background
==========

The following tutorial is used to cover **Tigerfish** functionality on a toy genome composed of the chr4 D4Z4 and chrX DXZ4 repeats which is derived from the latest version (v.2.0) of the `CHM13 human <https://github.com/marbl/CHM13>`_ genome. 

This reference file can be found within the project's Github repo `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/repeat_discovery_test/data/example.fa>`_. 

Users will download the full CHM13 V2.0 reference file later in this tutorial, to work on real-world data. 

.. image:: imgs/tutorials_summary.png
     :width: 500
     :alt: Tigerfish tutorials overview

For the test genome cases, there are two static directories that exist which will be referenced in this tutorial. First, we will walk through *Repeat Identification Mode* on a Test Genome and this directory is located `_here<https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/repeat_discovery_test>`_. 

Next, we will move through the *Probe Design Mode* directory which can be found `_here<https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/probe_design_test>`_. 

Following execution of both run-modes in the main pipeline, we will proceed with the *post-process pipeline* which exists `_here<https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/postprocess/dxz4_test>`_. 

Lastly, we will move into the real-world examples in the main pipeline with CHM13 V2.0 using *Probe Design Mode* and relevant post-process data on real oligo probes mined from the newly annotated chr9 HSAT repeat. 

It's important to note that these directories are static meaning that the config.yml files have been organized to run Tigerfish for expected behavior in each directory. In other words, the only thing that users will need to execute to validate their results in each directory is the run_pipeline.sh scripts independently. 

At UW Genome Sciences, users have access to the departmental Sun Grid Cluster Engine, which is used to deploy cluster compute jobs. Templates for how to run Tigerfish with such a file system can be found here. Even if you are not using the Sun Grid Engine, Tigerfish may be executed locally and all scripts are currently configured for generalizeable use in a command line environment.

To begin, we will describe some context for how Tigerfish is deployed by covering the config.yml file and appropriate run_pipeline.sh scripts found in each directory.

Config file
===========

**Tigerfish** is a command line workflow that is implemented using Snakemake. This `config.yml <https://github.com/beliveau-lab/TigerFISH/blob/master/example_run/main/main_pipeline/config.yml>`_ file summarizes parameters that users are able to modify for their own research. However, we provide default parameters summarized in our paper for recommended use. These parameters for probe design are also described within this documentation along with definitions on our command line definitions page. 

For the example test files and real-world data, these config.yml files in each of these directories have already been modified for use with example and template files. The only modification required will be to list the path to the CHM13 V2.0 genome in the real-world data example, but we walk through how to do this when we get to this step :). 

Repeat Identification Mode on a Test Genome
-------------------------------------------

*Repeat Identification Mode* is intended to be used when a user provides a given genome FASTA file and is perhaps unsure of where their target repeat regions of interest lie within the genomic sequence. Another valid use case for this option is if a user wants to perform genome-wide probe mining over all regions that Tigerfish deems as repetitive. In this case, the example genome FASTA contains a small subset of the DXZ4 and D4Z4 repeats. 

Here's a walkthrough of all the input files provided to get started with running Tigerfish in this example case:

**Pipeline input**

The test genome file described as **example.fa**: 

.. image:: imgs/repeat_disc_fasta.png
     :width: 500
     :alt: Tigerfish example genome FASTA
     
The test genome chrom.sizes file described as **test_chrom.sizes**:

.. image:: imgs/chrom_sizes_repeat_disc.png
     :width: 500
     :alt: Tigerfish example genome chrom.sizes file
     
The Bowtie2 directories for this test genome reference which are found in the path **data/bt2/** relative to the config.yml file:

.. image:: imgs/bt2_repeat_disc.png
     :width: 500
     :alt: Tigerfish Bowtie2 indices for example genome

**Pipeline output**

All expected output files can be found within `this directory <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/repeat_discovery_test/repeat_ID_output>`_. 

Here, a collection of probes for both repeats found on each scaffold are provided in independent directories.

**Pipeline executables**

The **config.yml** file which has preset parameters that **do not** need to be modified for proper execution:

.. image:: imgs/repeat_discovery_config.png
     :width: 500
     :alt: Tigerfish config.yml file for test genome
     
The **run_pipeline.sh** script is what is used to execute the pipeline:

.. image:: imgs/run_pipeline_repeat_disc.png
     :width: 500
     :alt: Tigerfish run pipeline executable shell script
     
     
To check if the expected output files match to what is generated after you run the pipeline you can use the script **run_check_repeatID.sh**:

.. image:: imgs/check_repeat_disc.png
     :width: 500
     :alt: Tigerfish check if repeat discovery mode outputs are as expected
     
     
**Let's walkthrough**

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

