
Tigerfish Command Line Tutorials
################################

insert graphical overview of the two methods

Overview
--------
These tutorials for **Tigerfish** are intended to provide users with an overview of the tool's features and functionality. Here, users will learn how to run **Tigerfish** using repeat identification mode and probe design mode. Additionally, all default parameters, config file usage, and parameters implemented in the workflow will be documented for further exploration. Similarly, the post process pipeline used to provide in silico predictions of probe binding in FISH experiments is documented for user implementation. Following this tutorial, users will have the resources to run **Tigerfish** for their own research applications.

Main Pipeline
=============

The following tutorial is used to cover **Tigerfish** functionality on a toy genome composed of the chr4 D4Z4 and chrX DXZ4 repeats in the latest version (ver. 2) of the `CHM13 human <https://github.com/marbl/CHM13>`_ genome. This reference file can be found within the project's Github repo `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/data/example.fa>`_. 

Insert graphical overview of the DXZ4 and D4Z4 tasks

Config file
===========

**Tigerfish** is a command line workflow that is implemented using Snakemake. This `config.yml <https://github.com/beliveau-lab/TigerFISH/blob/master/example_run/main/main_pipeline/config.yml>`_ file summarizes parameters that users are able to modify for their own research. However, we provide default parameters summarized in our paper for recommended use. These parameters for probe design are also described within this documentation along with definitions on our command line definitions page. 

Parameters for Implementation
-----------------------------

Insert table of parameters used in tutorial

Repeat Identification Mode
--------------------------

Summary/screenshot of how this differs from probe design mode

**Pipeline input**

Quick summary


**Pipeline output**

Quick screenshot and link to the expected output directory

Probe Design Mode
-----------------

Summary/screenshot of how this differs from repeat design mode

**Pipeline input**

Quick summary

**Pipeline output**

Quick screenshot of output probes



Postprocess Pipeline
====================

The Tigerfish postprocess pipeline is intended for analysis of specific oligo probes of interest after Tigerfish has been successfully run. Here, users may take selected probes directly from the final Tigerfish probe output file and generate plots of predicted thermodynamic binding sites for each scaffold. Maps of repeat location on each target scaffold are also generated using `chromoMap <https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html>`_. Output bedgraphs of normalized alignment pileup over 1Mb bins may be useful for other genomic analyses beyond Tigerfish use. Here, collections or individual designed probes are validated to check each probe(s) predicted binding behavior.

Config file
===========

Quick description of how this file differs from that of the previous config file and it's dependency on having main run first to test post process.

Parameters for Implementation
-----------------------------

Insert table of parameters used in the tutorial.

Implementing the workflow
=========================

**Pipeline input**

Summary of what's needed

**Pipeline output**

Basic screenshot of chromomaps, tables, etc for select probes of interest. And link to expected output directory.


