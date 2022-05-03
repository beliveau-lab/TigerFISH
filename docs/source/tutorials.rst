
Tigerfish Command Line Tutorials
################################

Overview
--------
These tutorials for **Tigerfish** are intended to provide users with an overview of the tool's features and functionality. Here, users will learn how to run **Tigerfish** using repeat identification mode and probe design mode. Additionally, all default parameters, config file usage, and parameters implemented in the workflow will be documented for further exploration. Similarly, the post process pipeline used to provide in silico predictions of probe binding in FISH experiments is documented for user implementation. Following this tutorial, users will have the resources to run **Tigerfish** for their own research applications.

Main Pipeline
#############

The following tutorial is used to cover **Tigerfish** functionality on a toy genome composed of the chr4 D4Z4 and chrX DXZ4 repeats in the latest version (ver. 2) of the `CHM13 human <https://github.com/marbl/CHM13>` genome. This reference file can be found within the project's Github repo `here <https://github.com/beliveau-lab/TigerFISH/tree/master/example_run/main/main_pipeline/data/example.fa>`. 

Config file
-----------

**Tigerfish** is a command line workflow that is implemented using Snakemake. This `config.yml <https://github.com/beliveau-lab/TigerFISH/blob/master/example_run/main/main_pipeline/config.yml>` file summarizes parameters that users are able to modify for their own research. However, we provide default parameters summarized in our paper for recommended use. These parameters for probe design are also described within this documentation here. Definitions of all parameters are shown below. 

**Parameters for Implementation**

`fasta_file`:

`chrom_sizes_file`:

`defined_coords`:

`repeat_discovery`:

`assembly`:

`samples`:

`window`:

`threshold`:

`composition`: 

`file_start`:

`min_length`:

`max_length`:

`min_temp`:

`max_temp`:

`mer_val`:

`c1_val`:

`c2_val`:

`enrich_score`:

`copy_num`:

`genome_windows`:

`target_sum`:

`off_bin_thresh`:

`mer_cutoff`:

`bt2_alignments`:

`max_pdups_binding`:

`seed_length`:

`model_temp`:

`min_on_target`:

`max_probe_return`:


Repeat Identification Mode
--------------------------

**Pipeline input**


**Pipeline output**


Probe Design Mode
-----------------


**Pipeline input**


**Pipeline output**


Postprocess Pipeline
####################

The Tigerfish postprocess pipeline is intended for analysis of specific oligo probes of interest after Tigerfish has been successfully run. Here, users may take selected probes directly from the final Tigerfish probe output file and generate plots of predicted thermodynamic binding sites for each scaffold. Maps of repeat location on each target scaffold are also generated using `chromoMap <https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html`>_. Output bedgraphs of normalized alignment pileup over 1Mb bins may be useful for other genomic analyses beyond Tigerfish use. Here, collections or individual designed probes are validated to check each probe(s) predicted binding behavior.

Config file
-----------

`genome_windows`:

`bowtie2_dir`:

`probe_file`:

`assembly`:

`samples`:

`bt2_alignments`:

`window_size`:

`model_temp`:

`seed_length`:

Implementing the workflow
-------------------------


**Pipeline input**


**Pipeline output**


