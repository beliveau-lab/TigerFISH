# Pipeline Input

Tigerfish pipeline input specification

### Overview

In the [config.yml](../example_run/config.yml), there are two options for pipeline implementation.

1. Defined coords - This setting may be run when a user is interested in designing oligo probes against a known and given repetitive region. To implement this setting, one must toggle the `defined_coords: "TRUE"` and `repeat_discovery: "FALSE"`.


2. Repeat discovery - This setting may be run when a user is interested in probe design on scaffolds where all identified regions with elevated k-mer counts undergo specific probe design. To implement this setting, one must toggle the `defined_coords: "FALSE"` and `repeat_discovery: "TRUE"`.

### Config file

The paths to the files is specified in the [config.yml](../example_run/config.yml). Sample files to run the Tigerfish using both implementation cases are included in the provided [example genome assembly](../example_run/data/).

For example, to run this pipeline on the chm13 human genome assembly, the [chm13v1.1](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) fasta file was downloaded from [this page](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) and decompressed.

The [config.yml](../example_run/config.yml) file sets all the pipeline parameters for a given run. A new config file should be created to accompany each pipeline run.

**Required parameters for `defined_coords` implementation:**

* `bed_file` is the path to the bed file containing the repeat coordinates of interest for probe design e.g. `dxz4.bed`
* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


**Required parameters for `repeat_discovery` implementation:**

* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


**Default parameters described for `repeat_discovery` and `defined_coords` implementations:**

