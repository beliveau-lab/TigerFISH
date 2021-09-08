# Pipeline Input

Tigerfish pipeline input specification

### Overview

In the [config.yml](../example_run/config.yml), there are two options for pipeline implementation.

1. Defined coords - This setting may be run when a user is interested in designing oligo probes against a known and given repetitive region. To implement this setting, one must toggle the `defined_coords: "TRUE"` and `repeat_discovery: "FALSE"`.


2. Repeat discovery - This setting may be run when a user is interested in probe design on scaffolds where all identified regions with elevated k-mer counts undergo specific probe design. To implement this setting, one must toggle the `defined_coords: "FALSE"` and `repeat_discovery: "TRUE"`.

### Config file

The paths to the files is specified in the [config.yml](../example_run/config.yml). Sample files to run the Tigerfish using both implementation cases are included in the provided [example directory](../example_run/data/).

For example, to run this pipeline on the chm13 human genome assembly, the [chm13v1.1](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) fasta file was downloaded from [this page](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) and decompressed.

The [config.yml](../example_run/config.yml) file sets all the pipeline parameters for a given run. A new config file should be created to accompany each pipeline run.

**Required parameters for `defined_coords` implementation:**

* `bed_file` is the path to the bed file containing the repeat coordinates of interest for probe design e.g. `dxz4.bed`
* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


**Required parameters for `repeat_discovery` implementation:**

* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


### Default parameters described for `repeat_discovery` and `defined_coords`

The following table summarizes which parameters provided in the [config.yml](../example_run/config.yml) are used in each implementation of the pipeline.

<div align="center">
    <table>
        <thead>
            <tr>
                <th align="center">Parameters</th>
                <th align="center">Used in `repeat_discovery`</th>
                <th align="center">Used in `defined_coords`</th>
                <th align="center">pipeline step</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td align="center">window</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:x:</td>
                <th align="center">repeat_identification</th>
            </tr>
            <tr>
                <td align="center">threshold</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:x:</td>
                <th align="center">repeat_identification</th>
            </tr>
            <tr>
                <td align="center">composition</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:x:</td>
                <th align="center">repeat_identification</th>
            </tr>
            <tr>
                <td align="center">file_start</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:x:</td>
                <th align="center">repeat_identification</th>
            </tr>
            <tr>
                <td align="center">mer_val</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">kmer_filter</th>
            </tr>
            <tr>
                <td align="center">c1_val</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">kmer_filter</th>
            </tr>
            <tr>
                <td align="center">c2_val</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">kmer_filter</th>
            </tr>
            <tr>
                <td align="center">enrich_score</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">probe_mer_filter</th>
            </tr>
            <tr>
                <td align="center">copy_num</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">probe_mer_filter</th>
            </tr>
            <tr>
                <td align="center">mer_cutoff</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">probe_mer_filter</th>
            </tr>
            <tr>
                <td align="center">p_dups</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">probe_count</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">target_sum</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">mer_cutoff</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">bt2_alignments</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">max_off_target</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
            <tr>
                <td align="center">max_pdups_binding</td>
                <td align="center">:heavy_check_mark:</td>
                <td align="center">:heavy_check_mark:</td>
                <th align="center">alignment_filter</th>
            </tr>
        </tbody>
    </table>
</div>

**Default parameter descriptions:**

Here, default parameters for window, threshold, and composition were determined to identify the largest number of repeat regions possible.

**Parameters used in repeat identification**

* `window`:

* `threshold`:

* `composition`:

* `file_start`:

**Parameters used in first probe filter; k-mer filter**

* `mer_val`:

* `c1_val`:

* `c2_val`:

* `enrich_score`:

* `copy_num`:

**Parameters used in second probe binding specificity filter**

* `p_dups`:

* `probe_count`:

* `target_sum`:

* `mer_cutoff`:

* `bt2_alignments`:

* `max_off_target`:

* `max_pdups_binding`:

