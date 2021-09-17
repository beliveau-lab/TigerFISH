# Pipeline Input

Tigerfish pipeline input specification

### Overview

In the config.yml, there are two options for pipeline implementation that are summarized with links to appropriate config file paths.

1. [Defined coords](../../example_run/main/defined_coords/config.yml) - This setting may be run when a user is interested in designing oligo probes against a known and given repetitive region. To implement this setting, one must toggle the `defined_coords: "TRUE"` and `repeat_discovery: "FALSE"`.

2. [Repeat discovery](../../example_run/main/repeat_ID/config.yml) - This setting may be run when a user is interested in probe design on scaffolds where all identified regions with elevated k-mer counts undergo specific probe design. To implement this setting, one must toggle the `defined_coords: "FALSE"` and `repeat_discovery: "TRUE"`.

### Config file

The paths to the files is specified in the config.yml for each run mode. Sample files to run the Tigerfish using both implementation cases are included in the provided in each [example directory](../../example_run/main/defined_coords/data/).

For example, to run this pipeline on the chm13 human genome assembly, the [chm13v1.1](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) fasta file was downloaded from [this page](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz) and decompressed.

The config.yml file in both forms of pipeline implementation sets all the pipeline parameters for a given run. A new config file should be created to accompany each pipeline run.

**Required parameters for `defined_coords` implementation:**

* `bed_file` is the path to the bed file containing the repeat coordinates of interest for probe design e.g. `dxz4.bed`
* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


**Required parameters for `repeat_discovery` implementation:**

* `samples` is the name of the scaffold that the bed coordinates corresponds to and is provided as a string e.g. `chrX`. This can also include a list of samples for multiple genome scaffolds.
* `fasta_file` is the path to the genome sequence file e.g. `chm13v_1_1.fa`


### Default parameters described for `repeat_discovery` and `defined_coords`

The following table summarizes which parameters provided in each config file are used in each implementation of the pipeline.

<div align="center">
    <table>
        <thead>
            <tr>
                <th align="center">Parameters</th>
                <th align="center">Used in repeat_discovery</th>
                <th align="center">Used in defined_coords</th>
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

**Parameters used in repeat identification**

Here, default parameters for window, threshold, and composition were determined to identify the largest number of repeat regions possible. For more information on how each pipeline step works, please read the [pipeline steps documentation](pipeline_steps.md).

* `window`: The size of this window is used to evaluate regions of elevated k-mer counts. The sum of all binarized threshold values within this window are taken in this window (default = 3000).

* `threshold`: Tigerfish binarizes the integer counts produced by the [Jellyfish](https://github.com/gmarcais/Jellyfish) query by applying this threshold (default = 10), where counts must be >= this value.

* `composition`: The summed counts within the window are divided by the length of the window variable. This value is then compared to the composision score (default = 0.50). When this proportion exceeds the defined composition value, start and end indices for regions of k-mer enriched genomic sequences are stored.

* `file_start`: This value is used to begin searching genomic regions at the start of scaffolds (default = 0).

**Parameters used in first probe filter; kmer_filter**

* `mer_val`: This value is used to generated Jellyfish query files, hash indices, and is used to split identified repeat regions and designed probe sequences into appropriate matching k-mer sizes to identify probes with high on-target sequence specificity (default = 18).

* `c1_val`: A constant used to derive a normalized specificity value which takes into account the sum of all k-mer matches for a given probe within an intended target repeat region (default = 5). This normalized value is used along with `c2_val` to rank probes that are more likely to be high, on-target, probe candidates.

* `c2_val`: A constant used to derive a normalized specificity value which takes into account the sum of all k-mer matches for a given probe within an intended target repeat region over all k-mer matches within the queried genome, this proportion is also known as the `enrich_score` (default = 1).

**Parameters used in probe_mer_filter**

* `enrich_score`: Here, this value is defined as the sum of all k-mer matches for a given probe within an intended target repeat region over all k-mer matches within the queried genome (default = 0.5). In order for probes to proceed with downstream filtering, they must exceed or be equivalent to this defined value.

* `copy_num`: This is the sum of all k-mer matches for a given probe within an intended target repeat region (default = 10). In order for probes to proceed with downstream filtering, they must exceed or be equivalent to this defined value.

* `mer_cutoff`: If a probe that shares at least this proportion of k-mers (default = 0.95) after rank ordering a probe set, lower ranking probes are culled.

**Parameters used in alignment_filter**

* `p_dups`: A determined proportion based on taking each candidate probe and identifying the predicted likelihood of binding using [NUPACK](http://www.nupack.org). This value is determined by taking the sum of all predicted duplexing (pdups) values for all on-target probe alignments within a candidate repeat region over all predicted duplexing values (default = 0.95). Here, values must be equal to or exceed default value for probes to be further evaluated in specificity.

* `probe_count`: To run the `alignment_filter` step, there are two main run settings which is based on total probe count per repeat region or by aggregate target binding sum for each repeat. If one is interested in returning up to X best probes for a given repeat region, you may adjust this value accordingly (default = 0).

* `target_sum`: Using this run setting, the aggregate on target binding sum based on alignments as determined by NUPACK is used as a runtime cutoff. When a user defined `target_sum` is set (default = 5000), the `alignment_filter` script will continue to search and store candidate probes until this threshold is exceeded.

* `bt2_alignments`: When generating alignments for candidate probes to the entire queried genome, users can specify the maximum number of alignments that bowtie2 can return (default=300000).

* `max_off_target`: This value is the total sum of all predicted duplexing values where a given probe has alignments to regions in the queried genome outside of the target repeat. Users can set a max off-target value when storing probe candidates (default = 100).

* `max_pdups_binding`: As probes undergo alignment, for probes that are stored, NUPACK is used to determined the likelihood that candidate probes will bind to other potential candidates. If a probe has a predicted binding score that is less than the `max_pdups_binding` value (default = 0.60), that probe is also stored in the final set.

