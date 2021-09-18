
# Tigerfish Output

Tigerfish postprocessing pipeline output specification.

### Overview

| Folder        | Description                                                       |
|---------------|-------------------------------------------------------------------|
| [01_reference_files/](../../example_run/postprocess/expected_pipeline_output/01_reference_files/)   | files that may be of use in other pipelines or analyses |
| [02_intermediate_files/](../../example_run/postprocess/expected_pipeline_output/02_intermediate_files)  | large intermediate files, useful when debugging, but otherwise disposable |
| [03_output_files/](../../example_run/postprocess/expected_pipeline_output/03_output_files) | all DNA probe alignment and chromomaps |

All unwanted files may be safely deleted once pipeline is run and completed.

### Output Analyses

#### File locations

#### In-silico probe binding prediction maps

| Item        | Location                                                       |
|---------------|-------------------------------------------------------------------|
| DNA FISH probe predicted target binding | [03_output_files/01_plot_alignment](../../example_run/postprocess/expected_pipeline_output/03_output_files/01_plot_alignment/) | 


#### Scaffold ideograms using chromoMap

| Item        | Location                                                       |
|---------------|-------------------------------------------------------------------|
| Target repeat chromoMap plots | [03_output_files/02_chromomap](../../example_run/postprocess/expected_pipeline_output/03_output_files/02_chromomap/) |

### Reporting

An HTML report with diagnostics and detailed pipeline information can by generated with the following command:

```
$ snakemake --snakefile path/to/Snakefile --configfile path/to/config.yml --report pipeline_output/report.html
```

An example report is available in the example output tutorial. For a visualization of the pipeline DAG structure, see: [pipeline.pdf](../../example_run/postprocess/expected_pipeline_output/pipeline.pdf) or [pipeline.svg](../../example_run/postprocess/expected_pipeline_output/pipeline.svg)

