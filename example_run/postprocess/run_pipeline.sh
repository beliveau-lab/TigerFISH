# configure file paths
CONFIG_FILE='config.yml'
SNAKE_FILE='../../workflow/postprocess/Snakefile'
CONDA_ENVS='../../../shared_conda_envs'

# activate conda env
source activate snakemake_env

# run pipeline
snakemake --snakefile  $SNAKE_FILE --configfile $CONFIG_FILE \
          --conda-prefix $CONDA_ENVS --use-conda --conda-frontend mamba \
          --jobs $NUM_JOBS --latency-wait $LATENCY_WAIT --restart-times $RESTART_ATTEMPTS 3

echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --dag | dot -Tsvg > pipeline_output/dag.svg
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --dag | dot -Tpdf > pipeline_output/dag.pdf
                
echo -e "Generating pipeline HTML report..."
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --report pipeline_output/report.html

#success
echo -e "\nDONE!\n"

