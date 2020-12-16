#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=60.0G
#$ -R y
#$ -l h_rt=72:0:0

source activate snakemake-tigerfish

python survey_pareto.py
