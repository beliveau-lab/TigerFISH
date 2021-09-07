#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=8.0G
#$ -R y
#$ -l h_rt=200:0:0

source activate tigerfish

qsub -cwd -l mfree=8G -l h_rt=200:0:0 -e ./log -o ./log snake_batch.sge

