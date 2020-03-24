#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=8.0G
#$ -R y
#$ -l h_rt=72:0:0

if [[ -e /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs
        module load bedtools/latest
fi 

qsub -cwd -l mfree=8G -l h_rt=72:0:0 -e ./log -o ./log snake_batch.sge
