#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=25.0G
#$ -pe serial 3
#$ -R y
#$ -l h_rt=120:0:0

if [[ -e /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs
        module load bedtools/latest
fi

python kmer_decomp_filter.py -f hg38_probes_filtered_2019_11_20.bed -k 10 -o hg38_probes_filter_20_20_02_k10
