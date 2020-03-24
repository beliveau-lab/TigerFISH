#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=40.0G
#$ -R y
#$ -l h_rt=90:0:0

if [[ -e /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs
        module load bedtools/latest
fi

python repeat_ID.py -j chrX_jf_temp.txt -f chrX.fa -chr "chrX" -s 3000 -t 10 -c 0.5 -st 0 -schr ../../../../reference/Assemblies/hg38_bp_noheader/chrX.fa

#for file in chr19_subset.fa; do
#        python repeat_ID.py -j i/net/beliveau/vol1/home/eaguil/beliveau/2019_10_28_hg38_run/hg38_build_18mers.jf -f $file -chr chr19 -s 3000 -t 10 -c 0.50 -schr $file -st 0
#done
