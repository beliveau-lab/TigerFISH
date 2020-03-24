#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=80.0G
#$ -R y
#$ -l h_rt=90:0:0

if [[ -e /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs
        module load python/3-anaconda
        module load bedtools/latest
fi

for file in *_regions.fa ; do
        chr=${file##*/}
        chr=${chr%_*}
	python design_probes.py -f $file -chr $chr
done

#python design_probes.py -f chr19_regions.fa 
