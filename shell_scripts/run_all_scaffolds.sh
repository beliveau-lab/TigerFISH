#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=90.0G
#$ -R y
#$ -l h_rt=90:0:0

if [[ -e /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs
        module load python/3-anaconda
	module load bedtools/latest
fi 

for file in ../../../beliveau/vol1/reference/Assemblies/hg38_bp_noheader/*.fa ; do
	chr=${file##*/}
	chr=${chr%.*}	
	python repeat_ID.py -j hg38_build_18mers.jf -f $file -chr $chr -s 3000 -t 10 -c 0.50 -schr $file -st 0
done

echo "all done"
