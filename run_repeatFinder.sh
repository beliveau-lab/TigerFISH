#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=30.0G

for file in trf_alt_output/trf_output/*.dat ;do
	basename=`basename $file`
        echo $basename
        python repeatFinder.py $file
done
