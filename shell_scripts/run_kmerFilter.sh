#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=10G

for file in /net/gs/vol1/home/eaguil/beliveau/2019_4_29_TRF_output/candidate_probes/*.bed ;do 
	basename=`basename $file`
	echo $basename
	python /net/beliveau/vol1/project/OligoMiner/kmerFilter.py -f $file -j /net/beliveau/vol1/reference/jf_files/hg38_18.jf -m 18 -k 490
done

