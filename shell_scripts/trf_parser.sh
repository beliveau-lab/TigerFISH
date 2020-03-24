#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=32.0G

###4/22/2019
###Robin Aguilar
###Beliveau Lab
###Script will implement TRF to run on all files in a particular directory. You need to specify the directory.
###The files that you should run this on are .fa files. This allows TRF to identify TRF in .fa files
###The output of this program will be a series of files with the identification coordinates of where the TR are located in the genome

for file in /net/beliveau/vol1/reference/Assemblies/hg38_arvispatch/*_alt.fa ;do
	basename=`basename $file`
	echo $basename
	trf409.legacylinux64 $file 2 7 7 80 10 50 500 -f -d &
done

wait
echo "Finished"
