#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=40.0G

jellyfish count -s 3300M -m 36 -o hg38_alt_patch_36.jf hg38_alt_patch.fa 
