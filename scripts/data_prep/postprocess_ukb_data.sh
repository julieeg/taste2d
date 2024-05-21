#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#

#$ -j y
#$ -cwd
#$ -o ./reports


ANC=$1


source /broad/software/scripts/useuse
use R-4.1


Rscript ../scripts/data_prep/postprocess_ukb_data.R ${ANC}

#EOF

