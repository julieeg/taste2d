#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=0:15:00

#$ -j y
#$ -cwd
#$ -o ./reports


## Run Rscript to compile ukb analysis dataset
Rscript ../scripts/data_prep/postprocess_ukb_data.R



#EOF

