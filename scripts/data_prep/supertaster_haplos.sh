#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=10:00:00

#$ -j y
#$ -cwd
#$ -N supertaster_haplos

#$ -o ./reports/


i=$SGE_TASK_ID


scratch=/broad/hptmp/gervis
opt=../../opt

plink2=${opt}/plink2  
shapeit=${opt}/shapeit



source /broad/software/scripts/useuse
use UGER
use R-4.1



# make bfile per chunk

${plink2} \
--bgen ${scratch}/supertaster.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--keep chunks/chunk${i} \
--make-bfile \
--memory 5000 \
--out ../data/temp/supertaster_chunk${i}

mv ../data/temp/supertaster_chunk${i}.log ./logs/plink2_supertaster_chunk${i}.log



# Run shapeit to phase haplotypes 

${shapeit} \
--input-bed ../data/temp/supertaster_chunk${i} \
--input-map ./supertaster_genetic_map \
--output-max ../data/temp/supertaster_chunk${i} \
--output-log ./logs/shapeit_supertaster_haplo.log \
--thread 8



# Prepare vcf with 0|1 notation

${shapeit} -convert \
--input-haps ../data/temp/supertaster_chunk${i}.haps ../data/temp/supertaster_chunk${i}.sample \
--output-vcf ../data/processed/genos/supertasters/supertaster_chunk${i}.vcf \
--output-log ./logs/shapeit_supertaster_haplos2vcf.log \
--thread 8



# Clean data/temp folder
rm ../data/temp/supertaster_chunk${i}*


#EOF

