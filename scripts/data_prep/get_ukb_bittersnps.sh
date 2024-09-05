#!/bin/bash
#$ -l h_vmem=30G#$ -l h_rt=3:00:00#$ -j y#$ -cwd#$ -o ./reports
ukb_bgen_dir=/broad/ukbb/imputed_v3ukb_sample_dir=/humgen/florezlab/UKBB_app27892 #[ukb27892_imp_chrAUT_v3_s487395.sample]scratch=/broad/hptmp/gervissource /broad/software/scripts/useusereuse -q Anaconda3source activate ../../opt/bgen## Grab SNPs: rs1726866 (PROP, chr7); rs10772420 (Quinine; chr12); rs2597979 (caffeine; chr12)

for i in 7 12; do
bgenix -g ${ukb_bgen_dir}/ukb_imp_chr${i}_v3.bgen -incl-rsids bitter_snps > ${scratch}/bitter_snps_chr${i}.bgen
done ; cat-bgen -g ${scratch}/bitter_snps*.bgen -og ../data/processed/bitter_snps.bgen -clobber \
&& rm ${scratch}/bitter_snps_chr*.bgen


## Export SNP dosage
../../opt/plink2 \--bgen ../data/processed/bitter_snps.bgen ref-first \--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \--export A \
--out ../data/processed/bitter_snps \
--memory 2000
&& rm ../data/processed/bitter_snps_chr*



# Covert .raw to .csv using R
R --vanilla <<EOF library(tidyverse) ; library(data.table) 
fread("../data/processed/bitter_snps.raw") %>% select(IID, starts_with("rs")) %>% fwrite("../data/processed/bitter_snps