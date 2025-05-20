#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=6:00:00

#$ -j y
#$ -cwd
#$ -o ./reports




plink2=../opt/plink2

geno_dir=../data/processed/supertasters
chr7_phased=${geno_dir}/chr7_2mb_phased.vcf.gz

chr7_bgen=/broad/ukbb/imputed_v3/ukb_imp_chr7_v3.bgen
scratch=/broad/hptmp/gervis


source /broad/software/scripts/useuse
use R-4.1

reuse -q Anaconda3
source activate ../../opt/bgen



## Extract phased supertaster alleles from vcf.gz
echo "" > ${geno_dir}/supertaster_phased.vcf
zcat ${chr7_phased} | head -10 > ${geno_dir}/supertaster_phased.vcf
zcat ${chr7_phased} | grep -E 'rs1726866|rs10246939|rs713598' >> ${geno_dir}/supertaster_phased.vcf



## Configure haplotypes in R
R --no-save <<EOF
library(tidyverse); library(data.table); library(vcfR)
vcf <- read.vcfR(paste0("${geno_dir}/supertaster_phased.vcf"))
as.data.frame(t(extract.haps(vcf, mask=F, unphased_as_NA=TRUE, verbose=T))) %>%
 mutate(parent=gsub(".*_", "", rownames(.)), id=gsub("_.*", "", rownames(.))) %>% 
 pivot_wider(values_from = c(rs713598, rs1726866, rs10246939),  names_from = parent, names_sep="_") %>% 
 mutate(haplo_0 = paste0(rs713598_0, rs1726866_0, rs10246939_0), 
        haplo_1 = paste0(rs713598_1, rs1726866_1, rs10246939_1)) %>%
  mutate(diplo = paste0(haplo_0, "/", haplo_1)) %>%
  fwrite(paste0("${geno_dir}/supertaster_phased_haps.csv"), row.names=F, col.names=T, quote=F, na=NA)
EOF



## Extract dosage for supertaster alleles
bgenix -g $chr7_bgen -incl-rsids ./supertaster_snps > ${geno_dir}/supertaster.bgen
${plink2} \
--bgen ${geno_dir}/supertaster.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--export A \
--out ${geno_dir}/supertaster



## Run Rscript to compile ukb analysis dataset
Rscript ../scripts/data_prep/postprocess_ukb_data.R



#EOF

