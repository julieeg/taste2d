#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=10:00:00
#$ -o reports/

#$ -j y
#$ -cwd
#$ -N supertaster_haplos


i=$SGE_TASK_ID


scratch=/broad/hptmp/gervis

plink2=../opt/plink2  
shapeit=../opt/shapeit


source /broad/software/scripts/useuse
use R-4.1


# chunk iids

awk '{ {FS = ","} if(NR > 1) {print $1, $1, 0, $4} }' ../data/processed/ukb_phenos_unrelated_ALL.csv > iids
mkdir -p chunks ; split -n 10 --numeric-suffixes=10 -a 2 iids chunks/chunk 


# make bfile per chunk
${plink2} \
--bgen ${scratch}/supertaster.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--keep chunks/chunk${i} \
--make-bfile \
--memory 10000 \
--out ../data/processed/supertasters/chunk${i}

mv ../data/processed/supertasters/chunk${i}.log logs/plink2_supertaster_chunk${i}.log



# Run shapeit to phase haplotypes 

${shapeit} \
--input-bed ../data/processed/supertasters/chunk${i} \
--input-map supertaster_genetic_map \
--output-max ../data/processed/supertasters/chunk${i} \
--output-log logs/shapeit_supertaster_haplo.log \
--thread 8 \
&& rm ../data/processed/supertasters/chunk${i}.b* \
&& rm ../data/processed/supertasters/chunk${i}.fam



# Prepare vcf with 0|1 notation

${shapeit} -convert \
--input-haps ../data/processed/supertasters/chunk${i}.haps ../data/processed/supertasters/chunk${i}.sample \
--output-vcf ../data/processed/supertasters/chunk${i}.vcf \
--output-log logs/shapeit_supertaster_haplos2vcf.log \
--thread 8 \
&& rm ../data/processed/supertasters/chunk${i}.sample \



# Convert to haplotype notation using R
R --vanilla <<EOF
library(tidyverse); library(data.table); library(vcfR)
vcf <- read.vcfR(paste0("../data/processed/supertasters/chunk${i}.vcf"))
hap <- as.data.frame(t(extract.haps(vcf, mask=F, unphased_as_NA=TRUE, verbose=T))) %>%
 mutate(parent=gsub(".*_", "", rownames(.)), id=gsub("_.*", "", rownames(.))) %>%
 pivot_wider(values_from = c(rs713598, rs1726866, rs10246939),  names_from = parent, names_sep="_") %>%
 mutate(haplo_0 = paste0(rs713598_0, rs1726866_0, rs10246939_0), 
        haplo_1 = paste0(rs713598_1, rs1726866_1, rs10246939_1)) %>%
  mutate(diplo = paste0(haplo_0, "/", haplo_1)) %>%
  fwrite(paste0("../data/processed/supertasters/chunk${i}.csv"), row.names=F, col.names=T, quote=F, na=NA)
EOF



##EOF

