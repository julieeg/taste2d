#!/bin/bash

#$ -l h_vmem=60G
#$ -l h_rt=48:00:00
#$ -o reports/

#$ -j y
#$ -cwd



scratch=/broad/hptmp/gervis
plink2=../opt/plink2


source /broad/software/scripts/useuse
use .shapeit4-4.2.1
use .bcftools-1.19
use .tabix-0.2.6


bgen_dir=/broad/ukbb/imputed_v3
chr7_bgen=${bgen_dir}/ukb_imp_chr7_v3.bgen

# set prefix for chr7 files
chr7=../data/processed/supertasters/chr7


reuse -q Anaconda3
source activate ../../opt/bgen


####################################################
# Run phasing for chr7 on common highqual variants
####################################################

# Build small bgen file with 2MB window around TAS2R38 variants
cat ${bgen_dir}/ukb_mfi_chr7_v3.txt | awk '{ if($3 > 140000000 && $3 < 142000000 && $6 > 0.05 && $7 > 0.4) {print $2} }' > ${chr7}_2mb
bgenix -g ${chr7_bgen} -incl-rsids ${chr7}_2mb > ${chr7}_2mb.bgen
bgenix -index -g ${chr7}_2mb.bgen


# Make list of unrelated samples for phasing
R --no-save <<EOF
library(data.table) ; library(tidyverse)
fread("../data/processed/ukb_phenos.csv") %>% filter(unrelated == T) %>% select('#FID'=id, IID=id) %>% write_tsv("../data/processed/ukb_unrelated.sample")
EOF


# Convert bgen to vcf.gz
../opt/plink2 \
  --bgen ${chr7}_2mb.bgen ref-first \
  --sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
  --keep ../data/processed/ukb_unrelated.sample \
  --export vcf vcf-dosage=DS \
  --out ${chr7}_2mb


# Zip vcf file to vcf.gz & create index file for phasing
bgzip -c ${chr7}_2mb.vcf > ${chr7}_2mb.vcf.gz
bcftools index -c ${chr7}_2mb.vcf.gz


# ===========================
## Run phasing with shapeit
# ===========================

# Get genetic map for b37 (save in run folder, from where you run commands)
wget -O chr7_b37_gmap ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140404_data_for_phase3_paper/shapeit2_scaffolds/genetic_map_chr7_combined_b37.20140701.txt 


shapeit4.2 --input ${chr7}_2mb.vcf.gz \
--map chr7_b37_gmap \
--region 7 \
--output ${chr7}_2mb_phased.vcf.gz \
--thread 10


#EOF

