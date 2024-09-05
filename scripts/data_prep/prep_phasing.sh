#!/bin/bash


#$ -l h_vmem=560G
#$ -l h_rt=6:00:00
#

#$ -j y
#$ -cwd
#$ -o ./reports


source /broad/software/scripts/useuse
use R-4.1

geno_dir=../data/processed/supertasters



# Extract phased supertaster alleles from vcf.gz
echo "" > ${geno_dir}/supertaster_phased.vcf
zcat ${geno_dir}/chr7_2mb_phased.vcf.gz | head -10 > ${geno_dir}/supertaster_phased.vcf
zcat ${geno_dir}/chr7_2mb_phased.vcf.gz | grep -E 'rs1726866|rs10246939|rs713598' >> ${geno_dir}/supertaster_phased.vcf








phased <- read.vcfR("../data/processed/supertasters/supertaster_phased.vcf") 
haps <- extract.haps(phased, unphased_as_NA = TRUE, verbose = TRUE)

fwrite("../data/processed/supertasters/chr7_2mg_phased_haps.txt")

EOF


#EOF - TEMPORARY
