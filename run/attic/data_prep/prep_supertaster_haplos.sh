#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=0:30:00
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

supertaster_dir=../data/processed/supertasters



## create bgen file for 100m window on chr7

reuse -q Anaconda3
source activate ../../opt/bgen

cat ${bgen_dir}/ukb_mfi_chr7_v3.txt | awk '{ if($3 > 100000000 && $3 <150000000) {print $2} }' > chr7_b37_50mb


bgenix -g ${chr7_bgen} -i ${chr7_bgen}.bgi -vcf -incl-rsids chr7_b37_50mb | bgzip -c > ${supertasters_dir}/chr7_50mb_unphased.vcf.gz




## run phasing with shapeit4

# get genetic map for b37
#wget -O chr7_b37_gmap ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140404_data_for_phase3_paper/shapeit2_scaffolds/genetic_map_chr7_combined_b37.20140701.txt 


# index the vcf.gz input file

bcftools index ${supertasters_dir}/chr7_50mb_unphased.vcf.gz -o ${supertasters_dir}/chr7_50mb_unphased.vcf.gz


shapeit4.2 --input ${supertasters_dir}/chr7_50mb_unphased.vcf.gz \
--map chr7_b37_gmap \
--region 7:100000000-150000000 \
--output ${supertaster}/chr7_50mb_phased.vcf.gz \
--thread 8

#&& rm ../data/processed/supertasters/chunk${i}.b* \
#&& rm ../data/processed/supertasters/chunk${i}.fam





#EOF
