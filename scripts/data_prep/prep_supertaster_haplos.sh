#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=0:30:00
#$ -o reports/

#$ -j y
#$ -cwd



scratch=/broad/hptmp/gervis
plink2=../opt/plink2

source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen



## extract snp dosage

echo "rs713598
rs1726866
rs10246939" > supertaster_snps

bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr7_v3.bgen -incl-rsids supertaster_snps > ${scratch}/supertaster.bgen

${plink2} --bgen ${scratch}/supertaster.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--extract supertaster_snps \
--export A \
--memory 30000 \
--out ../data/processed/supertasters/supertasters

mv ../data/processed/supertasters/supertasters.log logs/plink2_supertasters_dose.log



## Prep for haplotype phasing

# make genetic map (use GRCh37 build)
wget -O genetic_map_GRCh37 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140404_data_for_phase3_paper/shapeit2_scaffolds/genetic_map_chr7_combined_b37.20140701.txt 

cat /broad/ukbb/imputed/ukb_mfi_chr7_v2.txt | grep -f supertaster_snps | awk '{print $2}' > supertaster_pos 
cat genetic_map_GRCh37 | grep -f supertaster_pos > supertaster_genetic_map


#EOF

