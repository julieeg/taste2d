#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=10:00:00

#$ -j y
#$ -cwd

#$ -o ./reports/


scratch=/broad/hptmp/gervis

opt=../../opt
plink2=${opt}/plink2


source /broad/software/scripts/useuse
use UGER

reuse -q Anaconda3
source activate $opt/bgen




## Extract snp dosage

bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr7_v3.bgen -incl-rsids supertaster_snps > ${scratch}/supertaster.bgen

${plink2} --bgen ${scratch}/supertaster.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--extract supertaster_snps \
--export A \
--memory 30000 \
--out ../data/processed/genos/supertasters/supertasters

mv ../data/processed/genos/supertasters/supertasters.log ./plink2_supertasters_dose.log



## Prep for haplotype phasing

# make genetic map (Build: GRCh37 to match ukb)

wget -O ./genetic_map_GRCh37 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140404_data_for_phase3_paper/shapeit2_scaffolds/genetic_map_chr7_combined_b37.20140701.txt 

cat /broad/ukbb/imputed/ukb_mfi_chr7_v2.txt | grep -f supertaster_snps | awk '{print $2}' > ./supertaster_pos.txt 
cat ./genetic_map_GRCh37 | grep -f ./supertaster_pos.txt > ./supertaster_genetic_map



# chunk iids to run phasing as an array 

awk '{ {FS = ","} if(NR > 1) {print $1, $1, 0, $4} }' ../data/processed/phenos/ukb_phenos_unrelated.csv > iids
mkdir -p chunks ; split -n 30 --numeric-suffixes=10 -a 2 iids chunks/chunk 


# Run phasing as an array

qsub -t 10-29 ../scripts/data_prep/supertaster_haplos.sh



#EOF


