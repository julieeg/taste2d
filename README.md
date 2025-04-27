# taste2d
# Chef File 
# Last Update: 04-26-2025


#### Prepare ukb phenotypes & genotypes 
qsub ../scripts/data_prep/prep_ukb_phenos.sh
qsub ../scripts/data_prep/prep_phasing.sh

#### Run postprocessing to prepare analytical dataset
qsub ../scripts/data_prep/postprocess_ukb_data.sh 

#sync data to local R for analysis
rsync -aP uger:$florezlab_taste2d/data/processed/"ukb_analysis_EUR.rda" ../data/processed/
rsync -aP uger:$florezlab_taste2d/data/processed/diet/"ukb_EUR_dietPCs.rda" ../data/processed/diet/

#primary association analyses 
Rscript --vanilla ../scripts/analysis/run_lm.R EUR; done 

#secondary analysis with other bitter snps
qsub ../scripts/get_bitter_snps.sh
rsync uger:$florezlab_taste2d/data/processed/bitter_snps.raw ../data/procesed/bitter_snps.raw

