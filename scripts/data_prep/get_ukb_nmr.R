## Download & prepare NMR data using ukbnmr package

# load packages
library(tidyverse)
library(data.table)

library(ukbnmr)



ukb48298 <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb48298.tab.gz",
                  data.table=TRUE, stringsAsFactors = FALSE) ; head(ukb48298)

nmr <- remove_technical_variation(ukb48298)

# save files as .csvs
fwrite(nmr$biomarkers, file="../data/processed/nmr/ukb_nmr_processed.csv")
fwrite(nmr$biomarker_qc_flags, file="../data/processed/nmr/ukb_nmr_biomarker_qc_flags.csv")
fwrite(nmr$sample_processing, file="../data/processed/nmr/ukb_nmr_sample_qc_flags.csv")
fwrite(nmr$log_offset, file="../data/processed/nmr/ukb_nmr_biomarker_log_offset.csv")
fwrite(nmr$outlier_plate_detection, file="../data/processed/nmr/ukb_outlier_plate_info.csv")


##EOF