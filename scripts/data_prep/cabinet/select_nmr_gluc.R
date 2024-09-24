
library(tidyverse)
library(data.table)

# load nmr files
files <- list.files("../data/processed/nmr")
nmr.l <- lapply(files, function(f) fread(paste0("../data/processed/nmr/"), f))


# select glucose fields from each file for download
nmr_glucose.l <- list()

nmr_glucose.l[[1]] <- list(
  nmr_biomarker_qc_flags = nmr.l[[2]] %>%
    select(id=eid, visit_index, Glucose),
  nmr_processed = nmr.l[[3]] %>%
    select(id=eid, visit_index, Glucose),
  sample_qc_flags = nmr.l[[4]] %>%
    rename(id=eid),
  outlier_plate_info = nmr.l[[5]] %>%
    filter(Biomarker == "Glucose")
) 

saveRDS(nmr_glucose.l, "../data/processed/nmr/ukb_nmr_glucose.rda")
