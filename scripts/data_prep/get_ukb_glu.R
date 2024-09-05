

# load packages
library(tidyverse)
library(data.table)


# load basic functions for data preparation & cleaning
source("../scripts/pantry/pantry.R")


## Adding variables for glucose correction factors; NMR glucose & glucose at difference time.points


######################################
##  Glucose  ##
######################################

print("Preparing glucose biomarker data  ...")


## Biochemical measures ===============

#assay results ------
gluB_fields <- c(gluB = 30740)

gluB_vars <- setNames(
  paste0("f.", gluB_fields, c(".0.0", ".1.0")), paste0(names(gluB_fields), c(".0", ".1"))
  ) ; gluB_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz",
                         data.table=TRUE, stringsAsFactors = FALSE) %>%
      select(id = f.eid, all_of(gluB_vars))

    
# assay quality ------
gluB_assay_fields <- c(gluB_assay_date = 30741, gluB_aliquot = 30742, gluB_correction_type = 30743, 
  gluB_correction_reason = 30744, gluB_missing = 30745, gluB_reportable = 30746,
  dilution = 30897) ; gluB_assay_vars <- setNames(
    paste0("f.", rep(gluB_assay_fields, each=2), c(".0.0", ".1.0")), 
    paste0(rep(names(gluB_assay_fields), each=2), c(".0", ".1"))
    ) ; gluB_assay_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb46203.tab.gz",
                    data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(gluB_assay_vars))


## Fasting time at blood draw ============
fast_field <- c(fasting_hrs=74); fast_vars <- setNames(
  paste0("f.", fast_field, c(".0.0", ".1.0")), paste0(names(fast_field), c(".0", ".1"))
  ) ; fast_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                         data.table=FALSE, stringsAsFactors=FALSE) %>% 
  select(id=f.eid, all_of(fast_vars))


## NMR glucose =========================

gluM_fields <- c(gluM = 23470) ; gluM_vars <- setNames(
  paste0("f.", gluM_fields, c(".0.0", ".1.0")), paste0(names(gluM_fields), c(".0", ".1"))
  ) ; gluM_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb47805.tab.gz",
                        data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(gluM_vars))



#############################################
## Recode variables as factors & labels
#############################################

aliquot_labs <- list("Manual" = 0, "Aliquot 1" = 1, "Aliquot 2" = 2, "Aliquot 3" = 3, "Aliquot 4" = 4)
correction_type_labs <- list("None" = 0, "Date & aliquot" = 1, "Date only" = 2, "No data" = NA)
correction_reason_labs <- list("Normal for assay type" = 0, "No tip info" = 1, 
                               "No aliquot info" = 2, "Original assay >15 SDs below mean of manual assays" = 3, 
                               "Resurvey only:Insufficient N for aliquot dispensing machine" = 4)
missing_labs <- list("No data" = 1, "Value outside reportable range" = 2, "Aliquot problem - tip" = 3, 
                     "Aliquot problem - missclassifiation" = 4, "Aliquot 4 used" = 5, 
                     "Analyzer deemed result not reportable" = 7, "Sample dilution factor b/w 1-10%" = 8, 
                     "Sample dilution factor >10%")
reportable_labs <- list("Reportable" = 1, "Not reportable: too low after corrections" = 2, 
                        "Not reportable: too high after corrections" = 3, 
                        "Not reportable: original assay too low" = 4, 
                        "Not reportable: original assay too high" = 5)

gluB_assay_id <- gluB_assay_id %>% mutate(
  #across(contains("date"), ~ as.Date.numeric(.)),
  across(contains("aliquot"), ~as.factor(.)),
  across(contains(c("correction", "missing", "reportable")), ~as.factor(.))
) %>% mutate(
  gluB_aliquot.0.lab = descr_label_ordered.fun(., "gluB_aliquot.0", aliquot_labs),
  gluB_aliquot.1.lab = descr_label_ordered.fun(., "gluB_aliquot.1", aliquot_labs),
  gluB_correction_type.0.lab = descr_label_ordered.fun(., "gluB_correction_type.0", correction_type_labs),
  gluB_correction_type.1.lab = descr_label_ordered.fun(., "gluB_correction_type.1", correction_type_labs),
  gluB_correction_reason.0.lab = descr_label_ordered.fun(., "gluB_correction_reason.0", correction_reason_labs),
  gluB_correction_reason.1.lab = descr_label_ordered.fun(., "gluB_correction_reason.1", correction_reason_labs),
  gluB_missing.0.lab = descr_label_ordered.fun(., "gluB_missing.0", missing_labs),
  gluB_missing.1.lab = descr_label_ordered.fun(., "gluB_missing.1", missing_labs),
  gluB_reportable.0.lab = descr_label_ordered.fun(., "gluB_reportable.0", reportable_labs),
  gluB_reportable.1.lab = descr_label_ordered.fun(., "gluB_reportable.1", reportable_labs)
)


######################################
##  Merge & save datasets
######################################

## merge & combine with UKB_analysis dataset 
glu_merged_id <- gluB_id %>%
  full_join(gluB_assay_id, by = "id") %>%
  full_join(fast_id, by = "id") %>%
  full_join(gluM_id, by = "id")


ukb_index <- readRDS("../data/processed/ukb_analysis_ALL.rda") %>%
  select(id, unrelated, ancestry)

## Merge & save 
ukb_index %>% 
  left_join(glu_merged_id, by = "id") %>%
  filter(unrelated == T) %>%
  write_csv("../data/processed/ukb_glucose_assays.csv")


## Merge & save 
ukb_index %>% 
  left_join(glu_merged_id, by = "id") %>%
  filter(unrelated == T) %>%
  saveRDS("../data/processed/ukb_glucose_assays.rda")

## END_OF_SCRIPT


