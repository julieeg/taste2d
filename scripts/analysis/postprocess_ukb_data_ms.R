# Post-processing script for UKB data (phenotypes + genotypes) 
# Last Updated: 09-24-2024
# Author: Julie E. Gervis


################################################################################
## Set Up & Load Required Libraries 
################################################################################

## Load packages
library(tidyverse) ; library(data.table) ; library(vcfR)
#installed_packages <- c("table1", "paletteer", "circlize")


## Load basic functions from pantry GitHub repo
#system("git clone https://github.com/julieeg/pantry.git")
invisible(lapply(list.files("../../pantry/functions/", full.names = T), function(f) {
  suppressMessages(suppressWarnings(source(f, echo = FALSE, verbose = FALSE)))
})) #source("../scripts/pantry/pantry.R")


# ===================================
##  Define variable vectors  
# ===================================

# Genetic PCs ----------
gPC_vars <- c(paste0("gPC", 1:10))

## Basic phenotypes  ----------
pheno_vars <- c("unrelated", "ancestry", "t2d_case.f", "t2d_med_any", "med_mets", "med_code", "fasting_hrs",
                "ac", "ac_date", "age", "sex", gPC_vars, "smoke.lab", "alch_freq.lab", "ALC", 
                "physact_level", "physact_met_excess", "income", "educ_level.lab", "educ_isced.lab", "educ_years",
                "glu", "hba1c_max", "bmi", "sbp", "dbp", "crp", "chol", "tg", "ldl" , "hdl")

# FFQ variables ----------
ffq_vars <- c("cooked_veg", "raw_veg", "fresh_fruit", "dried_fruit", 
              "bread_intake", "bread_type_white_vs_brown_or_whole", "bread_type.lab",
              "cereal_intake", "cereal_type_sugar_vs_any_bran", "cereal_type.lab", 
              "cheese", "milk_type_full_vs_low_or_nonfat", "milk_type.lab", 
              "spread_type_butter_vs_any_other", "oily_fish", "nonoily_fish", "procmeat", "poultry", 
              "beef", "lamb", "pork", "tea", "coffee", "coffee_type_decaf_vs_regular",
              "water", "addsalt.lab", "addsalt_freq_QT", "addsalt_always_often_vs_nrs", 
              "hotdrink_temp", "hotdrink_temp_hot_or_vhot_vs_warm")



################################################################################
## Build primary genetic exposures
################################################################################

# ============================================
## Compile TAS2R38 haplotypes  & diplotypes
# ============================================

paste("Compiling supertaster haplotypes & SNP dosage")

## Compile haplotype data --------------------

haplos_id <- fread("../data/processed/supertasters/supertaster_phased_haps.csv")

haplos_id <- haplos_id %>%
  mutate(taster_status = case_when(
    haplo_0 == "GGC" & haplo_1 == "GGC" ~ "Supertaster",
    haplo_0 == "GGC" & haplo_1 == "CAT" | haplo_1 == "GGC" & haplo_0 == "CAT" ~ "Taster",
    haplo_0 == "CAT" & haplo_1 == "CAT" ~ "Nontaster",
    TRUE ~ as.character(NA))) %>% 
  mutate(taster_status = ifelse(is.na(taster_status), "Other", taster_status)) %>%
  mutate(taster_status = factor(taster_status, levels = c("Nontaster", "Taster", "Supertaster", "Other")))


## Merge haplotype & dosage data --------------------

genos_id <- haplos_id %>%
  select(c("id", ends_with("_0"), ends_with("_1"), "diplo", "taster_status")) %>%
  left_join(fread("../data/processed/supertasters/supertasters.raw") %>%
              select(id=IID, starts_with("rs")), by = "id") %>%
  # recode to match taster allele 
  mutate("rs713598_G"=2-rs713598_C, 
         "rs10246939_C"=2-rs10246939_T) %>%
  select(-c("rs713598_C", "rs10246939_T"))


## Clean diplotype variables --------------------

## Recode diplotypes so palindromic heterozygotes are coded the same
## example: CAT/GGC AND GGC/CAT are both coded as CAT/GGC

recode_diplos <- function(diplos) {
  haplo_pairs <- strsplit(diplos, "/")

  # Sort each pair and then combine them back into a string
  recoded <- sapply(haplo_pairs, function(pair) {
    sorted_pair <- sort(pair)
    paste(sorted_pair, collapse = "/")
  }) ; return(recoded) 
}

genos_id <-  genos_id %>% 
  mutate(diplos_pal=recode_diplos(diplo)) 


## Make variable for common diplotypes --------------------

diplos_common <- which(prop.table(table(genos_id$diplos_pal))>0.001)
genos_id <- genos_id %>% mutate(
  diplos_common = ifelse(diplos_pal %in% names(diplos_common), diplos_pal, NA))


## Recode alleles as amino acids  --------------------

genos_id <- genos_id %>% 
  mutate(diplos_common_AA=gsub("CAT", "AVI", gsub("CAC", "AVV", gsub("CGT", "AAI", 
                               gsub("CGC", "AAV", 
                                     gsub("GGC", "PAV", gsub("GGT", "PAI", 
                                          gsub("GAT", "PVI", gsub("GAC", "PVV", 
                                                diplos_common)))))))))


## Create variable for canonical diplotypes (AVI/AVI, AVI/PAV, PAV/PAV)

taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")

genos_id <- genos_id %>%
  mutate(taste_diplos = factor(ifelse(diplos_common_AA %in% taste_diplos, diplos_common_AA, NA),
                               levels=taste_diplos)) %>%
  mutate(
    taste_diplos.num = case_when(
      taste_diplos == "AVI/AVI" ~ 0,
      taste_diplos == "AVI/PAV" ~ 1,
      taste_diplos == "PAV/PAV" ~ 2,
      TRUE ~ as.numeric(NA)),
    taste_alleles = (rs713598_G+rs1726866_G+rs10246939_C)
  )
  

## Add haplo_0_aa and haplo_1_aa to analysis.l dataframe list
genos_id <- genos_id %>% left_join(geno_to_aa.fun(genos_id), by = "id")


paste("Done compiling genotype data! ")
head(genos_id)



################################################################################
## Run post-processing for ukb phenotype data
################################################################################

paste("Post-processing of UKB phenotypes ... ")

phenos_id <- phenos %>%   
select("id", all_of(pheno_vars), all_of(gPC_vars), all_of(nutrient_vars), all_of(ffq_vars)) %>%
  mutate(
    sex = case_when(sex == 1 ~ "Male", sex == 0 ~ "Female"),
    income_level.lab=case_when(
      income == "Prefer not to answer" ~ as.character(NA),
      income == "Do not know" ~ as.character(NA),
      income != "Prefer not to answer" & income != "Do not know" ~ as.character(income),
      TRUE ~ as.character(NA)),
    smoke_level.lab = case_when(
      smoke.lab == "No answer" ~ as.character(NA),
      smoke.lab != "No answer" ~ as.character(smoke.lab),
      TRUE ~ as.character(NA)),
    alch_freq.lab = case_when(
      alch_freq.lab == "Prefer not to answer" ~ as.character(NA),
      alch_freq.lab != "Prefer not to answer" ~ as.character(alch_freq.lab),
      TRUE ~ as.character(NA)),
    alch_heavydrinker = as.factor(ifelse(alch_freq.lab == "Daily or almost daily" | alch_freq.lab == "3-4 per week", 1, 0)),
    alch_lightdrinker = as.factor(ifelse(alch_freq.lab == "Special occasions only" | alch_freq.lab == "Never", 1, 0)),
    #Ref: https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=100026
    kcal_plaus = ifelse(sex == "Male" & TCALS >= 600 & TCALS <= 4800 | 
                          sex == "Female" & TCALS >= 600 & TCALS <= 4300, TCALS, NA),
    fasting_hrs_lt24 = ifelse(fasting_hrs <= 24, fasting_hrs, NA)) %>%
  mutate(
    income_level.lab = factor(income_level.lab, levels = c("Less than 18,000", "18,000 to 30,999",
                                                           "31,000 to 51,999", "52,000 to 100,000",  
                                                           "Greater than 100,000")),
    income_level = factor(income_level.lab, 
                              levels = c("Less than 18,000", "18,000 to 30,999",
                                         "31,000 to 51,999", "52,000 to 100,000",  
                                         "Greater than 100,000"),
                              labels = c("lt_18000", "from_18000_to_30999",
                                         "from_31000_to_51999", "from_52000_to_100000",  
                                         "gt_100000")),
    #Edu levels & yrs based on: https://www.nature.com/articles/s41380-019-0596-9#MOESM1)
    educ_isced.lab = factor(educ_isced.lab, levels=c(paste0("Level ", 1:5))),
    smoke_level.lab = factor(smoke_level.lab, levels = c("Never", "Previous", "Current")),
    alch_freq.lab = factor(alch_freq.lab, ordered = T,
                           levels = c("Daily or almost daily",
                                      "3-4 per week", "1-2 per week", "1-3 per month",
                                      "Special occasions only", "Never")),
    physact_level = factor(physact_level, ordered = T,
                               levels = c("Low", "Moderate", "High")),
    tg_log = log(tg),
    hba1c=hba1c_max,
    fast_cat = case_when(fasting_hrs_lt24 >= 0 & fasting_hrs_lt24 <= 2 ~ "0to2hr",
                         fasting_hrs_lt24 == 3 ~ "3hr",
                         fasting_hrs_lt24 == 4 ~ "4hr",
                         fasting_hrs_lt24 == 5 ~ "5hr",
                         fasting_hrs_lt24 >= 6 ~ "6+hr",
                         TRUE ~ as.character(NA))
  ) %>%
  ## Added descriptive variables for Table 1
  mutate(
    alch_1to4d.lab = ifelse(alch_freq.lab == "3-4 per week" | alch_freq.lab == "1-2 per week", "1-4 per week", "other"),
    addsalt_3lvl.lab = factor(ifelse(addsalt.lab == "Always" | addsalt.lab == "Often", "Always/Often", addsalt.lab),
                          levels=c("Always/Often", "Sometimes", "Never/Rarely")),
    alch_heavydrinker = as.factor(ifelse(alch_freq.lab == "Daily or almost daily" | alch_freq.lab == "3-4 per week", 1, 0)),
    alch_lightdrinker = as.factor(ifelse(alch_freq.lab == "Special occasions only" | alch_freq.lab == "Never", 1, 0))
    ) %>%
  #mutate(ALC_drinkers = ifelse(alch_freq.lab == "Never", NA, ALC)) %>%
  mutate(
    N_complete_kcal = ifelse(is.na(TCALS)==T, "Missing", "Complete"),
    plausible_kcal.f = ifelse(is.na(kcal_plaus)==F, "Plausible", "Implausible_Missing")
  ) %>%
  
  #add variable for 2-hr glucoster
  mutate(glu2hr=ifelse(fast_cat=="0to2hr",glu,NA)) %>%
  
  # Add glucose in mg/dL
  mutate(
    glu_mgdl = glu*18.018,
    glu2hr_mgdl = glu2hr*18.018
  )


  
# Post-process FFQ vars -------------

ffq_id <- phenos %>% select(
  id,
  cooked_veg_QT=cooked_veg, raw_veg_QT=raw_veg,
  fresh_fruit_QT=fresh_fruit, dried_fruit_QT=dried_fruit, 
  oily_fish_QT=oily_fish, nonoily_fish_QT=nonoily_fish,
  procmeat_QT=procmeat, poultry_QT=poultry, cheese_QT=cheese,
  beef_QT=beef, lamb_QT=lamb, pork_QT=pork, 
  bread_type_white_vs_brown_or_whole_BIN=bread_type_white_vs_brown_or_whole, 
  bread_intake_QT=bread_intake,
  milk_type_full_vs_low_or_nonfat_BIN=milk_type_full_vs_low_or_nonfat,
  cereal_type_sugar_vs_any_bran_BIN=cereal_type_sugar_vs_any_bran, cereal_intake_QT=cereal_intake,
  spread_type_butter_vs_any_other_BIN=spread_type_butter_vs_any_other,
  coffee_type_decaf_vs_regular_BIN=coffee_type_decaf_vs_regular, coffee_QT=coffee,
  tea_QT=tea, water_QT=water, 
  addsalt_always_often_vs_nrs_BIN=addsalt_always_often_vs_nrs,
  hotdrink_temp_hot_or_vhot_vs_warm_BIN=hotdrink_temp_hot_or_vhot_vs_warm) %>%
  
  # Replace missing values with medians
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x))  %>% 
  
  # Winsorize data to 5 SD
  mutate(across(where(is.numeric), function(i) winsorize(i, SDs=5))) %>%
  filter(complete.cases(.)==T)

cat("Dimensions of dietary data \n")
dim(ffq_postprocess)


print("Done compiling UKB phenotypes!")


#######################################
#### Merge data into ukb_processed ####
#######################################

phenos_processed_id <- phenos_id %>%
  left_join(genos_id, by = "id") %>%
  left_join(ffq_id, by = "id")



################################################################################
## Apply exclusion criteria to create analysis dataset
################################################################################

### Exclusion criteria
# - No diabetes (Eastwood algorithm + HbA1c <5.7%)
# - EUR genetic ancestry
# - Canonical TAS2R38 diplotype (and complete genetic data)
# - Complete data for random glucose & covariates  
#   - age, sex, fasting time, bmi, smoking, physical activity, alcohol, FFQ (plausible)
#   - sensitivity: Plausible total energy intakes (600-4300 [F] or 4800 [M] kcal/d)


# Function to identify missingness
find_complete_group.fun <- function(x, varname, data=phenos_processed_id) {
  return((data %>% select(c("id", all_of(x))) %>% mutate(complete = ifelse(complete.cases(.), 1, 0)))$complete) }

# add indicator variables for exclusion criteria
phenos_processed_id <- phenos_processed_id %>% 
  mutate(
    n_geno = ifelse(is.na(haplo_0) == F & is.na(haplo_1) == F, 1, 0),
    n_t2d_data = recode_na_as.fun(t2d_case.f),
    n_t2d_contrl = ifelse(!is.na(t2d_case.f) & t2d_case.f == "Control", 1, 0),
    n_glu = recode_na_as.fun(glu, 0),
    n_glu_5sd = find_outliers.fun(glu, SDs=5, recode_outliers = 0),
    n_diplo_taste = ifelse(is.na(taster_status) == T | taster_status == "Other", 0, 1),
    n_diplo_005 = recode_na_as.fun(diplos_common),
    n_fast = recode_na_as.fun(fasting_hrs),
    n_fast_24hr = ifelse(fasting_hrs >24, 0, 1),
    n_24hr = recode_na_as.fun(TCALS),
    n_24hr_plaus = recode_na_as.fun(kcal_plaus)
    ) %>%  mutate(n_covars = find_complete_group.fun(
      c("age", "sex", paste0("gPC", 1:10), "bmi", 
        "smoke_level.lab","physact_level", "alch_freq.lab", "raw_veg_QT"))) %>%
  mutate(n_complete=ifelse(n_t2d_data==1 & n_geno==1 & n_glu==1 & n_fast == 1 & n_covars==1,1,0)) %>%
  mutate(N_ancestry = ifelse(ancestry == "EUR", 1, 0)) %>%
  mutate(N_geno_contrl = ifelse(n_geno == 1 & n_t2d_contrl == 1, 1, 0),
         N_geno_contrl_diplo005 = ifelse(n_geno == 1 & n_t2d_contrl == 1 & n_diplo_005 == 1,1,0),
         N_geno_contrl_taste = ifelse(n_geno == 1 & n_t2d_contrl == 1 & n_diplo_taste == 1,1,0),
         N_geno_contrl_taste_compl = ifelse(n_geno == 1 & n_t2d_contrl == 1 & n_diplo_taste == 1 & n_covars==1,1,0),
         N_geno_contrl_taste_compl_24hrplaus = ifelse(n_geno==1 & n_t2d_contrl == 1 & n_diplo_taste == 1 & n_covars==1 & n_24hr_plaus == 1,1,0))


# ===============================================
## Filter out outliers for glucose levels (>5 SD)
# ===============================================

analysis <- phenos_processed_id %>%
  filter(n_geno==1 & n_t2d_contrl==1 & n_complete==1 & n_fast_24hr==1) %>%
  filter(find_outliers.fun(glu, SD=5)==0)


# ===============================================
## Compile & save ukb_analysis_EUR.rda
# ===============================================


#save postprocessed data as .rda
analysis %>% 
  saveRDS(paste0("../data/processed/ukb_analysis_vManuscript_EUR.rda"))


cat("Done compiling ukb_analysis_ANC.rda dataset!")

#EOF

