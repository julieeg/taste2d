## Post-processing


# load packages
library(tidyverse) 
library(data.table)

# load basic functions
source("../scripts/R/basic_functions.R", echo=F)


#####################
## Variable groups ##
#####################

# genetic PCs 
gPC_vars <- c(paste0("gPC", 1:10))

## phenotypes (base/ses/lifestyle/clinical; diet) 
pheno_vars <- c("age", "sex", "income", "educ_level.lab", "smoke.lab", "alch_freq.lab", "ALC",
                "pa_met_excess_lvl", "glu", "bmi", "waist2hip", "sbp", "dbp", 
                "chol", "tg", "ldl" , "hdl", "crp", "ancestry", "fasting_hrs")

ffq_vars <- c("cooked_veg", "raw_veg", "fresh_fruit", "dried_fruit", 
              "bread_intake", "bread_type_white_vs_brown_or_whole", "bread_type.lab",
              "cereal_intake", "cereal_type_sugar_vs_any_bran", "cereal_type.lab", 
              "cheese", "milk_type_full_vs_low_or_nonfat", "milk_type.lab", 
              "spread_type_butter_vs_any_other", "oily_fish", "nonoily_fish", "procmeat", "poultry", 
              "beef", "lamb", "pork", "tea", "coffee_intake", "coffee_type_decaf_vs_regular",
              "water", "addsalt.lab", "addsalt_freq", "addsalt_always_often_vs_nrs", 
              "hotdrink_temp", "hotdrink_temp_hot_or_vhot_vs_warm",
              "dietPC1", "dietPC2", "dietPC3", "dietPC4", "dietPC5")


nutrient_vars <- c("TCALS", "CHO_pct", "FIBER", "FAT_pct", "MUFA_pct", "PUFA_pct", 
                   "SFA_pct", "PRO_pct", "CHO2FIB", "FIB2CHO")



################################################################################
## Build ukb data
################################################################################


## Load & compile ukb genoypte data  -------------------------
geno_pre <- fread("../data/raw/ukb_bitter_haplotype.csv") %>% select(-"V1") %>%
  left_join(fread("../data/raw/supertaster_dose.csv") %>% select(-c("V1", "IID")), by = "id") %>% 
  #mutate(bitter_pts_std = zscore.fun(bitter_pts)) %>%
  mutate(taster_status = factor(taster_status, levels = c("nontaster", "taster", "supertaster", NA)),
         Taster_Status = factor(taster_status, levels = c("nontaster", "taster", "supertaster", NA),
                                labels=c("Nontaster", "Taster", "Supertaster")),
         taste_alleles = ifelse(rs713598_1 == "G", 1, 0) + ifelse(rs713598_0 == "G", 1, 0) + 
           ifelse(rs1726866_1 == "G", 1, 0) + ifelse(rs1726866_0 == "G", 1, 0) + ifelse(rs10246939_1 == "C", 1, 0) + 
           ifelse(rs10246939_0 == "C", 1, 0)) %>%
 
  # Add scaled score
  #mutate(bitter_pts_scale = ((bitter_pts / sum(bitter_sumstats$BETA_STD)) * 14 )) %>%
  
  select("id", "taster_status",  #"bitter_pts", "bitter_pts_std", "bitter_pts_scale",
         "haplo_0", "haplo_1", "diplo", "taste_alleles",
         "rs713598_C", "rs1726866_G", "rs10246939_T")


## Load & compile T2D data ---------------------------------

t2d <- fread("../data/raw/UKB_Diabetes.csv", data.table=FALSE, stringsAsFactors=FALSE) %>% 
  rename(id = f.eid) %>%
  select("id", starts_with(c("probable", "possible")), "agedm_ts_or_ni_all", 
         "meds_any_sr_ni_all", "meds_any_sr_ni_ts_all", "dm_unlikely_all") %>% 
  left_join(fread("../data/raw/UKB_HbA1c.csv") %>% rename(id=f.eid) %>%
              select("id", "hba1c.30750.NGSP.max"),
            by = "id") %>% 
  left_join(fread("../data/raw/UKB_ICD10DM.csv") %>% rename(id=f.eid), 
            by = "id")

t2d_id <- t2d %>% mutate(
  t2d_med_any = case_when(
    meds_any_sr_ni_all == 1 ~ 1, 
    meds_any_sr_ni_ts_all == 1 ~ 1,
    meds_any_sr_ni_all == 0 & meds_any_sr_ni_ts_all == 0 ~ 0
  )) %>% 
  mutate(
    t2d_case = case_when(
      possible_t2dm_all == 1 ~ 1, 
      probable_t2dm_all == 1 ~ 1,
      hba1c.30750.NGSP.max >= 6.5 & possible_t1dm_all != 1 & probable_t1dm_all != 1 & dm_main != 1 & dm_secondary != 1 ~ 1,
      hba1c.30750.NGSP.max < 5.7 & dm_unlikely_all == 1 & dm_main != 1 & dm_secondary != 1 & t2d_med_any != 1 ~ 0,
      TRUE ~ as.numeric(NA)
    )) %>%
  mutate(
    t2d_case.f = case_when(
      as.numeric(t2d_case) == 0 ~ "Control",
      as.numeric(t2d_case) == 1 ~ "Case"),
    t2d_age_diagnosis = agedm_ts_or_ni_all
  ) %>%
  select("id", "t2d_case", "t2d_case.f", "hba1c"="hba1c.30750.NGSP.max", 
         "t2d_med_any", "t2d_age_diagnosis") 


## Load & compile ukb phenotype data  -------------------------

phenos_pre <- fread("../data/raw/ukb_pheno_unrelated_EUR_rmsd4.csv")  %>%
  select("id", all_of(pheno_vars), all_of(gPC_vars), all_of(nutrient_vars)) %>%
  mutate(
    age2 = age*age,
    sex = case_when(sex == 1 ~ "Male", sex == 0 ~ "Female"),
    income_level.lab=case_when(
      income == "Prefer not to answer" ~ as.character(NA),
      income == "Do not know" ~ as.character(NA),
      income != "Prefer not to answer" & income != "Do not know" ~ as.character(income),
      TRUE ~ as.character(NA)),
    educ_level.lab = case_when(
      educ_level.lab == "Prefer not to answer" ~ as.character(NA),
      educ_level.lab != "Prefer not to answer" ~ as.character(educ_level.lab),
      TRUE ~ as.character(NA)),
    smoke_level.lab = case_when(
      smoke.lab == "No answer" ~ as.character(NA),
      smoke.lab != "No answer" ~ as.character(smoke.lab),
      TRUE ~ as.character(NA)),
    alch_freq.lab = case_when(
      alch_freq.lab == "Prefer not to answer" ~ as.character(NA),
      alch_freq.lab != "Prefer not to answer" ~ as.character(alch_freq.lab),
      TRUE ~ as.character(NA))) %>% 
  mutate(
    income_level.lab = factor(income_level.lab, 
                              levels = c("Less than 18,000", "18,000 to 30,999",
                                         "31,000 to 51,999", "52,000 to 100,000",  
                                         "Greater than 100,000"),
                              labels = c("lt_18000", "from_18000_to_30999",
                                         "from_31000_to_51999", "from_52000_to_100000",  
                                         "gt_100000")),
    #Edu levels & yrs based on: https://www.nature.com/articles/s41380-019-0596-9#MOESM1)
    educ_level.lab = factor(educ_level.lab, 
                            levels = c("College or university degree", # ~20yrs 
                                       "NVQ/HND or equivalent", # 2 of 3 years bachelor's degree ~19yrs
                                       "Other professional qualifications", # e.g., nursing degree, teaching degree ~ 15yrs
                                       "A/AS levels or equivalent", # 1 year bachelor's degree ~13yrs
                                       "O/GCSE levels or equivalent", # HS + Associates degree ~10yrs
                                       "CSEs or equivalent", # completed HS ~10yrs
                                       "None of the above")), #~7yrs
    alch_freq.lab = factor(alch_freq.lab, ordered = T,
                           levels = c("Daily or almost daily",
                                      "3-4 per week", "1-2 per week", "1-3 per month",
                                      "Special occasions only", "Never")),
    pa_met_excess_lvl = factor(pa_met_excess_lvl, ordered = T,
                               levels = c("Low", "Moderate", "High")),
    tg_log = log(tg),
    fast_cat = case_when(fasting_hrs > 0 & fasting_hrs <= 2 ~ "0to2hr",
                         fasting_hrs == 3 ~ "3hr",
                         fasting_hrs == 4 ~ "4hr",
                         fasting_hrs == 5 ~ "5hr",
                         fasting_hrs >= 6 ~ "6+hr",
                         TRUE ~ as.character(NA)),
    fasting_gt6hr = case_when(fasting_hrs >= 6 ~ 1,
                              fasting_hrs < 6 ~ 0,
                              TRUE ~ as.numeric(NA)),
    fasting_gt2hr = case_when(fasting_hrs >=2 ~ 1,
                              fasting_hrs < 2 ~ 0,
                              TRUE ~ as.numeric(NA)),
    fasting_gt3hr = case_when(fasting_hrs >=3 ~ 1,
                              fasting_hrs < 3 ~ 0,
                              TRUE ~ as.numeric(NA)),
    fasting_gt1hr = case_when(fasting_hrs >= 1 ~ 1,
                              fasting_hrs < 1 ~ 0,
                              TRUE ~ as.numeric(NA)),
    fasting_gt4hr = case_when(fasting_hrs >= 4 ~ 1,
                              fasting_hrs < 4 ~ 0,
                              TRUE ~ as.numeric(NA)),
    fasting_gt5hr = case_when(fasting_hrs >= 5 ~ 1,
                              fasting_hrs < 5 ~ 0,
                              TRUE ~ as.numeric(NA))
  ) 


# Load diet PCs -------------
source("../scripts/prep_data/build_dietPCs.R")
diet_id <- read.csv("../data/working/diet_pcs_prcomp_v3_noveg.rds")
diet_id <- diet_id %>% select(c(id, ffq_vars))


# Combine ukb data ---------------
ukb <- phenos_pre %>%
  left_join(geno_pre, by = "id") %>%
  left_join(t2d_id, by = "id") %>%
  left_join(diet_id, by = "id")


remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}

## REMOVE OUTLIERS FOR 5 SD OF CONTINUOUS VARIABLES****
ukb_rm5sd <- ukb %>% mutate_if(is.numeric, remove_outliers.fun) # no additional values removed


#############################################################################################
## Sample Exclusions
#############################################################################################

## Add variables for exclude/miss/keep to tabulate missingness
ukb_rm5sd <- ukb_rm5sd %>% mutate(
  exclude_taster = case_when(taster_status == "nontaster" | taster_status == "taster" |
                               taster_status == "supertaster" ~ 0,
                             taster_status != "nontaster" & taster_status != "taster" & taster_status != "supertaster" ~ 1,
                             is.na(taster_status) == T ~ 1),
  exclude_kcal = case_when(TCALS <600 | TCALS > 4800 ~ 1, 
                           TCALS >= 600 & TCALS <= 4800 ~ 0,
                           TRUE ~ as.numeric(NA)))

# Tabulate missing gPCs separately
miss_gPC <- ukb_rm5sd %>% select(paste0("gPC", 1:10), "id") %>% mutate(miss_gPC = ifelse(complete.cases(.), 0, 1)) %>%
  select(miss_gPC, "id")


## Tabulate miss/exclude/keep
ukb_miss <- ukb_rm5sd %>% mutate(
  miss_age = is.na(age),
  miss_sex = is.na(sex),
  miss_pa = is.na(pa_met_excess_lvl),
  miss_smoke.lab = is.na(smoke.lab),
  miss_ses = is.na(income_level.lab) | is.na(educ_level.lab),
  miss_dietPCs = is.na(dietPC1),
  miss_fib2cho = is.na(FIB2CHO),
  miss_t2d_case = is.na(t2d_case),
  miss_glu = is.na(glu),
  miss_alch = is.na(alch_freq.lab)
) %>% left_join(miss_gPC, by = "id") %>% 
  mutate(
    keep_primary = case_when(
      exclude_taster == 0 & miss_glu == F & exclude_kcal == F & miss_fib2cho == F ~ 1,
      exclude_taster == 1 | miss_glu == T | exclude_kcal == T & miss_fib2cho == T ~ 0,
      TRUE ~ as.numeric(NA))) %>% 
  mutate(
    keep_primary_controls = case_when(
      keep_primary == 1 & t2d_case == 0 ~ 1,
      keep_primary == 0 | t2d_case == 1 ~ 0,
      TRUE ~ as.numeric(NA))) %>%
  mutate(
    keep_base_controls = case_when(
      keep_primary_controls == 1 & miss_age == F & miss_sex == F & miss_gPC == F ~ 1,
      keep_primary_controls == 0 | miss_age == T | miss_sex == T | miss_gPC == T ~ 0,
      TRUE ~ as.numeric(NA))) %>% 
  mutate(
    keep_life_controls = case_when(
      keep_base_controls == 1 & miss_pa == F & miss_smoke.lab == F ~ 1,
      keep_base_controls == 0 & miss_pa == T | miss_smoke.lab == T ~ 0,
      TRUE ~ as.numeric(NA))) %>%
  mutate(keep_dietPC_controls = case_when(
    keep_life_controls == 1 & miss_dietPCs == F ~ 1,
    keep_life_controls == 0 | miss_dietPCs == T ~ 0,
    TRUE ~ as.numeric(NA))) %>% 
  select("id", starts_with(c("miss", "keep", "exclude")))

ukb_rm5sd <- ukb_rm5sd %>% left_join(ukb_miss %>% select("id", starts_with("keep")), by = "id")

exclude_table <- ukb_miss %>% 
  summarise(exclude_taster = table(exclude_taster),
          miss_age = table(miss_age),
          miss_sex = table(miss_sex),
          miss_gPC = table(miss_gPC),
          miss_pa = table(miss_pa),
          miss_smoke.lab = table(miss_smoke.lab),
          miss_alch = table(miss_alch),
          miss_t2d_case = table(miss_t2d_case),
          miss_glu = table(miss_glu)) %>% 
  t()

colnames(exclude_table) <- c("keep", "exclude")

keep_table <- ukb_miss %>% 
  summarise(
    N_for_primary = table(keep_primary),
    N_for_primary_controls = table(keep_primary_controls),
    N_for_base = table(keep_base_controls),
    N_for_lifestyle = table(keep_life_controls)) %>% 
  t() ; colnames(keep_table) <- c("exclude", "keep")


# write csvs with missing data tables
write.csv(exclude_table, "../output/tab_exclude.csv")
write.csv(keep_table, "../output/tab_keep.csv")


####################################
## Save rds & csv files
####################################

## Save data as .csv & RDS
write_csv(ukb_rm5sd, "../data/working/ukb_eur_unrel_processed_rm5sd.csv")
saveRDS(ukb_rm5sd, "../data/working/ukb_eur_unrel_processed_rm5sd.rds")


#EOF


