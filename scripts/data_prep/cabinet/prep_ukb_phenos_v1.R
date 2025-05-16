# prep_ukb_phenos_v1.R

# load packages
library(tidyverse)
library(data.table)


# load basic functions for data preparation & cleaning
source("../scripts/pantry/pantry.R")




##############################################
##  Load demographic & lifestyle variables  ##
##############################################

print("Preparing basic phenotypes ...")

base_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE)


### Basic phenotypes -----------------------------------------------

smoke_labs <- list("No answer"=-3, "Never"=0, "Previous"=1, "Current"=2)
smoking_num <- c("0" = 0, "1" = 1, "2" = 2, "-9" = -3)
med_mets_labs <- list("Cholesterol lowering" = 1, "Blood pressure" = 2, "Insulin" = 3,
                      "None of the above" = -7, "Do not know" = -1, "Prefer not to answer" = -3)

base_phenos_id <- base_phenos %>% 
  select(id = f.eid,
         ac = f.54.1.0,
         ac_date = f.53.1.0,
         bmi = f.21001.1.0,
         sbp = f.4080.1.0,
         dbp = f.4079.1.0,
         waist = f.48.1.0,
         fasting_hrs = f.74.1.0,
         smoking = f.20116.1.0) %>%
  mutate(
    smoke.lab = descr_label.fun(., "smoking", smoke_labs),
    smoking.num = descr_label.fun(., "smoking", smoking_num),
    meds.lab = descr_label.fun(., "med_mets", med_mets_labs)
  )

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw/withdraw27892_232_14_Nov_2022.txt", what=character())


### Education level ---------------------------------------------------

## Coding based on: Ge T., et al. Cerebral Cortex 2019;29(8): 3471-3481.
educ_level_labs <- list(
  "None of the above" = -7, 
  "Prefer not to answer"= -3, 
  "College or university degree" = 1, 
  "A/AS levels or equivalent" = 2, 
  "O/GCSE levels or equivalent" = 3, 
  "CSEs or equivalent" = 4,
  "NVQ/HND or equivalent" = 5, 
  "Other professional qualifications" = 6)

educ_isced_level_labs <- list("Level 5" = 1, "Level 3" = 2, "Level 2" = 3, 
                              "Level 2" = 4, "Level 5" = 5, "Level 4" = 6, "Level 1" = -7)  # NA = -3 or missing

educ_years_labs <- list("20"=1, "13"=2, "10"=3, "10"=4, "19"=5, "15"=6, "7"=-7)


educ_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672670.tab.gz",
                 data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, educ_level = f.6138.1.0) %>%
  mutate(educ_level.lab = descr_label.fun(., "educ_level", educ_level_labs),
         educ_isced.lab = descr_label.fun(., "educ_level", educ_isced_level_labs),
         educ_years = as.numeric(descr_label.fun(., "educ_level", educ_years_labs))
  )



### Alcohol intake frequency -------------------------------------------

alch_freq_labs <- list("Prefer not to answer" = -3, "Daily or almost daily" = 1, 
                       "3-4 per week" = 2, "1-2 per week" = 3, "1-3 per month" = 4, "Special occasions only" = 5, 
                       "Never" = 6) 

alch_num <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "-9" = -3)


alch_freq_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669173.tab.gz",
                      data.table=FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid,
         alch_freq = f.1558.1.0) %>%
  mutate(alch_freq.lab = descr_label.fun(., "alch_freq", alch_freq_labs),
         alch_freq.num = descr_label.fun(., "alch_freq", alch_num))



### Physical Activity -------------------------------------------

pa_fields <- c("walking_dur", "walking_frq", "moderate_dur", "moderate_frq",
               "vigorous_dur", "vigorous_frq")

pa_vars <- c("pa_met_excess", "pa_met_excess_lvl")


pa_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb671173.tab.gz", 
               data.table=FALSE, stringsAsFactors=FALSE) %>% select(
                 id = f.eid,
                 #iqpa_met = f.22040.1.0,
                 
                 # activity duration (mins) & frequency (days/week) 
                 walking_dur = f.874.1.0, walking_frq = f.864.1.0, #duration of walks (min) & days/week of walks + 10 min
                 moderate_dur = f.894.1.0, moderate_frq = f.884.1.0, #duration of moderate activity (min) & days/week of moderate activity
                 vigorous_dur = f.914.1.0, vigorous_frq = f.904.1.0, #duration of vigorous activity (min) & days/week of vigorous activity
                 pa_type = f.6164.1.0
               ) %>%  
  
  # if walking_frq = -2 (Unable to walk) --> Recode to 0
  mutate_at("walking_frq", ~ifelse(walking_frq == -2, 0, walking_frq)) %>%
  
  # replace Do not know (-3), Prefer not to answer (-3) or missing (NA) with median
  mutate(across(pa_fields, ~ median_imp.fun(.))) %>%
  
  # convert duration from min to hr
  mutate(across(c(walking_dur, moderate_dur, vigorous_dur), ~ . /60)) %>% 
  
  # calculate excess met, per activity type
  mutate(met_excess_walk = (walking_dur * walking_frq * 2.3),
         met_excess_mod = (moderate_dur * moderate_frq * 3.0),
         met_excess_vig = (vigorous_dur * vigorous_frq * 7.0)) %>%
  
  # calculate excess met-hr/wk
  mutate(pa_met_excess = met_excess_walk + met_excess_mod + met_excess_vig) %>% 
  
  # create variable for physical activity LEVEL
  mutate(
    pa_met_excess_lvl = ifelse(
      pa_met_excess < 10, "Low", ifelse(
        pa_met_excess > 10 & pa_met_excess < 49.8, "Moderate", ifelse(
          pa_met_excess >= 50, "High", NA)))
  ) %>% 
  select(id, all_of(pa_vars))



### Income level ---------------------------------------------------

income_fields <- c(income = 738)
income_coding <- c(
  "1" = "Less than 18,000",
  "2" = "18,000 to 30,999",
  "3" = "31,000 to 51,999",
  "4" = "52,000 to 100,000",
  "5" = "Greater than 100,000",
  "-1" = "Do not know",
  "-3" = "Prefer not to answer"
)

income_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                   data.table = FALSE, stringsAsFactors = FALSE) 
income_id <- income_id %>%
  mutate(income = income_coding[as.character(f.738.1.0)]) %>%
  select(id=f.eid, income) %>%
  mutate(income_level = case_when(income == "Do not know" ~ mode(NA),
                                  income == "Prefer not to answer" ~ as.character(NA),
                                  income != "Do not know" & income != "Prefer not to answer" ~ as.character(income)))


################################
##  Load Biomarker variables  ##
################################

print("Preparing biomarker & T2D phenotypes ...")

## Biochemical parameters ---------------------------------

biomark_fields <- c(
  chol = 30690, tg = 30870, ldl = 30780, hdl = 30760, hba1c = 30750, crp = 30710
) ; biomark_vars <- setNames(paste0("f.", biomark_fields, ".0.0"), names(biomark_fields))

biomark_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz",
                    data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, 
         all_of(biomark_vars))


print(paste0("DONE: Basic phenotypes prepared for", nrow(base_phenos_id), " participants"))


#################
#### MERGE  ####
#################

phenos_id <- base_phenos_id %>% 
  full_join(educ_id, by = "id") %>%
  full_join(alch_freq_id, by = "id") %>%
  full_join(pa_id, by = "id") %>%
  full_join(income_id, by = "id") %>%
  full_join(biomark_id, by = "id")


## ALL participants
phenos_id %>% rename_with(., ~paste0(., ".1"), -"id") %>% 
  write_csv("../data/processed/ukb_v1_phenos.csv")


## Subet to EUR/unrealated  & save ----------------
ukb_eur_unrelated <- readRDS("../data/processed/ukb_analysis_EUR.rda") %>%
  filter(unrelated == T)

phenos_id <- phenos_id %>% 
  filter(id %in% ukb_eur_unrelated$id)

phenos_id %>% write_csv("../data/processed/ukb_v1_phenos_EUR.csv")



########################################
## Postprocessing 
########################################


## basic phenotypes  ---------------
pheno_vars <- c("ac", "ac_date", "smoke.lab", "alch_freq.lab", "pa_met_excess_lvl", 
                "pa_met_excess", "income", "educ_level.lab", "educ_isced.lab", "educ_years",
                "glu", "hba1c_max", "bmi", "sbp", "dbp", "crp", "chol", "tg", "ldl" , "hdl")

phenos_id_postprocess <- phenos_id %>%   
  mutate(
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
      TRUE ~ as.character(NA))
  ) %>%
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
    smoke_level.lab = factor(smoke_level.lab, levels = c("Never", "Previous", "Current")),
    alch_freq.lab = factor(alch_freq.lab, ordered = T,
                           levels = c("Daily or almost daily",
                                      "3-4 per week", "1-2 per week", "1-3 per month",
                                      "Special occasions only", "Never")),
    physact_level = factor(pa_met_excess_lvl, ordered = T,
                           levels = c("Low", "Moderate", "High")),
    tg_log = log(tg),
    fast_cat = case_when(fasting_hrs >= 0 & fasting_hrs <= 2 ~ "0to2hr",
                         fasting_hrs == 3 ~ "3hr",
                         fasting_hrs == 4 ~ "4hr",
                         fasting_hrs == 5 ~ "5hr",
                         fasting_hrs >= 6 & fasting_hrs <= 24 ~ "6+hr",
                         TRUE ~ as.character(NA))
  )








##EOF
