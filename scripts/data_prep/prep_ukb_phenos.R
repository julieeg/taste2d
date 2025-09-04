# Prepare UKB phenotypes 
# Last Updated: 09-24-2024
# Author: Julie E. Gervis



################################################################################
## Set Up & Load Required Libraries 
################################################################################

# load packages
library(tidyverse) ; library(data.table)
#installed_packages <- c("table1", "paletteer", "circlize")


# Load basic functions from pantry GitHub repo
#system("git clone https://github.com/julieeg/pantry.git")
invisible(lapply(list.files("../../pantry/functions/", full.names = T), function(f) {
  suppressMessages(suppressWarnings(source(f, echo = FALSE, verbose = FALSE)))
})) #source("../scripts/pantry/pantry.R")



################################################################################
## Load demographic & behavioral UKB files
################################################################################

print("Preparing basic phenotypes ...")

base_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
                     data.table=FALSE, stringsAsFactors=FALSE)

### Basic phenotypes --------------------

female_labs <- list("Female" = 0, "Male" = 1)
smoke_labs <- list("No answer"=-3, "Never"=0, "Previous"=1, "Current"=2)
smoking_num <- c("0" = 0, "1" = 1, "2" = 2, "-9" = -3)
med_mets_labs <- list("Cholesterol lowering" = 1, "Blood pressure" = 2, "Insulin" = 3,
  "None of the above" = -7, "Do not know" = -1, "Prefer not to answer" = -3)

# Assessment Center 
ac_labs <- c("Barts"=11012, "Birmingham" = 11021, "Bristol" =	11011, "Bury" =	11008, 
             "Cardiff" =	11003, "Cheadle (revisit)" =	11024, "Croydon" =	11020, 
             "Edinburgh" =	11005, "Glasgow" = 11004, "Hounslow" = 11018, "Leeds" = 11010,
             "Liverpool"=11016, "Manchester"=11001, "Middlesborough"=11017, "Newcastle" =11009, 
             "Nottingham"=11013, "Oxford"=11002, "Reading"=11007, "Sheffield"=11014, "Stockport (pilot)"=10003,
             "Stoke"=11006, "Swansea"=	11022,"Wrexham" =11023, "Cheadle (imaging)"=11025,
             "Reading (imaging)"=11026, "Newcastle (imaging)" =11027, "Bristol (imaging)"=11028)


base_phenos_id <- base_phenos %>% 
  select(id = f.eid,
         ac = f.54.0.0,
         ac_date = f.53.0.0,
         sex = f.31.0.0,
         age = f.21022.0.0,
         bmi = f.21001.0.0,
         sbp = f.4080.0.0,
         dbp = f.4079.0.0,
         waist = f.48.0.0,
         fasting_hrs = f.74.0.0,
         med_code = f.20003.0.0,
         med_mets = f.6177.0.0,
         smoking = f.20116.0.0) %>%
  mutate(
    female.lab = descr_label.fun(., "sex", female_labs),
    smoke.lab = descr_label.fun(., "smoking", smoke_labs),
    smoking.num = descr_label.fun(., "smoking", smoking_num),
    meds.lab = descr_label.fun(., "med_mets", med_mets_labs),
    ac.lab = descr_label.fun(., "ac", ac_labs))


# Exclude participants that withdrew consent
withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw/withdraw27892_232_14_Nov_2022.txt", 
                          what=character())


## Additional smoking variables --------------------

addn_smoke_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb45575.tab.gz", 
                       data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid,
         cigarettes_per_day = f.6183.0.0)


### Education level --------------------

## Coding based on: Ge T., et al. Cerebral Cortex 2019;29(8): 3471-3481.
educ_isced_level_labs <- list("Level 5" = 1, "Level 3" = 2, "Level 2" = 3, 
  "Level 2" = 4, "Level 5" = 5, "Level 4" = 6, "Level 1" = -7)  # NA = -3 or missing

educ_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672670.tab.gz",
                  data.table = FALSE, stringsAsFactors = FALSE) %>%
              select(id = f.eid, educ_level = f.6138.0.0) %>%
  mutate(educ_isced.lab = descr_label.fun(., "educ_level", educ_isced_level_labs))


### Alcohol intake frequency --------------------

alch_freq_labs <- list("Prefer not to answer" = -3, "Daily or almost daily" = 1, 
                       "3-4 per week" = 2, "1-2 per week" = 3, "1-3 per month" = 4, 
                       "Special occasions only" = 5, "Never" = 6) 
alch_num <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "-9" = -3)
drinker_status_labs <- c("Never" = 0, "Previous" = 1, "Current" = 2)


alch_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669173.tab.gz",
                 data.table=FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid,
         alch_freq = f.1558.0.0,
         alch_drinker_status = f.20117.0.0,
         alch_redwine_wk = f.1568.0.0,
         alch_wht_chm_wk = f.1578.0.0,
         alch_beer_wk = f.1588.0.0,
         alch_spirit_wk = f.1598.0.0,
         alch_fortwin_wk = f.1608.0.0,
         alch_othr_wk = f.5364.0.0) %>%
  
  mutate_at("alch_drinker_status", ~ifelse(is.na(.) | . == -3, NA, .)) %>%
  mutate(alch_drinker_status.lab = descr_label.fun(., "alch_drinker_status", drinker_status_labs)) %>%
  
  mutate(across(ends_with("_wk"), ~ifelse(is.na(.)==T | . %in% c(-1,-3), 0, .))) %>%
  mutate(alch_currdrinker = ifelse(is.na(alch_drinker_status) | alch_drinker_status == -3, NA, 
                                   ifelse(alch_drinker_status==2, 1, 0) )) %>%
  mutate(alch_neverdrinker = ifelse(!is.na(alch_drinker_status) & alch_drinker_status == 0, 1, 0)) %>%
  mutate(tot_alch_gm_wk_uk = alch_redwine_wk*16.8 + alch_wht_chm_wk*16.8 + alch_beer_wk*16 + 
           alch_spirit_wk*8 + alch_fortwin_wk*14.08 + alch_othr_wk*12) %>%
  mutate(tot_alch_drinks_per_week = tot_alch_gm_wk_uk / 14 ) %>%
  mutate(alch_freq.lab = descr_label.fun(., "alch_freq", alch_freq_labs),
         alch_freq.num = descr_label.fun(., "alch_freq", alch_num)) %>%
  select(id, alch_freq.num, alch_freq.lab, alch_drinker_status, tot_alch_drinks_per_week, tot_alch_gm_wk_uk)


### Physical Activity --------------------

pa_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb671173.tab.gz", 
      data.table=FALSE, stringsAsFactors=FALSE) %>% select(
        id = f.eid,
        iqpa_met = f.22040.0.0,
        
        # activity duration (mins) & frequency (days/week) 
        walking_dur = f.874.0.0, walking_frq = f.864.0.0, #duration of walks (min) & days/week of walks + 10 min
        moderate_dur = f.894.0.0, moderate_frq = f.884.0.0, #duration of moderate activity (min) & days/week of moderate activity
        vigorous_dur = f.914.0.0, vigorous_frq = f.904.0.0, #duration of vigorous activity (min) & days/week of vigorous activity
        pa_type = f.6164.0.0) %>%  
  
  # if walking_frq = -2 (Unable to walk) --> Recode to 0
  mutate_at("walking_frq", ~ifelse(. == -2, 0, .)) %>%
  
  # if activity duration <10 min/day --> Recode to 0
  mutate(across(ends_with("dur"), ~ifelse(.<10,0,.))) %>%
  
  # replace Do not know (-1), Prefer not to answer (-3) or missing (NA) with median
  mutate(across(ends_with("dur") | ends_with("frq") , ~ ifelse(.<0 | is.na(.)==T, 0, .))) %>%
  
  # multiply minutes per day per activity by excess MET score, per activity type
  mutate(met_excess_walk = ((walking_dur * 2.3)/60)*walking_frq,
         met_excess_mod = ((moderate_dur * 3.0)/60)*moderate_frq,
         met_excess_vig = ((vigorous_dur * 7.0)/60)*vigorous_frq) %>%
  
  # calculate excess met-hr/wk
  mutate(physact_met_excess = met_excess_walk + met_excess_mod + met_excess_vig) %>%
  mutate(physact_met_excess = winsorize(physact_met_excess))
  
# Add quantile variable for PA metabolic equivalents
physact_met_excess_lvls <- quantile(pa_id$physact_met_excess, probs=seq(0,1,0.33), include.lowest=F)[-1]

pa_id <- pa_id %>% mutate(
  physact_level = case_when(
    physact_met_excess < physact_met_excess_lvls[1] ~ "1",
    physact_met_excess >= physact_met_excess_lvls[1] & 
      physact_met_excess < physact_met_excess_lvls[2] ~ "2",
    physact_met_excess >= physact_met_excess_lvls[2] ~ "3")) %>%
  mutate(physact_level.lab = factor(physact_level, levels=c("1", "2", "3"), labels=c("Low", "Moderate", "High"))) %>%
  select(id, physact_met_excess, physact_level, physact_level.lab, iqpa_met)


### Income level --------------------

income_coding <- c(
  "1" = "Less than 18,000", "2" = "18,000 to 30,999", "3" = "31,000 to 51,999",
  "4" = "52,000 to 100,000", "5" = "Greater than 100,000", "-1" = "Do not know",
  "-3" = "Prefer not to answer"
  )

income_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                   data.table = FALSE, stringsAsFactors = FALSE) %>% 
  select(id=f.eid, income) %>%
  mutate(income = income_coding[as.character(f.738.0.0)]) %>%
  mutate(income_level = case_when(income == "Do not know" ~ mode(NA),
                                  income == "Prefer not to answer" ~ as.character(NA),
                                  income != "Do not know" & income != "Prefer not to answer" ~ as.character(income)))


#################
#### MERGE 1 ####
#################

phenos_id <- base_phenos_id %>% 
  full_join(addn_smoke_id, by = "id") %>%
  full_join(educ_id, by = "id") %>%
  full_join(alch_freq_id, by = "id") %>%
  full_join(pa_id, by = "id") %>%
  full_join(income_id, by = "id")
  
print(paste0("DONE: Basic phenotypes prepared for", nrow(phenos_id), " participants"))



################################################################################
## Load Biomarker & T2D variables
################################################################################

print("Preparing biomarker & T2D phenotypes ...")

## Biochemical parameters --------------------

biomark_fields <- c(
  chol = 30690, tg = 30870, ldl = 30780, hdl = 30760, glu = 30740, hba1c = 30750, crp = 30710
) ; biomark_vars <- setNames(paste0("f.", biomark_fields, ".0.0"), names(biomark_fields))

biomark_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz",
                 data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, 
         all_of(biomark_vars))


## T2D (Eastwood algorithm + HbA1c <5.7) --------------------

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
    meds_any_sr_ni_all == 0 & meds_any_sr_ni_ts_all == 0 ~ 0)) %>% 
  mutate(
    t2d_case = case_when(
      possible_t2dm_all == 1 ~ 1, 
      probable_t2dm_all == 1 ~ 1,
      hba1c.30750.NGSP.max >= 6.5 & possible_t1dm_all != 1 & probable_t1dm_all != 1 & dm_main != 1 & dm_secondary != 1 ~ 1,
      hba1c.30750.NGSP.max < 5.7 & dm_unlikely_all == 1 & dm_main != 1 & dm_secondary != 1 & t2d_med_any != 1 ~ 0,
      TRUE ~ as.numeric(NA) )) %>%
  mutate(
    t2d_case.f = case_when(
      as.numeric(t2d_case) == 0 ~ "Control",
      as.numeric(t2d_case) == 1 ~ "Case"),
    t2d_age_diagnosis = agedm_ts_or_ni_all ) %>%
  select("id", "t2d_case", "t2d_case.f", "hba1c_max"="hba1c.30750.NGSP.max", 
         "t2d_med_any", "t2d_age_diagnosis") 


print(paste0("DONE: Biochemical & anthropometric phenotypes prepared for ", nrow(t2d_id), " participants: ",
             "T2D Cases = ",  paste0(table(t2d_id$t2d_case.f))[1], ";",
            "T2D Controls = ", paste0(table(t2d_id$t2d_case.f))[2]) )


#################
#### MERGE 2 ####
#################

bmt2d_id <- left_join(biomark_id, t2d_id, by = "id")



################################################################################
## Load diet variables 
################################################################################

print("Preparing dietary data ...")

# =============
## FFQ data 
# =============

## Functions to prepare FFQ data - dapted with permission from KE Westerman 

# Function to convert frequency values to servings/day
ffq_freq_to_sev <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}

# function to recode negative values as meaninginful (missing)
neg_to_num <- function(x) {
  #x <- as.double(x) #"double" required to add values with decimals (previously, integer)
  case_when(
    x >= 0 ~ as.numeric(x), # -1 = "Do not know" ; -3 = "Prefer not to answer"
    x == -10 ~ 0.5, # -10 = "Less than 1 serving/day"
    TRUE ~ as.numeric(NA)
  )
}


### Food Frequency Questionnaire data --------------------
## Note: raw data included in the base_phenos file loaded in step 1

intake_fields <- c(cooked_veg = 1289, raw_veg = 1299,
                   fresh_fruit = 1309, dried_fruit = 1319,
                   bread_intake = 1438, bread_type = 1448, 
                   water = 1528, milk_type = 1418, spread_type = 1428,
                   spread_type_nonbutter=2654,
                   cereal_intake = 1458, cereal_type = 1468,
                   addsalt=1478, tea=1488, coffee = 1498, coffee_type = 1508,
                   hotdrink_temp = 1518) #**non_butter_spread_type

freq_fields <- c(oily_fish = 1329, nonoily_fish = 1339,
                 procmeat = 1349, poultry = 1359, cheese = 1408,
                 beef = 1369, lamb = 1379, pork = 1389)

ffq_fields<-c(intake_fields, freq_fields)
ffq_vars <- setNames(paste0("f.", ffq_fields, ".0.0"), names(ffq_fields))


ffq_id <- base_phenos %>% select(id=f.eid, ffq_vars) %>%
  mutate(whole_bread = case_when(
    bread_type == 3 ~ 1,
    bread_type %in% c(1, 2, 4) ~ 0,
    TRUE ~ as.numeric(NA) )) %>%
  
  # Recode frequency vars as intakes (servings/wk)
  mutate(across(names(freq_fields), ffq_freq_to_sev)) %>%
  
  # Recode negative values for intake vars as meaningful/missing
  mutate(across(names(intake_fields), neg_to_num)) %>%
  
  # Calculate 
  mutate(bread_intake = bread_intake / 7, # bread intake was provided in slices/week
         cereal_intake = cereal_intake / 7 # cereal intake was provided in bowls/week)
         ) %>%
  
  # Recode categorical food preference vars as binary (based on Cole et al.)
  mutate(
    bread_type_white_vs_brown_or_whole = case_when(      
      bread_type == 1 ~ 1, bread_type == 2 | bread_type == 3 | 
        bread_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_full_vs_low_or_nonfat = case_when(       
      milk_type == 1 ~ 1, milk_type == 2 | milk_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_rare_never_BIN = case_when(
      milk_type == 6 ~ 1, milk_type != 6 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_butter_vs_any_other = case_when(          
      spread_type == 1 ~ 1, spread_type == 2 | spread_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_rare_never_BIN = case_when(
      spread_type == 0 ~ 1, spread_type != 0 ~ 0, TRUE ~ as.numeric(NA)),
    cereal_type_sugar_vs_any_bran = case_when(
      cereal_type == 5 ~ 1, cereal_type != 5 ~ 0, TRUE ~ as.numeric(NA)),
    coffee_type_decaf_vs_regular = case_when(              
      coffee_type == 1 ~ 1, coffee_type == 2 | coffee_type == 3 | 
        coffee_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    addsalt_freq_QT = addsalt,        # INCLUDE
    addsalt_always_often_vs_nrs = case_when(
      addsalt == 3 | addsalt == 4 ~ 1, addsalt == 1 | addsalt == 2 ~ 0, TRUE ~ as.numeric(NA)),
    hotdrink_temp_hot_or_vhot_vs_warm = case_when(
      hotdrink_temp == 1  ~ 1, hotdrink_temp == 3 | hotdrink_temp == 2 ~ 0, TRUE ~ as.numeric(NA))
    )


## Create categorical diet variables with descriptive levels

ffq_id <- ffq_id %>% 
  mutate(
    bread_type.lab = case_when(
      bread_type == 1 ~ "White", bread_type == 2 ~ "Brown", bread_type == 3 ~ "Wholemeal/Wholegrain",
      bread_type == 4 ~ "Other", bread_type == -1 ~ "Do not know", bread_type == -3 ~ "Prefer not to answer",
      TRUE ~ as.character(NA)),
    
    milk_type.lab = case_when(
      milk_type == 1 ~ "Full cream", milk_type == 2 ~ "Semi-skimmed", milk_type == 3 ~ "Skimmed",
      milk_type == 4 ~ "Soy", milk_type == 5 ~ "Other", milk_type == 6 ~ "Never/rarely have milk",
      milk_type == -1 ~ "Do not know", milk_type == -3 ~ "Prefer not to answer", 
      TRUE ~ as.character(NA)),
    
    spread_type.lab = case_when(
      spread_type == 1 ~ "Butter/spreadable butter", spread_type == 2 ~ "Flora Pro-Active/Benecol",
      spread_type == 3 ~ "Other spread/margarine", spread_type == 0 ~ "Never/rarely use spread",
      spread_type == -1 ~ "Do not know", spread_type == -3 ~ "Prefer not to answer"),
    
    spread_type_nonbutter.lab = case_when(
      spread_type_nonbutter == 4 ~ "Soft (tub) margarine", spread_type_nonbutter == 5 ~	"Hard (block) margarine",
      spread_type_nonbutter == 6 ~ "Olive oil based spread (eg: Bertolli)",
      spread_type_nonbutter == 7 ~ "Polyunsaturated/sunflower oil based spread (eg: Flora)",
      spread_type_nonbutter == 2 ~ "Flora Pro-Active or Benecol",
      spread_type_nonbutter == 8 ~ "Other low or reduced fat spread",
      spread_type_nonbutter == 9 ~ "Other type of spread/margarine", spread_type_nonbutter == -1 ~	"Do not know",
      spread_type_nonbutter == -3 ~ "Prefer not to answer"),
    
    cereal_type.lab = case_when(
      cereal_type == 1 ~ "Bran cereal (e.g. All Bran, Branflakes)", cereal_type == 2 ~ "Biscuit cereal (e.g. Weetabix)",
      cereal_type == 3 ~ "Oat cereal (e.g. Ready Brek, porridge)", cereal_type == 4 ~ "Muesli",
      cereal_type == 5 ~ "Other (e.g. Cornflakes, Frosties)", cereal_type == -1 ~ "Do not know",
      cereal_type == -3 ~ "Prefer not to answer"),
    
    coffee_type.lab = case_when(
      coffee_type == 1 ~ "Decaffeinated coffee (any type)", coffee_type == 2 ~ "Instant coffee",
      coffee_type == 3 ~ "Ground coffee (include espresso, filter etc)", coffee_type == 4 ~ "Other type of coffee",
      coffee_type == -1 ~ "Do not know", coffee_type == -3 ~ "Prefer not to answer"),
    
    addsalt.lab = case_when(
      addsalt == 1 ~ "Never/Rarely", addsalt == 2 ~ "Sometimes", addsalt == 3 ~ "Often",
      addsalt == 4 ~ "Always", TRUE ~ as.character(NA)),
    
    hotdrink_temp.lab = case_when(
      hotdrink_temp == 1 ~ "Very hot", hotdrink_temp == 2 ~ "Hot", hotdrink_temp == 3 ~ "Warm",
      TRUE ~ as.character(NA))
    )

print(paste0("Dietary data from UKB FFQs prepared for N = ", nrow(ffq_id), " participants"))



# =======================================
## 24HR data (for total energy intake)
# =======================================


### Functions to prepare 24HR data across multiple measurements (4/1 year) 
## Adapted with permission from KEW 

fetch_diet_fields <- function(fieldIDs, df, coding=FALSE) {
  # Given a list of fields constituting a food group:
  # - Determine the set of 24HR that are valid for that food group
  # - Recode the relevant variables based on their codings if necessary
  # - Sum over all fields for that food group within each instance
  # - Take the mean food group quantity over all instances
  diet_field_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.26002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    instance_fields <- paste0("f.", fieldIDs, ".", i, ".0")  # Variable names for all fields in instance "i"
    instance_df <- df[, instance_fields, drop=FALSE]
    if (coding) {  # Recode the variable if necessary (for food groups)
      instance_df <- mutate_all(instance_df, ~codings[as.character(.)])
    }
    ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
           rowSums(instance_df, na.rm=TRUE), NA) } ) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  diet_mean <- rowMeans(diet_field_df, na.rm=TRUE)
  ifelse(is.nan(diet_mean), NA, diet_mean)
}

check_num_valid_24hr <- function(df) {
  valid_24hr_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.26002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    valid_24hr
  }) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  rowSums(valid_24hr_df, na.rm=TRUE)
}


### Nutrient variables -------------------------------------------------
## Merge "Typical diet yesterday" 

diet_vars <- c("TCALS_raw", "TCALS", "CHO", "CHO_kcal", "FAT", "FAT_kcal", 
               "PRO", "PRO_kcal", "ALC")

diet_id <- left_join(fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb670995.tab.gz", 
                           data.table=FALSE, stringsAsFactors=FALSE),
                     fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz",
                           data.table=FALSE, stringsAsFactors = FALSE) %>% 
                       select(f.eid, starts_with("f.100020")),
                     by = "f.eid") %>%
  
  mutate(TCALS_raw = fetch_diet_fields("26002", .),
         CHO = fetch_diet_fields("26013", .),
         CHO_kcal = fetch_diet_fields("26013", .),
         FAT = fetch_diet_fields("26008", .),
         FAT_kcal = fetch_diet_fields("26008", .),
         SFA = fetch_diet_fields("26014", .),
         SFA_kcal = fetch_diet_fields("26014", .),
         PRO = fetch_diet_fields("26005", .),
         PRO_kcal = fetch_diet_fields("26005", .),
         ALC = fetch_diet_fields("26030", .),
         num_recalls = check_num_valid_24hr(.) ) %>%
  
  # Energy from kJ to kcal, likely includes alcohol
  mutate(TCALS_raw = TCALS_raw / 4.18) %>% 
  
  mutate_at(vars(CHO_kcal, PRO_kcal), ~. * 4, ) %>%  # Nutrients from g to kcals (other than alcohol)
  mutate_at(vars(FAT_kcal), ~. * 9) %>%
  mutate(TCALS = (CHO_kcal + PRO_kcal + FAT_kcal)) %>%
  
  #mutate_at(vars(all_of(diet_vars)), remove_outliers.fun) %>%
  select(id=f.eid, all_of(diet_vars))


#################
#### MERGE 3 ####
#################

diet_all_id <- left_join(ffq_id, diet_id, by = "id")



################################################################################
## Retrieve genetic ancestry (from Pan UKB)
################################################################################

### Compile Pan-UKBB genetic PCs (for European subset) --------------------

anc_rel_id <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                    data.table=FALSE, stringsAsFactors=FALSE) %>% 
  mutate(f.eid = as.character(f.eid),
         unrelated = !related_return2442) %>%
  select(id=f.eid, ancestry=pop_return2442, unrelated,
         one_of(paste0("PC", 1:20, "_return2442"))) %>%
  rename_at(vars(contains("PC")), ~gsub("PC", "gPC", gsub("_return2442", "", .))) %>%
  mutate(id=as.integer(id))

print("Breakdown of available data by relatedness & ancestry from PanUKBB: ")
print(table(anc_rel_id$unrelated, anc_rel_id$ancestry))



############################################
####  Merge all raw datasets --> phenos ####
############################################

phenos <- phenos_id %>%
  left_join(bmt2d_id, by="id") %>%
  left_join(diet_all_id, by="id") %>%
  left_join(anc_rel_id, by="id") %>%
  filter(!(id %in% withdrawn_consent)) %>% # Remove participants who withdrew consent
  mutate(id = format(id, scientific=FALSE)) 


## ALL participants
phenos %>% write_csv("../data/processed/ukb_phenos.csv")

# Filter to unrelated subjects of EUR ancestry 
phenos %>%
  filter(unrelated == TRUE & ancestry == "EUR") %>%
  write_csv("../data/processed/ukb_phenos_unrelated_EUR.csv")

print("Done preparing UKB Phenotype data. Datasets are ready for analysis.")

## Print Data Dictionary
sink("../data/processed/ukb_phenos_str.txt")
str(phenos)
sink()


##EOF


