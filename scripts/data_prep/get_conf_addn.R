
library(tidyverse)
library(data.table)


# ========================
## Winsorize data by SD
# ========================
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


# =========================
## Add descriptive labels 
# =========================

#e.g., labs_vals = female_labs <- list("Female" = 0, "Male" = 1)
descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}


descr_label_ordered.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; temp <- factor(temp, levels=names(labs_vals)) 
  return(temp)
}




### Alcohol intake frequency -------------------------------------------

alch_freq_labs <- list("Prefer not to answer" = -3, "Daily or almost daily" = 1, 
                       "3-4 per week" = 2, "1-2 per week" = 3, "1-3 per month" = 4, "Special occasions only" = 5, "Never" = 6) 
alch_num <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "-9" = -3)
drinker_status_labs <- c("Never" = 0, "Previous" = 1, "Current" = 2)


## updated: 09-24-2024
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
  
  ## Get alcohol drinker status (never/former/current)
  mutate_at("alch_drinker_status", ~ifelse(is.na(.) | . == -3, NA, .)) %>%
  mutate(alch_drinker_status.lab = descr_label.fun(., "alch_drinker_status", drinker_status_labs)) %>%
  mutate(alch_drinker_status.lab = factor(descr_label.fun(., "alch_drinker_status", c("Never" = 0, "Previous" = 1, "Current" = 2)), levels=c("Never", "Previous", "Current"))) %>%
  
  ## Recode #drinks per week: do not know/no answer asmissing & winsorise
  mutate(across(ends_with("_wk"), ~ifelse(is.na(.)==T | . %in% c(-1,-3), 0, .))) %>%
  mutate(across(ends_with("_wk"), ~winsorize(.))) %>%
  
  ## Make variable for curent/nondrinkers
  mutate(alch_currdrinker = ifelse(is.na(alch_drinker_status) | alch_drinker_status == -3, NA, 
                                   ifelse(alch_drinker_status==2, 1, 0) )) %>%
  mutate(alch_neverdrinker = ifelse(!is.na(alch_drinker_status) & alch_drinker_status == 0, 1, 0)) %>%
  mutate(alch_gm_per_wk = alch_redwine_wk*16.8 + alch_wht_chm_wk*16.8 + alch_beer_wk*16 + 
           alch_spirit_wk*8 + alch_fortwin_wk*14.08 + alch_othr_wk*12) %>%
  mutate(alch_drinks_per_week = alch_gm_per_wk / 14 ) %>%
  
  # prepare alcohol frequency variable as ordered factor & numeric
  mutate(alch_freq.lab = descr_label.fun(., "alch_freq", alch_freq_labs),
         alch_freq.num = descr_label.fun(., "alch_freq", alch_num)) %>%
  mutate(alch_freq.lab = factor(alch_freq.lab, levels= c(names(alch_freq_labs)[-1]) )) %>%
  select(id, alch_freq.num, alch_freq.lab, alch_drinker_status, alch_drinker_status.lab, alch_drinks_per_week, alch_gm_per_wk)




### Physical Activity -------------------------------------------

pa_fields <- c("walking_dur", "walking_frq", "moderate_dur", "moderate_frq",
               "vigorous_dur", "vigorous_frq")

pa_vars <- c("physact_met_excess", "physact_level")

pa_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb671173.tab.gz", 
               data.table=FALSE, stringsAsFactors=FALSE) %>% select(
                 id = f.eid,
                 iqpa_met = f.22040.0.0,
                 
                 # activity duration (mins) & frequency (days/week) 
                 walking_dur = f.874.0.0, walking_frq = f.864.0.0, #duration of walks (min) & days/week of walks + 10 min
                 moderate_dur = f.894.0.0, moderate_frq = f.884.0.0, #duration of moderate activity (min) & days/week of moderate activity
                 vigorous_dur = f.914.0.0, vigorous_frq = f.904.0.0, #duration of vigorous activity (min) & days/week of vigorous activity
                 pa_type = f.6164.0.0
               ) %>%  
  
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

physact_met_excess_lvls <- quantile(pa_id$physact_met_excess, probs=seq(0,1,0.33), include.lowest=F)[-1]

pa_id <- pa_id %>% mutate(
  physact_level = case_when(
    physact_met_excess < physact_met_excess_lvls[1] ~ "1",
    physact_met_excess >= physact_met_excess_lvls[1] & 
      physact_met_excess < physact_met_excess_lvls[2] ~ "2",
    physact_met_excess >= physact_met_excess_lvls[2] ~ "3")) %>%
  mutate(physact_level.lab = factor(physact_level, levels=c("1", "2", "3"), labels=c("Low", "Moderate", "High"))) %>%
  
  select(id, physact_met_excess, physact_level, physact_level.lab, iqpa_met)



### Additional & medication variables --------------

meds_female_labs = c("Cholesterol lowering" = 1, "Blood pressure" = 2, "Insulin" = 3, "Hormone replacement therapy" = 4, 
                     "Oral contraceptive pill" = 5, "None of the above" = -7)
meds_male_labs = c("Cholesterol lowering" = 1, "Blood pressure" = 2, "Insulin" = 3, 
                   "None of the above" = -7)

addn_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb45575.tab.gz", 
                       data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, 
         cigarettes_per_day = f.6183.0.0,
         meds_male = f.6177.0.0)

addn_med_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669176.tab.gz", 
                  data.table=FALSE, stringsAsFactors=FALSE) %>% 
  select(id = f.eid, meds_female = f.6153.0.0)

addn_all_id <- addn_id %>% full_join(addn_med_id, by = "id") %>%
  mutate(across(c(meds_male, meds_female), ~ifelse(. %in% c(-3, -1) | is.na(.) == T, NA, .))) %>%
  mutate(meds_male.lab = descr_label.fun(., "meds_male", meds_male_labs),
         meds_female.lab = descr_label.fun(., "meds_female", meds_female_labs)) %>%
  mutate_at("cigarettes_per_day", ~ifelse(.==-10, 0.5, ifelse(. == -1, NA, .))) %>%
  select(id, cigarettes_per_day, meds_male, meds_male.lab, meds_female, meds_female.lab)



### Combine & save

left_join(alch_id, pa_id, by = "id") %>% left_join(addn_all_id, by = "id") %>%
  saveRDS("../data/processed/ukb_conf_addn_09242024.rda")



