
## load packages & functions
lapply(c("tidyverse", "data.table", "paletteer", "ggpubr", "emmeans", "ComplexHeatmap",
         "knitr"),  library, character.only = TRUE)

# load summary table function
source("../scripts/pantry/print_summary_table_fun.R", echo=F)
source("../scripts/pantry.R")


#############################
# Build analytical dataset ##
#############################

## RUN ON CLUSTER TO CREATE DATAFRAMER ========================================= 
# load dataframes
df <- fread("/humgen/florezlab/users/msevilla/coffe_t2d/df_filtered.csv") %>%
  select(id=f.eid, 
         age=Age_at_recruitment, sex=Sex,  smoking_qualy, alcohol_four, smoking_quantity, 
         tea_intake_four,  PA_miss, MET_minutes , paste0("PC", 1:10, "_return2442"), tea=tea_intake_four,
         Body_mass_index_BMI.21001, coffee_intake_cat4, coffee_cat) 

taste <- readRDS("../data/processed/ukb_postprocessed_EUR_20241227.rda") %>%
  select(id, taste_diplos, taste_diplos_ADD = taste_diplos.num, rs713598_G, 
         rs1726866_G, rs10246939_C, starts_with("dietPC"), ancestry, t2d_case,
         smoke_level.lab, physact_level.lab, alch_freq.lab, bmi, starts_with("gPC"),
         coffee_QT)

## QC
df_taste <- df %>% left_join(taste, by = "id") %>%
  mutate(smoking_qualy = factor(smoking_qualy, levels=c("Current_everyday", "Current_occasionally", "Past", "Never or unknown"))) %>%
  mutate(smoking_quantity = as.numeric(smoking_quantity)) %>%
  mutate(alcohol_four = factor(alcohol_four, levels = c("Daily or almost daily", "Three or four times a week", 
                                                        "Once or twice a week", "One to three times a month", 
                                                        "Special occasions only", "Never"))) %>%
  mutate(tea = as.numeric(tea)) %>%
  mutate(coffee_intake_cat4 = factor(coffee_intake_cat4, levels=paste(c("0", "<=1", "2", "3", "4", ">=5")))) %>%
  filter(taste_diplos %in% c("AVI/AVI", "AVI/PAV", "PAV/PAV")) %>%
  mutate(taste_diplos = factor(taste_diplos, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV"))) %>%
  filter(!is.na(coffee_intake_cat4))

fwrite(df_taste, "../data/processed/asn_coffee_t2d_incident.csv")

## =============================================================================


#############################################
## Test associations with longitudinal T2D ##
#############################################

library(survival)
library(survminer)
library(ggsurvfit)

taste_diplotypes = c("AVI/AVI", "AVI/PAV", "PAV/PAV")
coffee_intake_cat4 = c(paste(c("0", "<=1", "2", "3", "4", "5+")))

exposures.l <- list("taste_diplos", "taste_diplos_ADD", "coffee_intake_cat4", "coffee_QT")
interactions.l <- list("taste_diplos*coffee_intake_cat4", "taste_diplos_ADD*coffee_intake_cat4", "taste_diplos*coffee_QT", "taste_diplos_ADD*coffee_QT")

df <- fread("../data/processed/asn_coffee_t2d_incident.csv") %>%
  mutate(coffee_intake_cat4 = as.factor(coffee_intake_cat4)) %>% 
  filter(coffee_intake_cat4 %in% coffee_intake_cat4) %>%
  mutate(coffee_intake_cat4 = relevel(coffee_intake_cat4, ref="0"))

# grab incident t2d data
t2d_incid <- fread("../data/t2d_incident_phenotype_07022024_SH.csv") %>%
  select(id=eid, baseline_date, t2d, t2d_date, followup_date, censor_date, death_date) %>%
  mutate(across(ends_with("date"), ~as.Date(.))) %>%
  mutate(
    time_to_t2d = t2d_date - baseline_date, # 
    time_to_censor = censor_date - baseline_date) %>%
  mutate(time_to_event = ifelse(!is.na(time_to_t2d), time_to_t2d, time_to_censor))

# create combined dataset with bitter/coffee phenotypes & incident phenotypes
incid = left_join(df, t2d_incid, by = "id") %>%
  filter(t2d_case == 0) %>%
  mutate(coffee_intake_cat4 = factor(coffee_intake_cat4, levels = c("0", "<=1", "2", "3", "4", ">=5")))


## Run CoxPH models for TAS2R38 and coffee main effects on T2D indicence

m_lifestyle <- paste0(c("age", "sex", paste0("PC", 1:10, "_return2442"), 
                    "smoking_qualy", "smoking_quantity", "alcohol_four", "tea",  
                    "PA_miss", "MET_minutes" , "tea"), collapse = "+")
m_lifestyle_bmi <- paste0(c(m_lifestyle, "bmi"), collapse = "+")

## Lifestyle model
coxph_maineffects.l <- lapply(exposures.l, function(exposure) {
  coxph(formula(paste0("Surv(time_to_event, t2d) ~", exposure, "+", m_lifestyle)), data = incid)
}) ; coxph_interacteffects.l <- lapply(interactions.l, function(interaction) {
  coxph(formula(paste0("Surv(time_to_event, t2d) ~", interaction, "+", m_lifestyle)), data = incid)
}) ; coxph_results.l <- list(coxph_maineffects.l, coxph_interacteffects.l)

## Lifestyle+BMI model
coxph_maineffects_adjbmi.l <- lapply(exposures.l, function(exposure) {
  coxph(formula(paste0("Surv(time_to_event, t2d) ~", exposure, "+bmi+", m_lifestyle_bmi)), data = incid)
}) ; coxph_interacteffects_adjbmi.l <- lapply(interactions.l, function(exposure) {
  coxph(formula(paste0("Surv(time_to_event, t2d) ~", exposure, "+bmi+", m_lifestyle_bmi)), data = incid)
}) ; coxph_results_adjbmi.l <- list(coxph_maineffects_adjbmi.l, coxph_interacteffects_adjbmi.l)


saveRDS(coxph_results.l, "../data/processed/archive/analysis/ASN_bitterxcoffee_t2d_incident_maininteraction_v2.rda")
saveRDS(coxph_results_adjbmi.l, "../data/processed/archive/analysis/ASN_bitterxcoffee_t2d_incident_maininteraction_adjbmi_v2.rda")
saveRDS(incid, "../data/processed/asn_coffee_incidt2d_data_postprocessed.rda")


## Additive interaction ========================

## Lifestyle+BMI model
add_interactions.l <- list("taste_diplos:coffee_intake_cat4", "taste_diplos_ADD:coffee_intake_cat4",
                           "taste_diplos:coffee_QT", "taste_diplos_ADD:coffee_QT")

## Additive Interactions
incid <- incid %>% 
  mutate(taste_diplo_AVIPAV=ifelse(taste_diplos == "AVI/PAV", 1, ifelse(taste_diplos=="AVI/AVI",0, NA)),
         taste_diplo_PAVPAV=ifelse(taste_diplos == "PAV/PAV", 1, ifelse(taste_diplos=="AVI/AVI",0, NA))) %>% 
  mutate(coffee_intake_cat4_1=ifelse(coffee_intake_cat4=="<=1",1, ifelse(coffee_intake_cat4=="0", 0, NA)),
         coffee_intake_cat4_2=ifelse(coffee_intake_cat4=="2",1, ifelse(coffee_intake_cat4=="0", 0, NA)),
         coffee_intake_cat4_3=ifelse(coffee_intake_cat4=="3",1, ifelse(coffee_intake_cat4=="0", 0, NA)),
         coffee_intake_cat4_4=ifelse(coffee_intake_cat4=="4",1, ifelse(coffee_intake_cat4=="0", 0, NA)),
         coffee_intake_cat4_5=ifelse(coffee_intake_cat4==">=5",1, ifelse(coffee_intake_cat4=="0", 0, NA))) %>%
  mutate(taste_diplo_AVIPAV_coffee_cat1 = ifelse(taste_diplo_AVIPAV == 1 & coffee_intake_cat4_1 == 1, 1, 0),
         taste_diplo_PAVPAV_coffee_cat1 = ifelse(taste_diplo_PAVPAV == 1 & coffee_intake_cat4_1 == 1, 1, 0),
         taste_diplo_AVIPAV_coffee_cat2 = ifelse(taste_diplo_AVIPAV == 1 & coffee_intake_cat4_2 == 1, 1, 0),
         taste_diplo_PAVPAV_coffee_cat2 = ifelse(taste_diplo_PAVPAV == 1 & coffee_intake_cat4_2 == 1, 1, 0),
         taste_diplo_AVIPAV_coffee_cat3 = ifelse(taste_diplo_AVIPAV == 1 & coffee_intake_cat4_3 == 1, 1, 0),
         taste_diplo_PAVPAV_coffee_cat3 = ifelse(taste_diplo_PAVPAV == 1 & coffee_intake_cat4_3 == 1, 1, 0),
         taste_diplo_AVIPAV_coffee_cat4 = ifelse(taste_diplo_AVIPAV == 1 & coffee_intake_cat4_4 == 1, 1, 0),
         taste_diplo_PAVPAV_coffee_cat4 = ifelse(taste_diplo_PAVPAV == 1 & coffee_intake_cat4_4 == 1, 1, 0),
         taste_diplo_AVIPAV_coffee_cat5 = ifelse(taste_diplo_AVIPAV == 1 & coffee_intake_cat4_5 == 1, 1, 0),
         taste_diplo_PAVPAV_coffee_cat5 = ifelse(taste_diplo_PAVPAV == 1 & coffee_intake_cat4_5 == 1, 1, 0)
         )

interaction_terms <- c("taste_diplo_AVIPAV_coffee_cat1", "taste_diplo_PAVPAV_coffee_cat1")

coxph(formula(paste0("Surv(time_to_event, t2d) ~", interaction_terms[2], "+bmi+", m_lifestyle_bmi)), data = incid)

coxph(formula(paste0("Surv(time_to_event, t2d) ~", add_interactions.l[[1]], "+bmi+", m_lifestyle_bmi)), data = incid) %>%
  saveRDS("../data/processed/asn_coffe_add_interaction.rda")


