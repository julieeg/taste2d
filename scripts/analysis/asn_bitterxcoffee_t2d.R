
## load packages & functions
lapply(c("tidyverse", "data.table", "paletteer", "ggpubr", "emmeans", "ComplexHeatmap",
         "knitr"),  library, character.only = TRUE)

# load summary table function
source("../scripts/pantry/print_summary_table_fun.R", echo=F)
source("../scripts/pantry.R")


# load raw data
base <- readRDS("../data/processed/ukb_postprocessed_EUR_20241227.rda") 


#############################
# Build analytical dataset ##
#############################

# build btc
btc <- base %>% 
  filter(n_geno==1) %>%
  select(id, age, sex, starts_with("gPC"), ancestry, diplos_common_AA,
         taste_diplos, taste_diplos.num, rs713598_G, rs1726866_G, rs10246939_C,
         bmi, physact_level.lab, smoke_level.lab, alch_freq.lab, fasting_hrs, fast_cat,
         income_level.lab, educ_level.lab, cigarettes_per_day,
         glu, glu2hr, t2d_case, t2d_case.f, t2d_med_any,
         coffee_QT, tea_QT, energy=TCALS, carb_pct=CHO_pct, fat_pct=FAT_pct, prot_pct=PRO_pct, fiber=FIBER,
         starts_with("dietPC")) 


# Add categorical coffee intake variable
btc <- btc %>% mutate(
  coffee_cat5 = case_when(
    coffee_QT == 0 ~ "0 cups/d",
    coffee_QT > 0 & coffee_QT <= 1 ~ "1 cups/d",
    coffee_QT == 2 ~ "2 cups/d",
    coffee_QT == 3 ~ "3 cups/d",
    coffee_QT == 4 ~ "4 cups/d",
    coffee_QT >= 5 ~ "5+ cups/d" 
  )) %>%
  mutate(coffee_cat5 = factor(coffee_cat5, levels=c(paste(c(0:4, "5+"), "cups/d")))) %>%
  mutate(coffee_cat4 = case_when(
    coffee_QT == 0 ~ "0 cups/d",
    coffee_QT > 0 & coffee_QT <= 1 ~ "<=1 cups/d",
    coffee_QT %in% c(2,3) ~ "2-3 cups/d",
    coffee_QT %in% c(4,5) ~ "4-5 cups/d",
    coffee_QT >= 6 ~ "6+ cups/d"
  )) %>%
  mutate(coffee_cat4 = factor(coffee_cat4, levels=c(paste(c("0", "<=1", "2-3", "4-5", "6+"), "cups/d")))) %>%
  mutate(diplos_commonAA = factor(diplos_common_AA, levels=c("AVI/AAV", "AVI/AVI", "AVI/PAV", "PAV/PAV", "AAV/PAV"))) %>%
  mutate(diplos_commonAA = relevel(diplos_commonAA, ref="AVI/AVI"))

taste_diplotypes = c("AVI/AVI", "AVI/PAV", "PAV/PAV")
coffee_cat4 = c(paste(c("0", "<=1", "2-3", "4-5", "6+"), "cups/d"))
coffee_cat5 = c(paste(c(0:4, "5+"), "cups/d"))


################################
## primary statistical models ##
################################

# build models
m.base=paste(c("age", "sex", paste0("gPC",1:10), "bmi"), collapse="+")
m.lifestyle=paste0(c(m.base, "smoke_level.lab", "physact_level.lab", "alch_freq.lab", "tea_QT"), collapse="+")
m.ses=paste0(c(m.lifestyle, "income_level.lab", "educ_level.lab"), collapse="+")
models.l = list(Base=m.base, Lifestyle=m.lifestyle, SES=m.ses)


# associations between TAS2R38 with T2D status
tab_bitter_t2d.l <- lapply(1:3, function(m) {
  print_glm(exposure = "taste_diplos", outcome="t2d_case", print_trend=T, covariates = models.l[[m]],
            label=paste0("T2D~TAS2R38+", names(models.l)[m]), data=btc) 
}) ; names(tab_bitter_t2d.l) = names(models.l)

tab_bitterADD_t2d.l <- lapply(1:3, function(m) {
  print_glm(exposure = "taste_diplos.num", outcome="t2d_case", print_trend=T, covariates = models.l[[m]],
            label=paste0("T2D~TAS2R38+", names(models.l)[m]), data=btc) 
}) ; names(tab_bitter_t2d.l) = names(models.l)


# associations of coffee intake (continuous) with T2D status
tab_coffeeQT_t2d.l <- lapply(1:3, function(m) {
  print_glm(exposure = "coffee_QT", outcome="t2d_case", print_trend=T, covariates = models.l[[m]],
            label=paste0("T2D~Coffee_QT+", names(models.l)[m]), data=btc)
}) ; names(tab_coffeeQT_t2d.l) = names(models.l)

# associations of coffee intake (5-level) with T2D status
tab_coffeeCAT5_t2d.l <- lapply(1:3, function(m) {
  print_glm(exposure = "coffee_cat5", outcome="t2d_case", print_trend=T, covariates = models.l[[m]],
            label=paste0("T2D~Coffee_5level+", names(models.l)[m]), data=btc)
}) ; names(tab_coffee5_t2d.l) = names(models.l)

tab_maineffects_t2d.l = list(
  bitter_t2d=tab_bitter_t2d.l, bitterADD_t2d=tab_bitternum_t2d.l,
  coffeeQT_t2d=tab_coffeeQT_t2d.l, coffeeCAT5_t2d=tab_coffee5_t2d.l
)


##################
## interactions ##
##################

get_interaction <- function(outcome="t2d_case", interaction, covariates, data) {
  m <- glm(formula(paste0(outcome, "~", interaction, "+", covariates)), data=btc, family=binomial("logit"))
  aov <- anova(m) ; aov.sel <- strsplit(interaction, "[*]")[[1]]
  return(list(glm=m, sum=summary(m), anova=as.data.frame(aov[c(aov.sel, gsub("[*]", ":", interaction)),])))
}

tab_bitterxcoffeeQT_t2d.l <-  lapply(1:3, function(m) {
  get_interaction(interaction = "taste_diplos*coffee_QT", covariates = models.l[[m]], data=btc) 
})

tab_bitterxcoffeeCAT5_t2d.l <-  lapply(1:3, function(m) {
  get_interaction(interaction = "taste_diplos*coffee_cat5", covariates = models.l[[m]], data=btc) 
}) 

tab_bitterADDxcoffeeQT_t2d.l <-  lapply(1:3, function(m) {
  get_interaction(interaction = "taste_diplos.num*coffee_QT", covariates = models.l[[m]], data=btc[1:1000,]) 
}) 

tab_bitterADDxcoffeeCAT5_t2d.l <-  lapply(1:3, function(m) {
  get_interaction(interaction = "taste_diplos.num*coffee_cat5", covariates = models.l[[m]], data=btc) 
}) 


tab_interactions_t2d.l = list(
 bitterxcoffee=tab_bitterxcoffeeQT_t2d.l, bitterxcoffeeCAT5=tab_bitterxcoffeeCAT5_t2d.l, 
 bitterADDxcoffeeQT=tab_bitterADDxcoffeeQT_t2d.l, bitterADDxcoffeeCAT5=tab_bitterADDxcoffeeCAT5_t2d.l 
)

## save preliminary data frames & .rdas
list(tab_maineffects_t2d.l, tab_interactions_t2d.l) %>% 
  saveRDS("../data/processed/archive/analysis/ASN_bitterxcoffee_t2d_main&interaction_v2.rda")

saveRDS(btc, file="../data/processed/ASN_btc_data.rda")
  

####################################
## Run in taste-stratified models ##
####################################

# are the protective effects of coffee, stronger among supertasters than nontasters?

# associations of coffee intake (5-level) with T2D status
tab_coffee5_strat_bitter_t2d.l <- lapply(1:2, function(m) {
  temp <- lapply(as.list(taste_diplotypes), function(t) {
    print_glm(exposure = "coffee_cat5", outcome="t2d_case", print_trend=T, covariates = models.l[[m]],
              label=paste0("T2D~Coffee_5level+", names(models.l)[m]), data=btc %>% filter(taste_diplos == t))
    }) ; names(temp) <- taste_diplotypes ; 
    return(temp)
  }) ; names(tab_coffee5_strat_bitter_t2d.l) = names(models.l)



############################################################
## Test associations with glycemic traits in T2D controls ##
############################################################

# build dataframe with T2D controls, only
btc_controls <- btc %>% filter(t2d_case == 0) %>%
  mutate(glu2hr_mgdl = glu2hr*18.018,
         glu_mgdl = glu*18.018)

# associations between TAS2R38 x coffee with random glucose
m_bitterxcoffeeCAT5_glu2hr <- print_lm_interaction(
  interaction="coffee_cat5", outcome="glu2hr_mgdl", model=paste0(models.l$Lifestyle, "+fasting_hrs"),
  data=btc_controls, label_model = "Lifestyle", label_interaction = "TAS2R38 x Coffee")

m_bitterxcoffeeCAT5_glu <- print_lm_interaction(
  interaction="coffee_cat5", outcome="glu_mgdl", model=paste0(models.l$Lifestyle, "+fasting_hrs"),
  data=btc_controls, label_model = "Lifestyle", label_interaction = "TAS2R38 x Coffee")


saveRDS(list(m_bitterxcoffeeCAT5_glu2hr, m_bitterxcoffeeCAT5_glu), 
        "../data/processed/archive/analysis/ASN_bitterxcoffee_glu_interaction.rda")



#############################################
## Test associations with longitudinal T2D ##
#############################################

library(survival)
library(survminer)
library(ggsurvfit)

# grab incident t2d data
t2d_incid <- fread("../data/t2d_incident_phenotype_07022024_SH.csv") %>%
  select(id=eid, baseline_date, t2d, t2d_date, followup_date, censor_date, death_date) %>%
  mutate(across(ends_with("date"), ~as.Date(.))) %>%
  mutate(
    time_to_t2d = t2d_date - baseline_date, # 
    time_to_censor = censor_date - baseline_date) %>%
  mutate(time_to_event = ifelse(!is.na(time_to_t2d), time_to_t2d, time_to_censor))

# create combined dataset with bitter/coffee phenotypes & incident phenotypes
incid = full_join(btc, t2d_incid, by = "id") %>%
  filter(t2d_case == 0)


## Run CoxPH models for TAS2R38 and coffee main effects on T2D indicence
cox_bitter <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos+", models.l$Lifestyle)), data = incid, exp)
cox_bitterADD <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos.num+", models.l$Lifestyle)), data = incid)
cox_coffeeCAT5 <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ coffee_cat5+", models.l$Lifestyle)), data = incid)
cox_coffeeQT <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ coffee_QT+", models.l$Lifestyle)), data = incid)


cox_bitterxcoffeeQT <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos*coffee_QT+", models.l$Lifestyle)), data = incid)
cox_bitterxcoffeeCAT5 <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos*coffee_cat5+", models.l$Lifestyle)), data = incid)
cox_bitterADDxcoffeeQT <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos.num*coffee_QT+", models.l$Lifestyle)), data = incid)
cox_bitterADDxcoffeeCAT5 <- coxph(formula(paste0("Surv(time_to_event, t2d) ~ taste_diplos.num*coffee_cat5+", models.l$Lifestyle)), data = incid)

results_cox.l = list(cox_bitter, cox_bitterADD, cox_coffeeCAT5, cox_coffeeQT,
                     cox_bitterxcoffeeQT, cox_bitterxcoffeeCAT5, cox_bitterADDxcoffeeQT, cox_bitterADDxcoffeeCAT5) 

saveRDS(results_cox.l, "../data/processed/archive/analysis/ASN_bitterxcoffee_t2d_incident_maininteraction.rda")


## survival models
surv_bitterxcoffeeCAT5 <- survfit(formula(paste0("Surv(time_to_event, t2d) ~ coffee_cat5+", models.l$Lifestyle)), data=incid)
print(surv_bitterxcoffeeCAT5) ; summary(surv_bitterxcoffeeCAT5)










