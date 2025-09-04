# Primary analysis for Taste2D
# Last Updated: 05-16-2025
# Author: Julie E. Gervis


################################################################################
## Set Up & Load Required Libraries 
################################################################################

setwd("~/Documents/GitHub/taste2d/run/")

## Load packages & basic functions --------------------

library(tidyverse) ; library(data.table) ; library(emmeans)

## Install pantry functions
#system("git clone https://github.com/julieeg/pantry.git")
#installed_packages <- c("table1", "paletteer", "circlize")
invisible(lapply(list.files("../../pantry/functions/", full.names = T), function(f) {
  suppressMessages(suppressWarnings(source(f, echo = FALSE, verbose = FALSE)))
}))


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_ms.rda")) %>%
  mutate(
    glu_mgdl = glu*18.018,
    glu2hr_mgdl = glu2hr*18.018) %>%
  # Add LDL, HDL and TG in mg/dL
  mutate(
    ldl_fasting_gt6_mgdl = ifelse(fasting_hrs >= 6, ldl*38.67, NA),
    hdl_fasting_gt6_mgdl = ifelse(fasting_hrs >= 6, hdl*38.67, NA),
    tg_fasting_gt6_mgdl = ifelse(fasting_hrs >= 6, tg*88.57, NA)
  )

dim(analysis)  # N=241378


# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common_AA))>0.005))
taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")

glycemic_vars_mgdl.l <- list("Glucose (mg/dL)"="glu_mgdl", 
                             "0-2hr Glucose (mg/dL)"="glu2hr_mgdl")

# create output directories
outDir=paste0("../data/processed/manuscript")



################################################################################
## Descriptive analysis of participant characteristics (Table 1)
################################################################################

# ============================================
## Distributions of haplotypes & diplotypes 
# ============================================

## Haplotypes
analysis %>%
  select(id, haplo_0_aa, haplo_1_aa) %>%
  pivot_longer(-id, values_to="haplo") %>%
  group_by(haplo) %>% reframe(n=n(), pct=(n()/nrow(.))*100) %>%
  arrange(n) %>% mutate(haplo_order = factor(haplo, levels=haplo[order(n, decreasing = T)])) %>%
  mutate(diplo_legend = factor(ifelse(haplo_order == "AVI", "Nontaster haplotype", 
                                      ifelse(haplo_order == "PAV", "Taster haplotype", "Non-canonical haplotype")), 
                               levels=paste(c("Nontaster", "Taster", "Non-canonical"), "haplotype"))) %>%
  #Filter to "common" haplotypes (>1% of the cohort; >0.005, frequency)
  filter(pct>0.05) %>%
  saveRDS("../data/processed/manuscript/tab_descr_haplos_common.rda")


## Diplotypes --------------
diplos_gt005 <- names(which(prop.table(table(analysis$diplos_common_AA))>0.005))
analysis <- analysis %>% mutate(
  diplos_AA = factor(ifelse(diplos_common_AA %in% diplos_gt005, diplos_common_AA, NA), levels=c(diplos_gt005)))


## Table
analysis %>% 
  filter(complete.cases(diplos_AA)) %>% #diplos_common_AA %in% diplos_common_AA) %>% 
  reframe(n=table(diplos_AA), pct=prop.table(table(diplos_AA)), diplos_AA=names(table(diplos_AA))) %>%
  arrange(-n) %>% 
  mutate(color = factor(ifelse(diplos_AA == "AVI/AVI", "Nontaster",
                               ifelse(diplos_AA == "AVI/PAV", "Taster",
                                      ifelse(diplos_AA == "PAV/PAV", "Supertaster", "Non-canonical"))),
                        levels=c("Nontaster","Taster","Supertaster","Non-canonical"))) %>%
  mutate(diplos_AA_order=factor(diplos_AA, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV", "AVI/AAV", "AAV/PAV")),
         diplo_legend = factor(paste(color, "diplotype"), levels=
                                 paste(c("Nontaster","Taster","Supertaster","Non-canonical"), "diplotype"))) %>%
  filter(pct>0.005) %>%
  saveRDS("../data/processed/manuscript/tab_descr_diplos_common.rda")


## TAS2R38 diplotypes by assessment center --------------------

print_summary_table(data=analysis, vars_to_summarize = c("taste_diplos" = "TAS2R38 Diplotypes"),
                    var_strata = "ac.f", p_types = "descriptive", factor_vars = "taste_diplos") %>%
  write.csv(paste0("../data/processed/manuscript/tab_descr_diplosByac.csv"))


## Demographic, behavioral & cardiometabolic risk factors 
vars_to_summarize <- c(
  age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking Status", cigarettes_per_day="Cigarettes per day (among smokers)", 
  physact_level.lab="Physical Activity Level", physact_met_excess="Physical Activity, excess MET/wk",
  alch_freq.lab="Alcohol Frequency", alch_1to4d.lab="Alcohol Frequency, 1-4 per week",
  alch_drinker_status.lab = "Alcohol Drinker", #alch_heavydrinker = "Heavy alcohol drinker",
  alch_drinks_per_week = "Drinks per week",
  kcal_plaus = "Energy intake, kcal/d",
  educ_isced.lab="Education Level", income_level.lab = "Income Level",
  sbp_adj="SBP, mmHg (med adjusted)", dbp="DBP, mmHg (med adjusted)",
  tg_fasting_gt6_mgdl="Triglyceride, mg/dL", ldl_fasting_gt6_mgdl="LDL, mg/dL", hdl_fasting_gt6_mgdl="HDL, mg/dL", hba1c = "HbA1c",
  raw_veg_QT="Raw vegetables", coffee_QT="Coffee, cups/day", tea_QT="Tea, cups/day", 
  addsalt_3lvl.lab="Added salt"
  )


# ================================================
## Summarise Table 1 over diplotypes with >0.005 
# ================================================

# Set of diplotypes for table 1 with frequency >0.005
diplo_sets <- c(diplos_common_AA="diplo005", taste_diplos="diplo_taste")
analysis <- analysis %>% mutate(
  diplos_common_AA = factor(
    diplos_common_AA, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV", "AVI/AAV", "AAV/PAV")))


## Confirm correct factor coding for categorical variables
lapply(c("smoke_level.lab", "alch_freq.lab", "physact_level.lab", "income_level.lab", "educ_isced.lab"), function(v) {
  str(analysis[[v]]) }) #--> Recode education as factor variable
analysis <- analysis %>% mutate_at("educ_isced.lab", ~factor(., levels=paste0("Level ", 1:5)))


## Run Table 1 descriptives over two sets of diplotypes 

for (i in 1:length(diplo_sets)) {
  
  # choose diplotype set
  diplo_strata_var = names(diplo_sets)[i]
  diplo_strata_order =  levels(analysis[[diplo_strata_var]])
  
  # Print summary table with digits: 1=means, 0=pct, 4=p-values
  print_summary_table(
    vars_to_summarize = vars_to_summarize, 
    var_strata = diplo_strata_var, var_strata_order = c(diplo_strata_order),
    factor_vars = c("smoke_level.lab", "physact_level.lab", "alch_freq.lab",
                    "income_level.lab", "educ_isced.lab", "addsalt_3lvl.lab"),
    p_print = T, p_adjust = "none", digits = c(1,1,3), data=analysis) %>%
    
    fwrite(paste0(outDir, "/tab_descr_basic_",diplo_sets[i], ".csv"), row.names = T)
}

## Additional pairwise comparisons for continuous outcomes 
summary(lm(alch_drinks_per_week~taste_diplos, data=analysis))$coef

# Add variable for 1 or two PAV haplotypes
analysis <- analysis %>% mutate(taste_pav = factor(ifelse(taste_diplos %in% c("AVI/PAV", "PAV/PAV"), "PAV carrier", "AVI/AVI"),
                                                   levels=c("AVI/AVI", "PAV carrier")))
summary(lm(cigarettes_per_day~taste_pav, data=analysis %>% filter(smoke.lab == "Current")))

#log-transformed TG
hist(log(analysis$tg_fasting_gt6_mgdl))
summary(lm(log(tg_fasting_gt6_mgdl)~taste_diplos.num, data=analysis))


# ==================================================
## Descriptive analysis of glucose homeostasis 
# ==================================================

# Primary glucose outcomes
glucose_vars_mgdl <- c(glu_mgdl = "Random glucose, mg/dL",
                       glu2hr_mgdl = "0-2hr Glucose, mg/dL") 

analysis <- analysis %>% mutate(glu_mgdl = glu*18.018,
                                glu2hr_mgdl = glu2hr*18.018)


## Summarise glucose outcomes by diplotype with age, sex, ac adjustment
print_summary_table(data=analysis, vars_to_summarize = glucose_vars_mgdl, 
                    var_strata = "taste_diplos", var_strata_order = c("AVI/AVI", "AVI/PAV", "PAV/PAV"),
                    p_adjust=c("age", "sex", "ac.f"), digits = c(3,0,6)) %>%
  mutate(P_adjust="Sex+Age+AC", .before="Total") %>%
  write.csv(paste0(outDir, "/tab_descr_glu_diplo_sexageac.csv"))



################################################################################
## Primary Analysis of glucose homeostasis by TAS2R38 diplotype
## Base and BMI-adjusted models
################################################################################

## Outcomes
glu_vars.l <- list(glu_mgdl = "Glucose, mg/dL", glu2hr_mgdl = "0-2hr Glucose, mg/dL")
names(glu_vars.l) <- c("glu_mgdl", "glu2hr_mgdl")

## Model covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+ac.f+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")

models.l = list("Base"= m1_base, "BMI"=m2_bmi)
models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)


# ================================================================
## Associations of TAS2R38 with Random & 0-2hr Glucose (mg/dL)
# ================================================================

## Random Glucose --------------------

# Categorical diplotype 
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome="Glucose (mg/dL)",
           label=names(models.l)[m], covariates = models.l[[m]], lm_trend = T, 
           data = analysis) } )) %>%
  write.csv(paste0(outDir, "/tab_res_rg_diplocat_basebmi.csv"), row.names = F)


## Continuous diplotypes --------------------
do.call(rbind.data.frame, lapply(1:length(models.l), function(m){
  print_lm(outcome = "glu_mgdl", exposure = "taste_diplos.num", label.outcome="Glucose (mg/dL)",
           covariates = models.l[[m]], label=names(models.l)[m]) })) %>% 
  write.csv(paste0(outDir, "/tab_res_rg_diplocont_basebmi.csv"), row.names = F)


## 0-2 hr Glucose --------------------

# Categorical diplotypes
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu2hr_mgdl", label.outcome="0-2hr Glucose (mg/dL)",
           label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], 
           lm_trend = T, data = analysis %>% filter(fast_cat == "0to2hr")) })) %>% 
  write.csv(paste0(outDir, "/tab_res_2hg_diplocat_basebmi.csv"), row.names = F)


## Continuous diplotypes
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m){
  print_lm(outcome = "glu2hr_mgdl", exposure = "taste_diplos.num", label.outcome="0-2hr Glucose (mg/dL)", 
           covariates = models.nofast.l[[m]], label=names(models.nofast.l)[m]) })) %>% 
  write.csv(paste0(outDir, "/tab_res_2hg_diplocont_basebmi.csv"), row.names = F)


# ===============================================================================
## Associations of TAS2R38 with Glucose (mg/dL) over subsequent fasting windows
# ===============================================================================

fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

# Categorical diplotypes
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
    print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = paste0(fast, "-Glucose (mg/dL)"), 
             label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
             data = analysis %>% filter(fast_cat == fast)) })) })) %>% 
  write.csv(paste0(outDir, "/tab_res_gluXfast_diplocat_basebmi.csv"), row.names = F)


# Continuous diplotypes
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
    print_lm(exposure = "taste_diplos.num", outcome = "glu_mgdl", label.outcome = paste0(fast, "-Glucose (mg/dL)"), 
             label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
             data = analysis %>% filter(fast_cat == fast)) })) })) %>% 
  write.csv(paste0(outDir, "/tab_res_gluXfast_diplocont_basebmi.csv"), row.names = F)

## Interaction in the diet model with Taste x Fasting category
int_pvals <- sapply(models.l, function(m) {
  anv<-anova(lm(formula(paste0("glu_mgdl~taste_diplos.num*fast_cat+", m)), data=analysis))
  anv[grep(":", rownames(anv)),5] # Interaction P-value 
  }) 

# Base        BMI 
# 0.02130727 0.02251386 


# ==============================================================================
## Calculate EMM random glucose & glucose levels at each fasting interval 
## in the BMI-adjusted model
# ==============================================================================

## Random glucose --------------------

get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$BMI, 
            label.outcome = "Glucose (mg/dL)", label="BMI", reference = "AVI/AVI", 
            data=analysis )$emm %>% 
  write.csv(paste0(outDir, "/tab_res_emm_glu_diplo_bmi.csv"), row.names = F)


## Glucose by fasting time (separating 6-12 and 12-24) --------------

fast_cat_with_6to12 <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6to12hr", "12+hr"))
analysis <- analysis %>% 
  mutate(fast_cat_with_6to12 = factor(ifelse(fasting_hrs >= 6 & fasting_hrs <= 12, "6to12hr", 
                                             ifelse(fasting_hrs >12, "12+hr", fast_cat)),
                                      levels=c("0to2hr", "3hr", "4hr", "5hr", "6to12hr", "12+hr")))

do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(f) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$BMI, 
              label.outcome = paste0(f, "-Glucose (mg/dL)"), label="BMI", reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat_with_6to12==f))$emm })) %>% 
  mutate(fast=gsub("[-].*","",outcome),.before=exposure) %>%
  write.csv(paste0(outDir, "/tab_res_emm_gluXfast_6to12_diplo_bmi.csv"), row.names = F)


## Calculate combined association of TAS2R38 diplotype with glucose >= 3 hr -----------

# Categorical 
make_pretty_lm(print_lm(exposure="taste_diplos", outcome="glu_mgdl", covariates = models.l$BMI,
         label="Glucose, >= 3 hr", data = analysis %>% filter(fast_cat != "0to2hr")))

# Continuous
make_pretty_lm(print_lm(exposure="taste_diplos.num", outcome="glu_mgdl", covariates = models.l$BMI,
         label="Glucose, >= 3 hr", data = analysis %>% filter(fast_cat != "0to2hr")))

## Estimated marginal mean
emm_gt3hr <- get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$BMI, 
            label.outcome = "3+hr-Glucose (mg/dL)", label="BMI", reference = "AVI/AVI", 
              data=analysis %>% filter(fasting_hrs >= 3 & fasting_hrs <= 12))
emm_gt3hr$emm %>% mutate(fast = "3+hr", .after=model) %>%
  fwrite("../archive/data/processed/manuscript/tab_sens_emm_glugt3hr_bmi.csv")



###############################
#### Post-hoc comparisons  ####
###############################

pct_beta_change.fun <- function(b_base, b_adj) {
  b_pct <- ((b_adj - b_base)/b_base)*100
  return(b_pct)
}

## Quantify the % confounding for BMI, lifestyle & food choices on models for glu
m_base <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[1]])), data=analysis))
m_bmi <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[2]])), data=analysis))

# BMI
pct_beta_change.fun(m_base$coef[2,1], m_bmi$coef[2,1]) #(1.063%)

## Quantify the % confounding for BMI, lifestyle & food choices on models for glu
m_base <- summary(lm(formula(paste0("glu2hr~taste_diplos.num+", models.nofast.l[[1]])), data=analysis))
m_bmi <- summary(lm(formula(paste0("glu2hr~taste_diplos.num+", models.nofast.l[[2]])), data=analysis))

# BMI
pct_beta_change.fun(m_base$coef[2,1], m_bmi$coef[2,1]) 



################################################################################
## Independent Replication in published GWAS (of 2-hr glucose from MAGIC)
################################################################################

tas2r38_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C")
tas2r38 <- c("taste_diplos.num", tas2r38_snps)

# ===========================================================================
## Summarize SNP-level associations with random glucose & 0-2hr glucose
# ===========================================================================

# Variant-level associations with 0-2 hr glucose 
do.call(rbind.data.frame, lapply(tas2r38_snps, function(snp) {
  print_lm(exposure = snp, outcome = "glu2hr_mgdl", label.outcome="0-2hr Glucose (mg/dL)",
           label="BMI", covariates = models.nofast.l$BMI, 
           lm_trend = T, data = analysis %>% filter(fast_cat == "0to2hr")) 
  })) %>%
  write.csv(paste0(outDir, "/tab_res_2hg_tas2r38snps_bmi.csv"), row.names = F)


# Variant-level associations with 0-2 hr glucose 
do.call(rbind.data.frame, lapply(tas2r38_snps, function(snp) {
  print_lm(exposure = snp, outcome = "glu2hr_mgdl", label.outcome="0-2hr Glucose (mg/dL)",
           label="BMI", covariates = models.nofast.l$BMI, 
           lm_trend = T, data = analysis %>% filter(fast_cat == "0to2hr")) 
})) %>%
  write.csv(paste0(outDir, "/tab_res_2hg_tas2r38snps_bmi.csv"), row.names = F)


# ===========================================================
## Calculate EMM 0-2 hr glucose for each TAS2R38 variant
# ===========================================================

tas2r38_snps_cat <- paste0(tas2r38_snps, ".cat")

analysis <- analysis %>% mutate(
  rs713598_G.cat = case_when(
    rs713598_G == 0 ~ "0_0", rs713598_G == 1 ~ "0_1", rs713598_G == 2 ~ "1_1"
  ),
  rs1726866_G.cat = case_when(
    rs1726866_G == 0 ~ "0_0",
    rs1726866_G == 1 ~ "0_1",
    rs1726866_G == 2 ~ "1_1"
  ),
  rs10246939_C.cat = case_when(
    rs10246939_C == 0 ~ "0_0",
    rs10246939_C == 1 ~ "0_1",
    rs10246939_C == 2 ~ "1_1")
)

# Calculate EMM 
do.call(rbind.data.frame, lapply(tas2r38_snps_cat, function(snp) {
  get_emm.fun(exposure = snp, outcome = "glu2hr_mgdl", covars = models.nofast.l$BMI, 
              label="BMI", reference = "0_0", data=analysis %>% filter(fast_cat == "0to2hr"))$emm 
  })) %>%
  write.csv("../data/processed/manuscript/tab_sens_emm_gluXfast_tas2r38snps_bmi.csv"), row.names = F)
getwd()

# ==============================================
## Summarize published GWAS summary statistics
# ==============================================

# Make data.frame with MAGIC meta-analysis results
cbind.data.frame(
  Study=rep("Saxena_MAGIC_2hG_AdjBMI", 3), 
  SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
  Beta=c(-0.0790, -0.0510, -0.0570), SE=c(0.022, 0.021, 0.021), P_value=c(0.0002571, 0.01228, 0.006042)) %>%
  mutate(Beta_mgdl=Beta*18.018, SE_mgdl=SE*18.018) %>%
  fwrite(paste0(outDir, "/tab_repl_glu_tas2r38_magic.csv"))



################################################################################
## Testing for behaviorally-mediated associations (adj for diet+lifestyle)
################################################################################

## Add binary lifestyle variables --------------------

analysis <- analysis %>%
  mutate(smoker_current_vs_never_BIN = ifelse(
    smoke_level.lab=="Current", 1, ifelse(smoke_level.lab=="Never",0,NA)),
    smoker_former_vs_never_BIN = ifelse(
      smoke_level.lab=="Previous", 1, ifelse(smoke_level.lab=="Never",0,NA)),
    physact_low_vs_high_BIN = ifelse(
      physact_level.lab=="Low", 1, ifelse(physact_level.lab=="High",0,NA)),
    physact_mod_vs_high_BIN = ifelse(
      physact_level.lab=="Moderate", 1, ifelse(physact_level.lab=="High",0,NA)),
    alch_daily_rarelynever_BIN = ifelse(
      alch_freq.lab=="Daily or almost daily", 1, 
      ifelse(alch_freq.lab %in% c("Never", "Special occasions only"),0, NA)),
    alch_weekly_rarelynever_BIN = ifelse(
      alch_freq.lab %in% c("3-4 per week", "1-2 per week"), 1,
      ifelse(alch_freq.lab %in% c("Never", "Special occasions only"),0, NA)),
    alch_monthly_rarelynever_BIN = ifelse(
      alch_freq.lab=="1-3 per month", 1, ifelse(
        alch_freq.lab %in% c("Never", "Special occasions only"),0, NA))
  )

## Create vector of additional diet/lifestyle traits 
behav_covars.labs <- c(
  raw_veg_QT="Raw vegetable intake",
  cooked_veg_QT="Cooked vegetable intake",
  fresh_fruit_QT="Fresh fruit intake",
  dried_fruit_QT="Dried fruit intake",
  oily_fish_QT="Oily fish intake",
  nonoily_fish_QT="Non-oily fish intake",
  procmeat_QT="Processed meat intake",
  poultry_QT="Poultry intake",
  cheese_QT="Cheese intake",
  beef_QT="Beef intake",
  lamb_QT="Lamb intake",
  pork_QT="Pork intake",
  bread_type_white_vs_brown_or_whole_BIN="Prefer white>brown/whole bread",
  bread_intake_QT="Bread intake",
  milk_type_full_vs_low_or_nonfat_BIN="Prefer full>low/no-fat milk",
  cereal_type_sugar_vs_any_bran_BIN="Prefer sugary>bran cereal",
  coffee_type_decaf_vs_regular_BIN="Prefer decaf>caffeinated coffee",
  cereal_intake_QT="Cereal intake",
  spread_type_butter_vs_any_other_BIN="Prefer butter>other spreads",
  coffee_QT="Coffee intake",
  tea_QT="Tea intake",
  water_QT="Water intake",
  addsalt_always_often_vs_nrs_BIN="Always/often add salt",
  hotdrink_temp_hot_or_vhot_vs_warm_BIN="Prefer hot>warm drink temps",
  smoker_current_vs_never_BIN = "Current (vs never) smoker",
  smoker_former_vs_never_BIN = "Former (vs never) smoker",
  physact_low_vs_high_BIN = "Low (vs high) PA level",
  physact_mod_vs_high_BIN = "Moderate (vs high) PA level",
  alch_daily_rarelynever_BIN = "Daily (vs rarely/never) drink alcohol",
  alch_weekly_rarelynever_BIN = "Weekly (vs rarely/never) drink alcohol",
  alch_monthly_rarelynever_BIN = "Monthly (vs rarely/never) drink alcohol"
) ; #behav_covars.labs %>% saveRDS(paste0(outDir, "/behav_covars_labs.rda"))

diet_vars <- c("raw_veg_QT", "cooked_veg_QT", "fresh_fruit_QT", "dried_fruit_QT", 
               "oily_fish_QT", "nonoily_fish_QT", "procmeat_QT", "poultry_QT", 
               "cheese_QT", "beef_QT", "lamb_QT",  "pork_QT", "bread_type_white_vs_brown_or_whole_BIN",
               "bread_intake_QT", "milk_type_full_vs_low_or_nonfat_BIN",
               "cereal_type_sugar_vs_any_bran_BIN", "cereal_intake_QT", 
               "spread_type_butter_vs_any_other_BIN", "coffee_type_decaf_vs_regular_BIN",  
               "coffee_QT", "tea_QT", "water_QT", "addsalt_always_often_vs_nrs_BIN",  
               "hotdrink_temp_hot_or_vhot_vs_warm_BIN")


## Build Models --------------------

adj_life <- c("alch_freq.lab", "smoke_level.lab", "physact_level.lab")
adj_diet <- c(names(behav_covars.labs)[1:24])
adj_dietlife <- c(adj_life, adj_diet)

model_life <- paste0(models.nofast.l$BMI, "+", paste0(adj_life, collapse = "+"))
model_diet <- paste0(models.nofast.l$BMI, "+", paste0(adj_diet, collapse = "+"))
model_lifediet <- paste0(model_life, "+", paste0(adj_diet, collapse = "+"))

adj_models <- list("BMI Model\n(Primary)"=models.nofast.l$BMI, "Lifestyle\nAdjusted"=model_life, 
                   "Diet\nAdjusted"=model_diet, "Lifestyle+Diet\nAdjusted"=model_lifediet,
                   "Lifestyle+Diet+\nTotal Energy Adjusted"=paste0(model_lifediet,"+kcal_plaus"))


# =============================================================================
## Check quality of behavioral covars for capturing variability in 0-2 hr glu
# =============================================================================

library(lmtest)
do.call(rbind.data.frame, lapply(2:length(adj_models), function(m) {
  complete <- analysis %>% select("glu2hr_mgdl", all_of(strsplit(adj_models[[m]], split="[+]")[[1]])) %>% 
    filter(complete.cases(.))
  mRed <- lm(formula(paste0("glu2hr_mgdl~", adj_models[[1]])), data=complete)
  mFull <- lm(formula(paste0("glu2hr_mgdl~", adj_models[[m]] )), data=complete)
  lrt <- lmtest::lrtest(mRed, mFull) ; cbind.data.frame(Model=names(adj_models)[m], pLRT = lrt$`Pr(>Chisq)`[2])
  })) %>% 
  fwrite(paste0("../archive/data/processed/manuscript//tab_res_behav_covars_bmi_R1.csv"))


# ===========================================================
## Select lifestyle and dietary covariates for comparison
# ===========================================================

## Confirm TAS2R38 variant asssoc with behavioral covars ----------

do.call(rbind.data.frame, lapply(1:length(behav_covars.labs), function(d) {
  print_lm(exposure="taste_diplos.num", outcome=names(behav_covars.labs)[d], 
           covariates = models.nofast.l$BMI, label=behav_covars.labs[d])
})) %>% fwrite(paste0(outDir, "/tab_sens_behavcovars_diplo_bmi.csv"))
  

## Check associations of behavioral covars with 0-2 hr glucose ----------

do.call(rbind.data.frame, lapply(1:length(behav_covars.labs), function(d) {
  print_lm(outcome="glu2hr_mgdl", exposure=names(behav_covars.labs)[d],
           covariates = models.nofast.l$BMI, label=behav_covars.labs[d])
  })) %>% fwrite(paste0(outDir, "/tab_sens_glu2hr_behavcovars_bmi.csv"))


# ==========================================================================
## Run sequentially adjusted associations (with continuous TAS2R38 )
# ==========================================================================

do.call(rbind.data.frame, lapply(1:length(adj_models), function(m) {
  mBase <- summary(lm(formula(paste0("glu2hr_mgdl~", adj_models[[m]])),data=analysis))
  mTaste <- summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", adj_models[[m]] )), data=analysis))
  cbind.data.frame(Model=names(adj_models)[m], 
                   baseR2=mBase$adj.r.squared, tasteR2=mTaste$adj.r.squared,
                   beta=mTaste$coef[2,1], se=mTaste$coef[2,2], p=mTaste$coef[2,4] ) 
  })) %>% 
  fwrite(paste0("../archive/data/processed/manuscript/tab_res_glu_diplo_behav_R1.csv"))


################################################################################
## Negative Control Experiments
################################################################################

behav_covars_sel <- c("tea_QT", "coffee_QT", "addsalt_3lvl.lab", "alch_freq.lab",
                      "cereal_intake_QT", "spread_type_butter_vs_any_other_BIN", 
                      "bread_type_white_vs_brown_or_whole_BIN", "fresh_fruit_QT", 
                      "procmeat_QT")

# =============================================================
## Load & prepare additional bitter taste receptor variants
# =============================================================

# load dosage data for bitter snps
bitter_snps.df <- read.table("../archive/data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 

## Add variants to analysis dataframe
analysis <- analysis %>% left_join(bitter_snps.df, by="id")
bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")
analysis <- analysis %>% mutate(rs10772420_A = 2-rs10772420_G)


# =============================================================================
## R1: Genotype tables
# =============================================================================

table(analysis$rs2597979_G)
table(analysis$rs10772420_A)


# =============================================================================
## Associations of negative control variants with diet/lifestyle covariates
# =============================================================================

## Unadjusted P-trend tests (for Figure)
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(bitter_foods), function(v) {
    dat <- analysis %>% select(SNP=snp, var=bitter_foods[v])
    if(is.numeric(dat$var) ==T) {
      cbind.data.frame(Var=bitter_foods[[v]], Pvalue=summary(lm(var~SNP, data=dat))$coef[2,4])
    } else {
      cbind.data.frame(Var=bitter_foods[[v]], Pvalue=chisq.test(dat$SNP, dat$var)$p.value)
    }
  })) %>% mutate(P_formatted = format_p.fun(Pvalue)) %>% mutate(SNP=snp)
  }))  %>% fwrite(paste0(outDir, "tab_sens_dietlife_negcntrl_snps_unadj.csv"))


## BMI-adjusted models (for sensitivity)
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(behav_covars.labs), function(d) {
    print_lm(exposure=snp, outcome=names(behav_covars.labs)[d],
             covariates = models.nofast.l$BMI, label=behav_covars.labs[d]) }))
  })) %>% fwrite(paste0(outDir, "/tab_sens_dietlife_negcntrl_snps_bmi.csv"))



# ======================================================================
## Associations of TAS2R38 and Negative Control variants with 0-2hr glucose, 
## before and after adjusting for diet/lifestyle
# ======================================================================

## EMM in the fully-adjusted model
emm_bittersnps_glu2hr <- do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(adj_models), function(m) {
    get_emm.fun(exposure = snp, outcome = "glu2hr_mgdl", covars = adj_models[[m]], 
                reference = "0_0", data=analysis, set.rg.limit = 900000)$emm })) %>%
    mutate(Model = adj_models[m])
})) ; emm_bittersnps_glu2hr <- emm_bittersnps_glu2hr %>% mutate(
  Model = rep(rep(adj_models, each=3), length(bitter_snps_cat))
) %>% fwrite(paste0(outDir, "/tab_sens_emm_glu_negcntrl_snps.csv"))

## --> Double check required rg.limit 



################################################################################
## Sex-stratified exploratory analyses
################################################################################

sex.l <- list(c("Female", "Male"))

## NO significanat sex interactions
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*sex+", adj_models[[1]])), data=analysis))
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*sex+", adj_models[[4]])), data=analysis))


## Sex-stratified associations with 0-2hr glucose ------------------
# BMI model
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[1]]),
         label="Female", data=analysis %>% filter(sex=="Female"))
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[1]]),
         label="Male", data=analysis %>% filter(sex=="Male"))

# Diet+Lifestyle adjusted model ------------------
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[4]]),
         label="Female", data=analysis %>% filter(sex=="Female"))
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[4]]),
         label="Female", data=analysis %>% filter(sex=="Male"))


# ======================================================================
## Re-code TAS2R38 as PAV carriers to re-test association in females
# ======================================================================

analysis <- analysis %>% mutate(
  diplo_pav = factor(case_when(taste_diplos == "PAV/PAV" | taste_diplos == "AVI/PAV" ~ "PAV Carrier", 
                               taste_diplos == "AVI/AVI" ~ "AVI/AVI",
                               TRUE ~ as.character(NA)
  ), levels = c("AVI/AVI", "PAV Carrier")))


# BMI model ------------------

make_pretty_lm(print_lm(exposure = "diplo_pav", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[1]]),
         label="Female", data=analysis %>% filter(sex=="Female")))
anova(lm(formula(paste0("glu2hr_mgdl~diplo_pav+", gsub("[+]sex", "", adj_models[[1]]))), data=analysis %>% filter(sex=="Female")))

make_pretty_lm(print_lm(exposure = "diplo_pav", outcome="glu2hr_mgdl", covariates = gsub("[+]sex","", adj_models[[1]]),
         label="Male", data=analysis %>% filter(sex=="Male")))


## Compare r2 for additive vs. dominant model ----------

# Additive (numeric)
summary(lm(glu2hr_mgdl~taste_diplos.num, data=analysis %>% filter(sex=="Female")))
  # R2: 4.059e-05 ; Adj R2: 8.635e-06  ; F-statistic:  1.27; p-value: 0.2597

summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", gsub("[+]sex","", adj_models[[1]]) )), 
           data=analysis %>% filter(sex=="Female"))) 
  # Multiple R-squared:  0.05994,	Adjusted R-squared:  0.05895 
  # F-statistic: 60.41 on 33 and 31262 DF,  p-value: < 2.2e-16


# Dominant (PAV carriers vs AVI/AVI)
summary(lm(glu2hr_mgdl~diplo_pav, data=analysis %>% filter(sex=="Female")))
  # R2: 7.446e-05,	Adj R2: 4.251e-05 ; F-statistic: 2.33 ; p-value: 0.1269

summary(lm(formula(paste0("glu2hr_mgdl~diplo_pav+", gsub("[+]sex","", adj_models[[1]]) )), 
           data=analysis %>% filter(sex=="Female"))) 
  # R2: 7.446e-05,	Adj R2: 4.251e-05 ; F-statistic: 2.33 ; p-value: 0.1269


####################################################################################
## Added after first round of review: Age and BMI-stratified exploratory analyses
####################################################################################

# ================================
## Age: above/below median
# ================================

analysis <- analysis %>% mutate(age_gtmed = ifelse(age>=57,"gte_57", "lt_57"))

## NO significanat Age interactions (above/below median)
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*age_gtmed+", adj_models[[1]])), data=analysis))
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*age_gtmed+", adj_models[[4]])), data=analysis))


## Age-stratified associations with 0-2hr glucose ------------------
# BMI model: slight differences in significance but effect estimates in consistent direction & magnitude 
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = adj_models[[1]],
         label="Above 57 years", data=analysis %>% filter(age_gtmed=="gte_57")) ## P=0.02
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = adj_models[[1]],
         label="Below 57 years", data=analysis %>% filter(age_gtmed=="lt_57")) ## P=0.06

# Diet+Lifestyle adjusted model ------------------
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = adj_models[[4]],
         label="Above 57 years", data=analysis %>% filter(age_gtmed=="gte_57"))
print_lm(exposure = "taste_diplos", outcome="glu2hr_mgdl", covariates = adj_models[[4]],
         label="Below 57 years", data=analysis %>% filter(age_gtmed=="lt_57"))


# ------------------------------------------------
## Continuous age interaction (& by quantiles)
# ------------------------------------------------

anova(lm(glu2hr_mgdl~taste_diplos.num*age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10, data=analysis))
anova(lm(glu2hr_mgdl~taste_diplos.num*age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+bmi, data=analysis))


## Derive categorical variable for age

quantile(analysis$age) 
# 0%  25%  50%  75% 100% 
# 39   49   57   63   72 

analysis <- analysis %>% mutate(
  age_cat = case_when(age<=49~"Q1", age>49 & age<=57~"Q2", age>57 & age<=63~"Q3", age>63~"Q4")
) ; prop.table(table(analysis$age_cat))


# Categorical age interaction 
anova(lm(glu2hr_mgdl~taste_diplos.num*age_cat+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10, data=analysis))
anova(lm(glu2hr_mgdl~taste_diplos.num*age_cat+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+bmi, data=analysis))


# ================================================================
## Testing for differences in total energy by TAS2R38 diplotype
# ================================================================

# N with taste diplo & energy intake data
analysis %>% filter(is.na(taste_diplos)==F) %>% reframe(has_kcal=table(is.na(kcal_plaus))) # n = 129615 (% = 59.27)

# Descriptive analyses of energy intake and TAS2R38 diplotypes
tapply(analysis$kcal_plaus, analysis$taste_diplos, mean_sd, d=1)
summary(lm(kcal_plaus~taste_diplos.num, data=analysis))
summary(lm(kcal_plaus~taste_diplos.num+bmi, data=analysis))

# Differences in energy intake by TAS2R38 (with confounder adjustments)
rbind.data.frame(
  print_lm(exposure="taste_diplos", outcome="kcal_plaus", covariates="taste_diplos", label="Unadjusted"),
  print_lm(exposure="taste_diplos", outcome="kcal_plaus", covariates="bmi", label="BMI only"),
  print_lm(exposure="taste_diplos", outcome="kcal_plaus", covariates="sex", label="Age+sex"),
  print_lm(exposure="taste_diplos", outcome="kcal_plaus", covariates="age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10", label="Age+sex+gPCs"),
  print_lm(exposure="taste_diplos", outcome="kcal_plaus", covariates="age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+bmi", label="Age+sex+gPCs+BMI")
)


# ------------------------------------------------------------------------
## Check whether other bitter taste alleles associate with total energy?
# ------------------------------------------------------------------------

# Summary of genotype distributions
analysis %>% filter(is.na(taste_diplos)==F) %>%
  reframe(rs10772420_A=n_pct(rs10772420_A),
          rs2597979_G=n_pct(rs2597979_G),
          taste_diplos = n_pct(taste_diplos))

# restricted to those with 0-2 hr glucose data
analysis %>% filter(is.na(taste_diplos)==F) %>%
  filter(fasting_hrs <= 2)%>%
  reframe(rs10772420_A=n_pct(rs10772420_A),
          rs2597979_G=n_pct(rs2597979_G),
          taste_diplos = n_pct(taste_diplos))

tapply(analysis$kcal_plaus, analysis$rs2597979_G, mean_sd, d=1)
tapply(analysis$kcal_plaus, analysis$rs10772420_A, mean_sd, d=1)

# Summarise TAS2R38 & energy intake
lapply(c("rs2597979_G", "rs10772420_A"), function(snp) {
  rbind.data.frame(
    print_lm(exposure=snp, outcome="kcal_plaus", covariates="taste_diplos", label="Unadjusted"),
    print_lm(exposure=snp, outcome="kcal_plaus", covariates="age+sex", label="Age+sex"),
    print_lm(exposure=snp, outcome="kcal_plaus", covariates="age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10", label="Age+sex+gPCs")
  )
})

summary(lm(kcal_plaus~rs2597979_G, data=analysis)) #P=0.117
summary(lm(kcal_plaus~rs10772420_A, data=analysis)) #P=0.36


# ====================================================================
# Add adjustment for total energy in analyses of 0-2 h glucose
# ====================================================================

adj_models$EnergyAdj = paste0(model_lifediet,'+kcal_plaus')
do.call(rbind.data.frame, lapply(1:length(adj_models), function(m) {
  mBase <- summary(lm(formula(paste0("glu2hr_mgdl~", adj_models[[m]])),data=analysis))
  mTaste <- summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", adj_models[[m]] )), data=analysis))
  cbind.data.frame(Model=names(adj_models)[m], 
                   baseR2=mBase$adj.r.squared, tasteR2=mTaste$adj.r.squared,
                   beta=mTaste$coef[2,1], se=mTaste$coef[2,2], p=mTaste$coef[2,4]) 
  })) 

print_lm(exposure="taste_diplos", outcome="glu2hr_mgdl", covariates=adj_models$EnergyAdj, 
         label="Lifestyle+Diet+EnergyAdj") %>% make_pretty_lm()
print_lm(exposure="taste_diplos.num", outcome="glu2hr_mgdl", covariates=adj_models$EnergyAdj,
         label="Lifestyle+Diet+EnergyAdj") %>%
  make_pretty_lm()


##EOF

# Last Updated: Sept 3, 2025 


