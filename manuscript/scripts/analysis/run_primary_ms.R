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
analysis <- readRDS(paste0("../data/processed/ukb_analysis_vManuscript_EUR.rda"))

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
## Distributions of haplotyeps & diplotypes 
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
  write.csv(paste0("../data/processed/descr_tab_tastediplo_ac_",tag,".csv"))


## Demographic, behavioral & cardiometabolic risk factors 
vars_to_summarize <- c(
  age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking Status", cigarettes_per_day="Cigarettes per day (among smokers)", 
  physact_level.lab="Physical Activity Level", physact_met_excess="Physical Activity, excess MET/wk",
  alch_freq.lab="Alcohol Frequency", alch_1to4d.lab="Alcohol Frequency, 1-4 per week",
  alch_drinker_status.lab = "Alcohol Drinker", #alch_heavydrinker = "Heavy alcohol drinker",
  alch_drinks_per_week = "Drinks per week",
  educ_isced.lab="Education Level", income_level.lab = "Income Level",
  #sbp_adj="SBP, mmHg (med adjusted)", dbp="DBP, mmHg (med adjusted)",
  #tg="Triglyceride, mmol/L", ldl="LDL, mmol/L", hdl="HDL, mmol/L", 
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


#do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(fast) {
#  do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
#    print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = paste0(fast, "-Glucose (mg/dL)"), 
#             label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
#             data = analysis %>% filter(fast_cat_with_6to12 == fast)) })) })) %>% 
#  write.csv(paste0(outDir, "res_tab_lm_gluXfast12to24_tastediplo_mgdl_",tag,".csv"), row.names = F)


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

#  Outcome            Model Exposure     N              Beta_95CI P_t.test P_f.test P_trend.test
# glu_mgdl Glucose, >= 3 hr  AVI/AVI 54035                      -        -    0.383        0.629
# glu_mgdl Glucose, >= 3 hr  AVI/PAV 78481  0.025 (-0.075, 0.124)    0.629        -            -
# glu_mgdl Glucose, >= 3 hr  PAV/PAV 28520 -0.046 (-0.177, 0.084)    0.485        -            -
  
  
# Continuous
make_pretty_lm(print_lm(exposure="taste_diplos.num", outcome="glu_mgdl", covariates = models.l$BMI,
         label="Glucose, >= 3 hr", data = analysis %>% filter(fast_cat != "0to2hr")))

#  Outcome            Model         Exposure      N              Beta_95CI P_t.test
# glu_mgdl Glucose, >= 3 hr taste_diplos.num 161036 -0.016 (-0.079, 0.048)    0.629


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
pct_beta_change.fun(m_base$coef[2,1], m_bmi$coef[2,1]) #1.063%



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
  write.csv(paste0(outDir, "/tab_sens_emm_gluXfast_tas2r38snps_bmi.csv"), row.names = F)


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
)

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
                   "Diet\nAdjusted"=model_diet, "Lifestyle+Diet\nAdjusted"=model_lifediet)


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
  fwrite(paste0(outDir, "/tab_res_behav_covars_bmi.csv"))


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
  fwrite(paste0(outDir, "/tab_res_glu_diplo_behav.csv"))



################################################################################
## Negative Control Experiments
################################################################################

# =============================================================================
## Associations of negative control variants with diet/lifestyle covariates
# =============================================================================

do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(behav_covars.labs), function(d) {
    print_lm(exposure=snp, outcome=names(behav_covars.labs)[d],
             covariates = models.nofast.l$BMI, label=behav_covars.labs[d]) }))
  })) %>% fwrite(paste0(outDir, "/tab_sens_dietlife_negcntrl_snps_bmi.csv"))


# =============================================================
## Load & prepare additional bitter taste receptor variants
# =============================================================

# load dosage data for bitter snps
bitter_snps.df <- read.table("../data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 

## Add variants to analysis dataframe
analysis <- analysis %>% left_join(bitter_snps.df, by="id")
bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")
analysis <- analysis %>% mutate(rs10772420_A = 2-rs10772420_G)











## P-values for unadjusted P-trend tests (with continuous variants) 
negcontrol_pvals <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(bitter_foods), function(v) {
    dat <- analysis %>% select(SNP=snp, var=bitter_foods[v])
    if(is.numeric(dat$var) ==T) {
      cbind.data.frame(Var=bitter_foods[[v]], Pvalue=summary(lm(var~SNP, data=dat))$coef[2,4])
    } else {
      cbind.data.frame(Var=bitter_foods[[v]], Pvalue=chisq.test(dat$SNP, dat$var)$p.value)
    }
  })) %>% mutate(P_formatted = format_p.fun(Pvalue)) %>% mutate(SNP=snp)
}))


negcontrl_palette <- c(diplo_palette[[2]],rgb2hex(109,165,103), rgb2hex(193,121,99))
names(negcontrl_palette) <- snp_labels[c(1,4:5)]

p
##EOF


