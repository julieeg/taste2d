# Title: Sensitivity analysis - plausible energy intake 
# Date Updated: 09-05-2024

source("../scripts/pantry.R")

# create output directories
ANC="EUR"
outDir=paste0("../data/processed/")
tag="v6"


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  mutate(n_24hr.lab = ifelse(n_24hr==1, "Has 24HR", "No 24HR")) %>%
  mutate(n_24hr_plaus.lab = factor(ifelse(n_24hr_plaus, "Plausible 24HR", "No 24HR"), levels=c("Plausible 24HR", "No 24HR")),
         glu_mgdl = glu*18.018, glu2hr_mgdl = glu2hr*18.018)

#Set variables for workflow
n24HR.l <- c("Has 24HR", "No 24HR")
n24HR_plaus.l <- c("Plausible 24HR", "No 24HR")

# Load technical covariates 
tech <- readRDS("../data/processed/ukb_glucose_assays.rda") %>% select(id, ends_with(".0") | ends_with(".0.lab"))
analysis <- analysis %>% left_join(tech, by="id")# %>% 

dim(analysis)


# ====================
## Build Models
# ====================
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+ac.f+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level.lab+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC3+dietPC4+dietPC6+dietPC11+dietPC12+dietPC13+dietPC15+dietPC16+dietPC17")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l))
names(models.nofast.l) <- names(models.l)


# sensitivity models
ms_ses = paste0(m4_diet, "+income_level.lab+educ_isced.lab")
ms_diet_all = paste0(m3_lifestyle, "+", paste0("dietPC",1:24,collapse = "+"))
ms_f2c = paste0(m4_diet, "+FIB2CHO")
ms_f2c_ses = paste0(m4_diet, "+FIB2CHO+income_level.lab+educ_isced.lab")
ms_tech = paste0(m4_diet, "+gluB_aliquot.0+gluB_assay_date.0+dilution.0")
ms_altcov = gsub("alch_freq.lab", "alch_drinks_per_week", gsub("physact_level.lab", "physact_met_excess", m3_lifestyle))
ms_diet_altcov = gsub("alch_freq.lab", "alch_drinks_per_week", gsub("physact_level.lab", "physact_met_excess", m4_diet))

models.sensitivity.l = list("Diet.Patterns"=m4_diet, "SES"=ms_ses, "All.Diet.PCs"=ms_diet_all, "Fib2Carb"=ms_f2c, "Fib2Carb+SES"=ms_f2c_ses, 
                            "Technical"=ms_tech, "Alt.Covariates"=ms_altcov, "Diet+Alt.Covariates"=ms_diet_altcov)
models.sensitivity.nofast.l <- as.list(gsub("[+]fast_cat", "", models.sensitivity.l)) ; names(models.sensitivity.nofast.l) <- names(models.sensitivity.l)


## List of descriptive variables 
vars_to_summarize_abbrev <- c(
  age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking Status", physact_level.lab="Physical Activity Level", 
  alch_freq.lab="Alcohol Frequency", alch_drinks_per_week = "Drinks per week",
  educ_isced.lab="Education Level", income_level.lab = "Income Level",
  hba1c_max="HbA1C, %", glu="Glucose, mmol/L", glu2hr="0-2hr Glucose, mmol/L"
)

#Glycemic
glycemic <- c(glu="Glucose", glu2hr="0-2hr Glucose", hba1c="HbA1c")
glycemic_vars.l <- list("Glucose"="glu", "0-2hr Glucose"="glu2hr", "HbA1c"="hba1c")
glycemic_vars_mgdl.l <- list("Glucose (mg/dL)"="glu_mgdl", "0-2hr Glucose (mg/dL)"="glu2hr_mgdl", "HbA1c (%)"="hba1c")



###################################################################
###  Check that covariates are capturing the expected variance  ###
## Suggested by JC (9/17)
###################################################################

lifestyle_ses <- c(smoke_level.lab="Smoking", alch_freq.lab="Alcohol", physact_level.lab="PA Level",
                   income_level.lab="Income Level", educ_isced.lab.f="Education Level")

# add factor variable for education
analysis <- analysis %>% mutate(
  educ_isced.lab.f = factor(educ_isced.lab, levels = c(paste0("Level ",1:5))))

## Summarize lm associations of lifestyle/ses vars with glycemic outcomes
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(lifestyle_ses), function(v) { 
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {cov<-models.nofast.l$BMI ; dat<-analysis %>% filter(fast_cat=="0to2hr")} else {
      cov <- models.l$BMI ; dat <- analysis} ; return(
        print_lm(exposure = names(lifestyle_ses)[v], outcome = glycemic_vars_mgdl.l[g], label.outcome=names(glycemic_vars_mgdl.l)[g],
                 covariates = cov, label=paste0("BMI_",lifestyle_ses[v]), data=dat, digits=c(1,1,6)) ) }))  })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_lifestyle_ses_mBMI_mgdl_",tag,".csv"), row.names = F)


## diet traits & 2hr glucose
diet_vars <- c(dietPC3="Diet PC3", dietPC4 = "Diet PC4", dietPC6="Diet PC6", 
               dietPC11="Diet PC11", dietPC12="Diet PC12", dietPC13 = "Diet PC13",
               dietPC15="Diet PC15", dietPC16="Diet PC16", dietPC17 = "Diet PC17")


## Summarize lm associations of lifestyle/ses vars with glycemic outcomes
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(diet_vars), function(v) { 
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {cov<-models.nofast.l$BMI ; dat<-analysis %>% filter(fast_cat=="0to2hr")} else {
      cov <- models.l$BMI ; dat <- analysis} ; return(
        print_lm(exposure = names(diet_vars)[v], outcome = glycemic_vars_mgdl.l[g], label.outcome=names(glycemic_vars_mgdl.l)[g],
                 covariates = cov, label=paste0("BMI_",diet_vars[v]), data=dat) ) }))  })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_covars_dietPCs_mgdl_",tag,".csv"), row.names = F)


do.call(rbind.data.frame, lapply(1:length(diet_vars), function(d) { 
   print_lm(exposure = names(diet_vars)[d], outcome = glycemic_vars_mgdl.l[[2]], label.outcome=names(glycemic_vars_mgdl.l)[2],
                 covariates = models.nofast.l$BMI, label=paste0("BMI_", diet_vars[d]), data=dat) })) %>% 
  #write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_lifestyle_ses_mBMI_mgdl_",tag,".csv"), row.names = F)


## Continuous (more powered?/precise?) covariates
lifestyle_ses.cont <- c(
  cigarettes_per_day="Cigarettes per day", alch_drinks_per_week="Drinks per week",
  physact_met_excess="Physical Activity, eMET/wk", educ_years="Education years")

## summarize associations of continuous (more precise) covariates with glycemic outcomes
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(lifestyle_ses.cont), function(v) { 
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {cov<-models.nofast.l$BMI ; dat<-analysis %>% filter(fast_cat=="0to2hr")} else {
      cov <- models.l$BMI ; dat <- analysis} ; return(
        print_lm(exposure = names(lifestyle_ses.cont)[v], outcome = glycemic_vars_mgdl.l[[g]], 
                 label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = cov, data=dat,
                 label=paste0("BMI_",lifestyle_ses[v]), lm_trend = T)) })) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_lifestyle_ses_cont_mBMI_mgdl_",tag,".csv"), row.names = F)


# Calculate adjusted R-squared for each outcome 
sapply(1:length(lifestyle_ses), function(v) {
  summary(lm(formula(paste0("hba1c~", names(lifestyle_ses)[v])),data=analysis))$adj.r.squared 
}) ; sapply(1:length(lifestyle_ses), function(v) {
  summary(lm(formula(paste0("glu2hr~", names(lifestyle_ses)[v])),data=analysis))$adj.r.squared 
}) ; sapply(1:length(lifestyle_ses), function(v) {
  summary(lm(formula(paste0("glu~fasting_hrs+", names(lifestyle_ses)[v])),data=analysis))$adj.r.squared 
})

## =======================================
# Adjust for More "precise" covariates
## =======================================

do.call(rbind.data.frame, lapply(1:3, function(g) {
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {models.l.use<-models.sensitivity.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr")} else {
    models.l.use <- models.sensitivity.l ; dat.use <- analysis } ; return(
  do.call(rbind.data.frame, lapply(c(7,8), function(v) { 
        print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], label.outcome = names(glycemic_vars_mgdl.l)[[g]],
                 covariates = models.l.use[[v]], label=names(models.l.use)[v], lm_trend = T, data=dat.use) })) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_mAltCovars_mDietAltCovars_mgdl_",tag,".csv"), row.names = F)


# ============================================
## Stratified analysis based diet extremes
# ============================================
dpc1_median=quantile(analysis$dietPC1, probs=seq(0,1,0.5), include.lowest=F)[-1]
dpc1_tertiles=quantile(analysis$dietPC1, probs=seq(0,1,0.33), include.lowest=F)[-1]

## Based on dietPC1 (generally healthy/unhealthy indicators)
analysis <- analysis %>% mutate(
  dpc1_med = factor(ifelse(dietPC1 >= dpc1_median[1],"Above","Below"), levels=c("Below", "Above")),
  dpc1_tert = factor(case_when(dietPC1<dpc1_tertiles[1] ~ "Low",
                          dietPC1>=dpc1_tertiles[1] & dietPC1<dpc1_tertiles[2] ~ "Moderate",
                          dietPC1>dpc1_tertiles[2] ~ "High"), levels=c("Low", "Moderate", "High"))
)

## Run models in each "extreme"
print_summary_table(data=analysis, vars_to_summarize = vars_to_summarize_abbrev, 
                    var_strata = "dpc1_med", var_strata_order = c("Below", "Above"), digits = c(3,1,10)) %>%
    write.csv(paste0("../data/processed/sens_tab_descr_dPC1median_",tag,".csv"))

# In dietPC median (above vs below)
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(c("Below", "Above"), function(i) {
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {models.l.use<-models.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr")} else {
      models.l.use <- models.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
                   label.outcome = paste0(i, " dietPC1.", names(glycemic_vars_mgdl.l)[g]),
                   covariates = paste0(models.l.use[[v]], "+income_level.lab+educ_isced.lab"), data=dat.use %>% filter(dpc1_med == i),
                   label=names(models.l.use)[v]) })) ) })) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_dPC1median_mSES_mgdl_",tag,".csv"), row.names = F)

# In dietPC tertiles
print_summary_table(data=analysis, vars_to_summarize = vars_to_summarize_abbrev, 
                    var_strata = "dpc1_tert", var_strata_order = c("Low", "Moderate","High"), digits = c(3,1,10)) %>%
  write.csv(paste0("../data/processed/sens_tab_descr_dPC1tertile_",tag,".csv"))

do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(c("Low", "Moderate", "High"), function(i) {
    if(glycemic_vars.l[[g]]=="glu2hr") {models.l.use<-models.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr")} else {
      models.l.use <- models.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos", outcome = glycemic_vars.l[g], 
                   label.outcome = paste0(i, " dietPC1.", names(glycemic_vars.l)[g]),
                   covariates = models.l.use[[v]], data=dat.use %>% filter(dpc1_tert == i),
                   label=names(models.l.use)[v], lm_trend = F) })) ) })) })) %>% 
    write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_dPC1tertile_",tag,".csv"), row.names = F)


############################################
###  Adjust for continuous fasting time  ###
############################################

models.fasthrs.l <- lapply(models.l, function(m) gsub("fast_cat", "fasting_hrs", m))

do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(models.fasthrs.l), function(v) { 
        print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
                 label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = models.fasthrs.l[[v]], 
                 label=names(models.fasthrs.l)[v], lm_trend = T, data=analysis) })) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_mFastHrs_mgdl_",tag,".csv"), row.names = F)


## =========================================================
# Restrict to 0-12 hours; and toggle fasting categories?
## =========================================================

do.call(rbind.data.frame, lapply(1:length(models.l), function(v) { 
  print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l$Glucose, 
           label.outcome = "Glucose (mg/dL)", covariates = models.l[[v]], 
           label=names(models.l)[v], lm_trend = T, data=analysis %>% 
             filter(fasting_hrs <= 12 ) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_glu_tastediplos_mFastHrs_mgdl_",tag,".csv"), row.names = F)


# ===================
## 1-hr Glucose
# ===================

do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome="1hr Glucose (mg/dL)",
           label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
           data = analysis %>% filter(fasting_hrs == 1)) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_glu1hr_tastediplos_mgdl_",tag,".csv"))


########################
###  Adjust for SES  ###
########################
rbind.data.frame(
  print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = "Glucose", label="SES", 
           covariates = models.sensitivity.l$SES, lm_trend = T, data = analysis),
  print_lm(exposure = "taste_diplos", outcome = "glu2hr_mgdl", label.outcome = "0-2hr Glucose", label="SES", 
           covariates = models.sensitivity.nofast.l$SES, lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")),
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label.outcome="HbA1c", label="SES", 
           covariates = models.sensitivity.l$SES, lm_trend = T, data = analysis)) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_rg_2hg_tastediplos_mSES_mgdl_",tag,".csv"), row.names = T)


####################################
###  Adjust for all 24 diet PCs  ###
####################################
# Categorical diplotypes 
rbind.data.frame(
  print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = "Glucose (mg/dL)", label="All.Diet.PCs", 
           covariates = models.sensitivity.l$All.Diet.PCs, lm_trend = T, data = analysis),
  print_lm(exposure = "taste_diplos", outcome = "glu2hr_mgdl", label.outcome = "0-2hr Glucose (mg/dL)", label="All.Diet.PCs", 
           covariates = models.sensitivity.nofast.l$All.Diet.PCs, lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")),
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label.outcome="HbA1c", label="All.Diet.PCs", 
           covariates = models.sensitivity.l$All.Diet.PCs, lm_trend = T, data = analysis)) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_rg_2hg_tastediplos_mAllDietPC_mgdl_",tag,".csv"), row.names = T)


# continuous diplotypes 
rbind.data.frame(
  print_lm(exposure = "taste_diplos.num", outcome = "glu_mgdl", label.outcome = "Glucose (mg/dL)", label="All.Diet.PCs", 
           covariates = models.sensitivity.l$All.Diet.PCs, lm_trend = T, data = analysis),
  print_lm(exposure = "taste_diplos.num", outcome = "glu2hr_mgdl", label.outcome = "0-2hr Glucose (mg/dL)", label="All.Diet.PCs", 
           covariates = models.sensitivity.nofast.l$All.Diet.PCs, lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")),
  print_lm(exposure = "taste_diplos.num", outcome = "hba1c", label.outcome="HbA1c", label="All.Diet.PCs", 
           covariates = models.sensitivity.l$All.Diet.PCs, lm_trend = T, data = analysis)) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_num_mAllDietPC_mgdl_",tag,".csv"), row.names = T)


## % change in Betas
print_lm(exposure = "taste_diplos", outcome = "glu2hr_mgdl", label="0-2hr.Glucose ~ All Diet Patterns", 
         covariates=models.sensitivity.nofast.l$Diet.Patterns, lm_trend = T, 
         data = analysis %>% filter(fast_cat=="0to2hr"))

print_lm(exposure = "taste_diplos", outcome = "glu2hr_mgdl", label="0-2hr.Glucose ~ All Diet Patterns", 
         covariates=models.sensitivity.nofast.l$All.Diet.PCs, lm_trend = T, 
         data = analysis %>% filter(fast_cat=="0to2hr"))


#############################################
###  Plausible 24HR + Fib2CHO adjustment  ###
#############################################

# ===================================================================
## Subsetting to individuals with **available 24HR data** to evaluate
## TAS2R38-diet associations; and additionally adjust for fiber/CHO
## Only include individuals with plausible enrgy: >= 600 & 
##    -for M: <20000 kJ (4780.2 kcal) 
##    -for F: <18000 kj (4302.1 kcal)
# ===================================================================

## describe subsample ==============
decr_vars <- c(
  taste_diplos = "TAS2R38 diplotype", age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking", physact_level="Physical Activity, MET/wk", alch_freq.lab="Alcohol Frequency", 
  alch_drinks_per_week="Drinks per week", educ_isced.lab="Education Level", income_level.lab = "Income Level")
descr_vars.f <-  c("taste_diplos", "smoke_level.lab", "physact_level",
                   "alch_freq.lab", "educ_isced.lab", "income_level.lab")

print_summary_table(data=analysis, var_strata = "n_24hr_plaus.lab", vars_to_summarize = decr_vars,
                    var_strata_order = n24HR_plaus.l, factor_vars = descr_vars.f, p_print=T) %>% 
  write.csv(paste0("../data/processed/sens_descr_tab_24HRplaus_",tag,".csv"))

print_summary_table(data=analysis, var_strata = "n_24hr.lab", var_strata_order = n24HR.l,
                    vars_to_summarize = c(ac.f="Assessment Center"), factor_vars = c("ac.f"), p_print=T) %>% 
  write.csv(paste0("../data/processed/sens_descr_tab_24HR_center_",tag,".csv"))

## subsample without 24HR data, for comparison ----  
models.sensitivity.f2c.l = models.sensitivity.l[c(1:2,4:5)]
models.sensitivity.f2c.nofast.l = models.sensitivity.nofast.l[c(1:2,4:5)]

models.sensitivity.nof2c.l = as.list(gsub("[+]FIB2CHO", "", models.sensitivity.f2c.l[1:2]))
names(models.sensitivity.nof2c.l)=names(models.sensitivity.f2c.l[1:2])
models.sensitivity.nof2c.nofast.l = as.list(gsub("[+]FIB2CHO", "", models.sensitivity.f2c.nofast.l[1:2]))
names(models.sensitivity.nof2c.nofast.l)=names(models.sensitivity.f2c.nofast.l[1:2])

## Add F2C adjustment in subsample with 24HR  ----
do.call(rbind.data.frame, lapply(1:3, function(g) {
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {
    models.l.use<-models.sensitivity.f2c.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr") } else {
    models.l.use <- models.sensitivity.f2c.l ; dat.use <- analysis} ; return(
      do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
        print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
             label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = models.l.use[[v]], 
             label=names(models.l.use)[v], lm_trend = T, data=analysis %>% filter(n_24hr.lab == "Has 24HR")) })) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplo_mF2C_mgdl_",tag,".csv"), row.names = F)

# numeric taste diplotype
do.call(rbind.data.frame, lapply(1:3, function(g) {
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {
    models.l.use<-models.sensitivity.f2c.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr") } else {
      models.l.use <- models.sensitivity.f2c.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos.num", outcome = glycemic_vars_mgdl.l[[g]], 
                   label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = models.l.use[[v]], 
                   label=names(models.l.use)[v], lm_trend = T, data=analysis %>% filter(n_24hr.lab == "Has 24HR")) })) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplo_num_mF2C_mgdl_",tag,".csv"), row.names = F)


## Add F2C adjusting in subsample without 24HR  ----
do.call(rbind.data.frame, lapply(1:3, function(g) {
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {
    models.l.use<-models.sensitivity.nof2c.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr") } else {
      models.l.use <- models.sensitivity.nof2c.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
                   label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = models.l.use[[v]], 
                   label=names(models.l.use)[v], lm_trend = T, data=analysis %>% filter(n_24hr.lab == "No 24HR")) })) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_rg_2hg_tastediplos_mNoF2C_mgdl_",tag,".csv"), row.names = F)

# numeric taste diplotype
do.call(rbind.data.frame, lapply(1:3, function(g) {
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {
    models.l.use<-models.sensitivity.nof2c.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr") } else {
      models.l.use <- models.sensitivity.nof2c.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos.num", outcome = glycemic_vars_mgdl.l[[g]], 
                   label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = models.l.use[[v]], 
                   label=names(models.l.use)[v], lm_trend = T, data=analysis %>% filter(n_24hr.lab == "No 24HR")) })) ) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_rg_2hg_tastediplos_num_mNoF2C_mgdl_",tag,".csv"), row.names = F)


###############################################################
##  Alternative variants for bitter taste: quinine/caffeine  ## 
###############################################################

# load dosage data for bitter snps
bitter_snps.df <- read.table("../data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 

## Add OTHER bitter SNPs
analysis <- analysis %>% left_join(bitter_snps.df, by="id")
bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")
analysis <- analysis %>% mutate(rs10772420_A = 2-rs10772420_G)

# categorical, per dominant allele 
analysis <- analysis %>% mutate(
  rs713598.a = descr_label_ordered.fun(., "rs713598_G", c("CC"=0, "CG"=1, "GG"=2)),
  rs1726866.a = descr_label_ordered.fun(., "rs1726866_G", c("AA"=0, "AG"=1, "GG"=2)),
  rs10246939.a = descr_label_ordered.fun(., "rs10246939_C", c("TT"=0, "TC"=1, "CC"=2)),
  rs2597979.a = descr_label_ordered.fun(., "rs2597979_G", c("AA"=0, "AG"=1, "GG"=2)),
  rs10772420.a = descr_label_ordered.fun(., "rs10772420_A", c("GG"=0, "GA"=1, "AA"=2))
) ; snps.a <- c("rs713598.a", "rs1726866.a", "rs10246939.a", "rs2597979.a","rs10772420.a")

# ===========================================================
## Correlation of bitter SNPs with bitter diet traits & PCs
# ===========================================================

bind_rows(
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    rbind.data.frame(
      print_lm(exposure = snp, outcome="coffee_QT", covariates = "age+sex+ac.f", label="Coffee", data=analysis),
      print_lm(exposure = snp, outcome="tea_QT", covariates = "age+sex+ac.f", label="Tea", data=analysis),
      print_lm(exposure = snp, outcome="raw_veg", covariates = "age+sex+ac.f", label="Raw vegetables", data=analysis),
      print_glm(exposure = snp, outcome="alch_heavydrinker", covariates = "age+sex+ac.f", label="Heavy drinker", data=analysis),
      print_lm(exposure = snp, outcome="addsalt_freq_QT", covariates = "age+sex+ac.f", label="Add Salt (cont)", data=analysis) ) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })),
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    do.call(rbind.data.frame, lapply(dietPCs, function(PC) {
      print_lm(exposure = snp, outcome = PC, covariates = c("age", "sex", "ac.f"), label = PC, data=analysis) })) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })) 
) %>% write.csv(paste0("../data/processed/sens_descr_tab_bitterDiet_bitterSNPs_",tag,".csv"), row.names=F)


# ===========================================================
## Primary analysis with RG, A1c & Glu x fasting
# ===========================================================

## Bitter SNPs & RG -----------------------
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {models.l.use<-models.nofast.l ; dat.use<-analysis %>% filter(fast_cat=="0to2hr")} else {
      models.l.use <- models.l ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = snp, outcome = glycemic_vars_mgdl.l[g], label.outcome = names(glycemic_vars_mgdl.l)[g],
                   covariates = models.l.use[[v]], label=names(models.l.use)[v], data=dat.use) })) ) 
    })) })) %>%
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_bittersnps_mgdl_",tag,".csv"), row.names = F)

#0-2hr glucose  
do.call(rbind.data.frame, lapply(snps.a, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu2hr_mgdl", covariates = models.nofast.l[[i]], 
             label=names(models.nofast.l)[i], label.outcome = "0-2hr Glucose (mg/dl)",
             data = analysis %>% filter(fast_cat == "0to2hr"), lm_trend = T)  })) %>%
    mutate(SNP=snp, .before=outcome) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_2hg_bittersnps_add_mgdl_",tag,".csv"), row.names = T)


# 6+hrs (~fasting glucose)
do.call(rbind.data.frame, lapply(snps.a, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu_mgdl", covariates = models.nofast.l[[i]], 
             label=names(models.nofast.l)[i], label.outcome = snp,
             data = analysis %>% filter(fasting_hrs >= 6), lm_trend = T)  })) })) %>%
  filter(model=="BMI")


## Bitter SNPs & Glu x fasting time -----------------------
fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l))
names(models.nofast.l) <- names(models.l)

## BMI Model: Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
      as.data.frame(print_lm(exposure = snp, outcome = "glu_mgdl", covariates = models.nofast.l[[i]], 
                             label=names(models.nofast.l)[i], label.outcome = paste0(fast," Glucose (mg/dL)"), 
                             data=analysis %>% filter(fast_cat == fast))) })) %>% 
      mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) })) %>%
  write.csv(paste0("../data/processed/sens_tab_lm_gluXfast_bittersnps_mBMI_mgdl_",tag,".csv"), row.names = F)




## Diet Model: Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
      as.data.frame(print_lm(exposure = snp, outcome = "glu_mgdl", covariates = models.nofast.l[[i]], 
                             label=names(models.nofast.l)[i], label.outcome = paste0(fast," Glucose (mg/dL)"), 
                             data=analysis %>% filter(fast_cat == fast))) })) %>% 
      mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) })) %>%
  write.csv(paste0("../data/processed/sens_tab_lm_gluXfast_bittersnps_mDiet_mgdl_",tag,".csv"), row.names = F)


# ==========================================================================
## Estimated marginal mean: bitter snps & Glu x fasting time (diet model)
# =========================================================================

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
    rs10246939_C == 2 ~ "1_1"
  ),
  rs10772420_A.cat = case_when(
    rs10772420_A == 0 ~ "0_0",
    rs10772420_A == 1 ~ "0_1",
    rs10772420_A == 2 ~ "1_1"
  ),
  rs2597979_G.cat = case_when(
    rs2597979_G == 0 ~ "0_0",
    rs2597979_G == 1 ~ "0_1",
    rs2597979_G == 2 ~ "1_1"
  )
)

bitter_snps_cat <- paste0(bitter_snps, ".cat")

# BMI
do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    get_emm.fun(exposure = snp, outcome = "glu_mgdl", covars = models.nofast.l$BMI, 
                reference = "0_0", data=analysis %>% filter(fast_cat == fast))$emm })) %>%
    mutate(fast = rep(fast_cat, each=3), .before=emmean) })) %>%
  write.csv(paste0("../data/processed/sens_tab_emm_gluXfast_bittersnps_mBMI_mgdl_",tag,".csv"), row.names = F)

# Diet Patterns
do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    get_emm.fun(exposure = snp, outcome = "glu_mgdl", covars = models.nofast.l$Diet.Patterns, 
                reference = "0_0",
                data=analysis %>% filter(fast_cat == fast))$emm })) %>%
    mutate(fast = rep(fast_cat, each=3), .before=emmean) })) %>%
  write.csv(paste0("../data/processed/sens_tab_emm_gluXfast_bittersnps_mDiet_mgdl_",tag,".csv"), row.names = F)

do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu2hr_mgdl", covariates = models.nofast.l[[i]], label=paste0(snp, "_", names(models.nofast.l)[i]), 
             data = analysis %>% filter(fast_cat == "0to2hr"), lm_trend = F)  })) 
  })) %>% write.csv(paste0("../data/processed/sens_tab_lm_2hg_bittersnps_cat_mgdl_",tag,".csv"), row.names = T)



#############################################
## Independent Replication & Meta-Analyses ##
#############################################

library(meta)
tas2r38_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C")

# ------------------------------------------------------
## Summarize TAS2R38 SNP - RG/2hG associations & save 
# ------------------------------------------------------

tas2r38 <- c("taste_diplos.num", tas2r38_snps)

meta_snps_pp_glu_ukbb <- matrix(NA, 12, 6)

do.call(rbind.data.frame, lapply(tas2r38, function(g) {
  do.call(rbind.data.frame, lapply(1:2, function(m) {
    rbind.data.frame(
      print_lm(exposure = g, outcome = "glu", covariates = models.nofast.l[[m]],
             label=paste0("0-2hr ~ ", g, " + ", names(models.nofast.l)[m]), 
             data=analysis %>% filter(fast_cat == "0to2hr")),
      print_lm(exposure = g, outcome = "glu", covariates = models.nofast.l[[m]],
               label=paste0("2hr ~ ", g, " + ", names(models.nofast.l)[m]), 
               data=analysis %>% filter(fasting_hrs == "2"))
    ) }))
}))


for(i in 1:length(tas2r38_snps)) {
  
  rg <- summary(lm(formula(paste0("glu ~ ", tas2r38_snps[i], "+", models.l$Base)), data = analysis))
  rgBMI <- summary(lm(formula(paste0("glu ~ ", tas2r38_snps[i], "+", models.l$BMI)), data = analysis))
  
  g2hr <- summary(lm(formula(paste0("glu ~ ", tas2r38_snps[i], "+", models.nofast.l$Base)), 
                    data = analysis %>% filter(fast_cat == "0to2hr")))
  g2hrBMI <- summary(lm(formula(paste0("glu ~ ", tas2r38_snps[i], "+", models.nofast.l$BMI)), 
                       data = analysis %>% filter(fast_cat == "0to2hr")))
  
  meta_snps_pp_glu_ukbb[i,1:6] <- c("UKBB_RG", tas2r38_snps[i], rbind(round(rg$coef[2,c(1:2,4)], 4)), length(rg$residuals))
  meta_snps_pp_glu_ukbb[i+3,1:6] <- c("UKBB_RG_AdjBMI", tas2r38_snps[i], rbind(round(rgBMI$coef[2,c(1:2,4)], 4)), length(rgBMI$residuals))
  meta_snps_pp_glu_ukbb[i+6,1:6] <- c("UKBB_02hG", tas2r38_snps[i],  rbind(round(g2hr$coef[2,c(1:2,4)], 4)), length(g2hr$residuals))
  meta_snps_pp_glu_ukbb[i+9,1:6] <- c("UKBB_02hG_AdjBMI" , tas2r38_snps[i],  rbind(round(g2hrBMI$coef[2,c(1:2,4)], 4)), length(g2hrBMI$residuals))
} ; meta_snps_pp_glu_ukbb <- as.data.frame(meta_snps_pp_glu_ukbb)

colnames(meta_snps_pp_glu_ukbb) <- c("Study", "SNP", "BETA", "SE", "P", "N")
meta_snps_pp_glu_ukbb

# ------------------------------------------------------
## Run IVW on MAGIC summary stats from
## Saxena et al., Nature genetics 2010;42;2;142-8
## Trait: 2-hr glucose from OGTT adjusted for BMI 
## plus (age, sex, & study-specific covariates)
# ------------------------------------------------------

meta_2hg<-read.csv("../data/meta_analysis_2hg.csv") %>%
  mutate(Variance=(SE^2), weights = 1/(SE^2)) %>% 
  mutate(W = c(1/Variance)) %>%
  mutate(IVW=BETA*W) %>%
  mutate(weights_pct = round(weights/sum(weights)*100, 1)) %>%
  filter(STUDY != "UKBB_0to2hG") %>%
  filter(STUDY != "Saxena_2010_MAGIC_2hG_OGTT_AdjBMI")

# Make data.frame with MAGIC meta-analysis results
mr_tas2r38_glu <- cbind.data.frame(
    Study=rep("Saxena_MAGIC_2hG_AdjBMI", 3), 
    SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
    Beta=c(-0.0790, -0.0510, -0.0570), SE=c(0.022, 0.021, 0.021), P_value=c(0.0002571, 0.01228, 0.006042)
    ) %>%
  mutate(Beta_mgdl=Beta*18.018, SE_mgdl=SE*18.018)

## Add Variance, Weights & Inverse Weights
mr_tas2r38_glu <- mr_tas2r38_glu %>% mutate(
  Variance=(SE_mgdl^2)) %>%
  mutate(W = c(1/Variance)) %>%
  mutate(IVW=Beta_mgdl*W) %>%
  mutate(weights = 1/(SE^2)) %>% 
  mutate(weights_pct = round(weights/sum(weights)*100, 1)) %>%
  mutate(SNP_alleles=c("rs713598 (C>G)", "rs1726866 (G>A)", "rs10246939 (T>C)")) %>%
  mutate(Bitter_Allele=c("G", "G", "C")) %>%
  mutate(Diplotype=1:3)


## Calculate Pooled Beta & SE  
m.gen <- metagen(TE = Beta_mgdl, seTE=SE_mgdl, studlab = SNP_alleles, method.random.ci = "classic",
                 fixed=TRUE, random = FALSE, sm= "SMD", method.tau = "REML", data=mr_tas2r38_glu)
z=-5.02
pnorm(q=z, mean=0, sd=1)

pdf("../output/meta_bittersnps_02hr_magic_v3.pdf", height = 3, width=8)
forest(m.gen, layout = "RevMan5",
       sortvar = Diplotype, 
       #hetstat = F, 
       digits=2, 
       at=seq(-2,2,1),
       colgap.forest=c("10 mm"),
       col.square =  "#96A0B3", col.square.lines = "#96A0B3", col.inside = "black",
       col.diamond = "#435269", col.diamond.lines = "#435269",
       #leftcols = c("SNP", "Bitter_Allele",  "w.fixed", "effect.ci"),
       leftcols = c("SNP", "Bitter_Allele", "effect.ci"),
       #leftlabs = c("Variant (Bitter\nTaste Allele)", "Taster\nAllele", "Wt\n(%)", "IVW Beta [95% CI]\nper Taster Allele"),
       leftlabs = c("Variant (Bitter\nTaste Allele)", "Taster\nAllele", "IVW Beta [95% CI]\nper Taster Allele"),
       smlab = "2-hr Glucose (mg/dL)\nadjusted for BMI", text.fixed = "Pooled Effect Estimate",
       just.addcols.left = "left", just.studlab = "left", just = "left",
       fontsize = 10,
       squaresize = 0.75, lwd.diamond = 5, lwd.square = 1.25, 
       colgap.left = c("5 mm"))
       # fs.axis = 6, fs.heading = 8, fs.smlab = 8,
dev.off()


## meta-analysis of Random Glucose
mr_tas2r38_rg <- cbind.data.frame(
  Study=rep("Lagou_2023_MAGIC_RG_MA_EUR_BMIadj", 3), 
  SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
  Beta=c(-0.000678723, -0.000684952, -0.000698865), SE=c(0.000299836, 0.000295672, 0.000295632), P_value=c(0.0236, 0.02053, 0.01808)) %>%
  mutate(Beta_mgdl=Beta*18.018, SE_mgdl=SE*18.018) %>%
  mutate(Variance=(SE_mgdl^2)) %>%
  mutate(W = c(1/Variance)) %>%
  mutate(IVW=Beta_mgdl*W) %>%
  mutate(weights = 1/(SE^2)) %>% 
  mutate(weights_pct = round(weights/sum(weights)*100, 1)) %>%
  mutate(SNP_alleles=c("rs713598 (C>G)", "rs1726866 (G>A)", "rs10246939 (T>C)")) %>%
  mutate(Bitter_Allele=c("G", "G", "C")) %>%
  mutate(Diplotype=1:3)

m.gen_rg <- metagen(TE = Beta_mgdl, seTE=SE_mgdl, studlab = SNP_alleles, method.random.ci = "classic",
                 fixed=TRUE, random = FALSE, sm= "SMD", method.tau = "REML", data=mr_tas2r38_rg)

pdf("../output/meta_bittersnps_rg_magic_v3.pdf", height = 3, width=8)
forest(m.gen, layout = "RevMan5",
       sortvar = Diplotype, 
       #hetstat = F, 
       digits=2, 
       at=seq(-2,2,1),
       colgap.forest=c("10 mm"),
       col.square =  "#96A0B3", col.square.lines = "#96A0B3", col.inside = "black",
       col.diamond = "#435269", col.diamond.lines = "#435269",
       #leftcols = c("SNP", "Bitter_Allele",  "w.fixed", "effect.ci"),
       leftcols = c("SNP", "Bitter_Allele", "effect.ci"),
       #leftlabs = c("Variant (Bitter\nTaste Allele)", "Taster\nAllele", "Wt\n(%)", "IVW Beta [95% CI]\nper Taster Allele"),
       leftlabs = c("Variant (Bitter\nTaste Allele)", "Taster\nAllele", "IVW Beta [95% CI]\nper Taster Allele"),
       smlab = "2-hr Glucose (mg/dL)\nadjusted for BMI", text.fixed = "Pooled Effect Estimate",
       just.addcols.left = "left", just.studlab = "left", just = "left",
       fontsize = 10,
       squaresize = 0.75, lwd.diamond = 5, lwd.square = 1.25, 
       colgap.left = c("5 mm"))
# fs.axis = 6, fs.heading = 8, fs.smlab = 8,
dev.off()

## Meta-analysis of FG
mr_tas2r38_fg <- cbind.data.frame(
  Study=rep("Manning_2012_MAGIC_FG_BMIadj", 3), 
  SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
  Beta=c(-0.005300000, -0.007300000, -0.006700000 ), SE=c(0.0038, 0.0035, 0.0035), P_value=c(0.1695, 0.03786, 0.05369)) %>%
  mutate(Beta_mgdl=Beta*18.018, SE_mgdl=SE*18.018) %>%
  mutate(Variance=(SE_mgdl^2)) %>%
  mutate(W = c(1/Variance)) %>%
  mutate(IVW=Beta_mgdl*W) %>%
  mutate(weights = 1/(SE^2)) %>% 
  mutate(weights_pct = round(weights/sum(weights)*100, 1)) %>%
  mutate(SNP_alleles=c("rs713598 (C>G)", "rs1726866 (G>A)", "rs10246939 (T>C)")) %>%
  mutate(Bitter_Allele=c("G", "G", "C")) %>%
  mutate(Diplotype=1:3)

m.gen_fg <- metagen(TE = Beta_mgdl, seTE=SE_mgdl, studlab = SNP_alleles, method.random.ci = "classic",
                    common=TRUE, random = FALSE, sm= "SMD", method.tau = "REML", data=mr_tas2r38_fg)

forest(m.gen_fg, layout = "RevMan5",
       sortvar = Diplotype, 
       digits=2, #at=seq(-2,2,1),
       colgap.forest=c("10 mm"),
       col.square =  "#96A0B3", col.square.lines = "#96A0B3", col.inside = "black",
       col.diamond = "#435269", col.diamond.lines = "#435269",
       leftcols = c("SNP", "Bitter_Allele", "effect.ci"),
       leftlabs = c("Variant (Bitter\nTaste Allele)", "Taster\nAllele", "IVW Beta [95% CI]\nper Taster Allele"),
       smlab = "Fasting Glucose (mg/dL)\nadjusted for BMI", text.fixed = "Pooled Effect Estimate",
       just.addcols.left = "left", just.studlab = "left", just = "left",
       fontsize = 10,
       squaresize = 0.75, lwd.diamond = 5, lwd.square = 1.25, 
       colgap.left = c("5 mm"))


######################################
## Summarize SNP-level associations ##
######################################
bitter_snps
print_lm(exposure = "taste_diplos.num", outcome="glu_mgdl", covariates = models.l$BMI, label)
summary(lm("glu_mgdl~"))

tas2r38_glu <- rbind.data.frame(
  
cbind.data.frame(
  Study=rep("Saxena_MAGIC_2hG_AdjBMI", 3), 
  SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
  Beta=c(-0.0790, -0.0510, -0.0570), SE=c(0.022, 0.021, 0.021), P_value=c(0.0002571, 0.01228, 0.006042)) %>%
  mutate(Beta_mgdl=Beta*18.018, SE_mgdl=SE*18.018),

  Study=rep("Manning_2012_MAGIC_FG_BMIadj", 3), 
  SNP=c("rs713598_G", "rs1726866_G", "rs10246939_C"), Effect_Allele=c("G", "G", "C"),
  Beta=c(-0.005300000, -0.007300000, -0.006700000 ), SE=c(0.0038, 0.0035, 0.0035), P_value=c(0.1695, 0.03786, 0.05369)) %>%

  
  
###################################
###  Run Sex-Stratified Models  ###
###################################

sex.l <- list(c("Female", "Male"))

## NO significanat sex interactions
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*sex+", models.nofast.l$Base)), data=analysis))
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*sex+", models.nofast.l$Diet)), data=analysis))
summary(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num*sex+", models.sensitivity.nofast.l$SES)), data=analysis))

## Run sex-stratified models, 
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(sex.l, function(x) {
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {
      models.l.use<-gsub("sex","",models.sensitivity.nofast.l$Diet.Patterns) ; dat.use<-analysis %>% filter(fast_cat=="0to2hr")} else {
      models.l.use <- gsub("sex","",models.sensitivity.l$Diet.Patterns) ; dat.use <- analysis} ; return(
        do.call(rbind.data.frame, lapply(1:length(models.l.use), function(v) { 
          print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
                   label.outcome = "temp" ,##paste0(x, ".", names(glycemic_vars_mgdl.l)[g]),
                   covariates = models.l.use, label="Diet.Patterns", data=dat.use %>% 
                     filter(sex == x), lm_trend = F) })) ) })) })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_tastediplos_sexStrat_mDiet_mgdl_",tag,".csv"), row.names = F)


confint(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", gsub("sex", "", models.sensitivity.nofast.l$Diet.Patterns))), data=analysis %>% filter(sex=="Female")))
confint(lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", gsub("sex", "", models.sensitivity.nofast.l$Diet.Patterns))), data=analysis %>% filter(sex=="Male")))


#####################################
###  Mediation Analysis: Attempt  ###
#####################################

library(mediation)
mediator <- "dietPC3"
mediate.mod <- lm(formula(paste0(mediator, "~taste_diplos.num+", models.nofast.l$Lifestyle)), data=analysis %>% filter(complete.cases(glu2hr_mgdl)))
summary(mediate.mod)
full.mod <- lm(formula(paste0("glu2hr_mgdl~taste_diplos.num+", models.nofast.l$Lifestyle, "+", mediator)), data=analysis)
summary(full.mod)
results <- mediate(mediate.mod, full.mod, treat = 'taste_diplos.num', mediator = mediator, boot=T, sims = 100)
summary(results)


################################
## Mendelian Randomization ?? ##
################################

#library(TwoSampleMR)
#ACCESS_TOKEN="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqdWxpZS5nZXJ2aXNAdHVmdHMuZWR1IiwiaWF0IjoxNzM1NjA2MzAzLCJleHAiOjE3MzY4MTU5MDN9.WytdJyfXGJug1-iod65Y8pmN1fUBIruW78t4zPxZTL0H8KdtsyuABB0KzZ-fbzIQZ5zx58yNx_5Xlf2Zo6btho53MNsI_VlZ3pOE0qQ3_cNMiJBDgHeGvYEQs2IVv5Mk0VEL5aSav4gYNoMAO1aeH_hHGX4rK57T2tN_ur7kBFLmsYNX--g1-EoGGi7hq1wiedR9Wm0YEwBCBotDudbCtOPlVbrsz4D8-LP-nT8Rq8Acz7kYm2oT84B_ldXWwkS5xcpzaMewrOEtmPc4rggbpROn2TN2YdSv7iu5ESz82CIqgL5vttWc9_2kPPIH8G9scbpDN-mlLW7hoFVMjtU_6A"
#twohrgluc <- extract_instruments(outcomes='ebi-a-GCST000569')

library(MendelianRandomization)
mr.glu <- as.data.frame(
  mr_tas2r38_rg %>% mutate(snps = gsub("_.*", "", SNP), outcome="glu2hr") %>%
    rename(beta.outcome=Beta, se.outcome=SE, p.outcome=P_value) %>%
    dplyr::select("snps", "outcome", "beta.outcome", "se.outcome", "p.outcome") %>%
    left_join(
      as.data.frame(cbind(snps=c("rs1726866","rs10246939"), exposure="Bitter taste perception", 
                          beta.exposure=c(0.534, 0.968), se.exposure=c(0.0332, 0.028), 
                          p.exposure=c(3e-59, 3e-199))) %>%
        mutate(across(c("beta.exposure", "se.exposure", "p.exposure"), ~as.numeric(.))), 
      by="snps"
    )
) %>% filter(complete.cases(.))


obj <- MendelianRandomization::mr_input(bx=mr.glu$beta.exposure, bxse = mr.glu$se.exposure,
                                 by=mr.glu$beta.outcome, byse = mr.glu$se.outcome,
                                 exposure="exposure",  outcome="outcome",
                                 snps=mr.glu$snps)

res <- MendelianRandomization::mr_ivw(obj)



#######################################################################
## REVISED sensitivity analyses: adjusting for individual food items ##
#######################################################################

## Associatino of diplotype summarized over all subsequent fasting hours
print_lm(exposure="taste_diplos", outcome="glu_mgdl", covariates = models.l$BMI, label="Glucose > 2hr",
         data=analysis %>% filter(fasting_hrs>2))

#   outcome         model exposure     n        beta         se                 p                 f               f_p           trend_p
#1 glu_mgdl Glucose > 2hr  AVI/AVI 54035          NA         NA              <NA> 0.961017511132983 0.382505680552789 0.629298841886853
#2 glu_mgdl Glucose > 2hr  AVI/PAV 78481  0.02452375 0.05071102 0.628671803632185              <NA>              <NA>                  
#3 glu_mgdl Glucose > 2hr  PAV/PAV 28520 -0.04637137 0.06640132 0.484959248814611              <NA>              <NA>                  
  
print_lm(exposure="taste_diplos.num", outcome="glu_mgdl", covariates = models.l$BMI, label="Glucose > 2hr",
         data=analysis %>% filter(fasting_hrs>2))  
#   outcome         model         exposure      n        beta         se         p
#1 glu_mgdl Glucose > 2hr taste_diplos.num 161036 -0.01562937 0.03237805 0.6292988

get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = "glu_mgdl", 
            covars = models.l$BMI, data=analysis %>% filter(fasting_hrs>2)) 
#$emm
#. outcome model      exposure   level     n   emmean         SE     df    lowCI     upCI     anv.p
#1 glu_mgdl  <NA> taste_diplos AVI/AVI 54035 88.04504 0.05282505 160997 87.94150 88.14858 0.3825057
#2 glu_mgdl  <NA> taste_diplos AVI/PAV 78481 88.06956 0.04803828 160997 87.97541 88.16372 0.3825057
#3 glu_mgdl  <NA> taste_diplos PAV/PAV 28520 87.99867 0.06427470 160997 87.87269 88.12465 0.3825057

#$anv
#Analysis of Variance Table

#Response: glu_mgdl
#Df        Sum Sq  Mean Sq F value    Pr(>F)    
#exp            2      158      79    0.9610 0.3825057    

## BItter Genotypes and food groups --------------------------------------------

diet_vars.labs <- c(
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

diet_vars <- c("raw_veg_QT", "cooked_veg_QT", "fresh_fruit_QT", "dried_fruit_QT", "oily_fish_QT", "nonoily_fish_QT", 
              "procmeat_QT", "poultry_QT", "cheese_QT", "beef_QT", "lamb_QT",  "pork_QT", "bread_type_white_vs_brown_or_whole_BIN",
              "bread_intake_QT", "milk_type_full_vs_low_or_nonfat_BIN", "cereal_type_sugar_vs_any_bran_BIN", "cereal_intake_QT", 
              "spread_type_butter_vs_any_other_BIN", "coffee_type_decaf_vs_regular_BIN",  "coffee_QT", "tea_QT", "water_QT",
              "addsalt_always_often_vs_nrs_BIN",  "hotdrink_temp_hot_or_vhot_vs_warm_BIN")

## Add binary lifestyle variables (??)
analysis <- analysis %>%
  mutate(smoker_current_vs_never_BIN = ifelse(smoke_level.lab=="Current", 1, ifelse(smoke_level.lab=="Never",0,NA)),
         smoker_former_vs_never_BIN = ifelse(smoke_level.lab=="Previous", 1, ifelse(smoke_level.lab=="Never",0,NA)),
         physact_low_vs_high_BIN = ifelse(physact_level.lab=="Low", 1, ifelse(physact_level.lab=="High",0,NA)),
         physact_mod_vs_high_BIN = ifelse(physact_level.lab=="Moderate", 1, ifelse(physact_level.lab=="High",0,NA)),
         alch_daily_rarelynever_BIN = ifelse(alch_freq.lab=="Daily or almost daily", 1, 
                                             ifelse(alch_freq.lab %in% c("Never", "Special occasions only"),0, NA)),
         alch_weekly_rarelynever_BIN = ifelse(alch_freq.lab %in% c("3-4 per week", "1-2 per week"), 1,
                                              ifelse(alch_freq.lab %in% c("Never", "Special occasions only"),0, NA)),
         alch_monthly_rarelynever_BIN = ifelse(alch_freq.lab=="1-3 per month", 1,
                                               ifelse(alch_freq.lab %in% c("Never", "Special occasions only"),0, NA))
         )

## tas2r38 and food and lifestyle traits 
diplo_diet <- do.call(rbind.data.frame, lapply(1:length(diet_vars.labs), function(d) {
  print_lm(exposure="taste_diplos.num", outcome=names(diet_vars.labs)[d], covariates = models.nofast.l$Base, label=diet_vars.labs[d])
})) ; diplo_diet %>% fwrite("../data/processed/sens_tab_lm_diplo_dietlife_mBMI_v6.csv")

## food ad lifestyle traits and 2hr glucose
diplo_2hrglu <- do.call(rbind.data.frame, lapply(1:length(diet_vars.labs), function(d) {
  print_lm(outcome="glu2hr_mgdl", exposure=names(diet_vars.labs)[d], covariates = models.nofast.l$BMI, label=diet_vars.labs[d])
})) ;  diplo_2hrglu %>% fwrite("../data/processed/sens_tab_lm_dietlife_glu2hr_mBMI_v6.csv")





# bitter SNPs-diet
snps_diet <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(diet_vars.labs), function(d) {
    print_lm(exposure=snp, outcome=names(diet_vars.labs)[d], covariates = models.nofast.l$BMI, label=diet_vars.labs[d]) }))
}))

snps_diet %>% fwrite("../data/processed/sens_tab_lm_bittersnps_dietlife_mBMI_v6.csv")

# bitter SNPs-diet
snps_diet <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(diet_vars.labs), function(d) {
    print_lm(exposure=snp, outcome=names(diet_vars.labs)[d], covariates = models.nofast.l$BMI, label=diet_vars.labs[d]) }))
}))

snps_diet %>% fwrite("../data/processed/sens_tab_lm_bittersnps_diet_mBMI_v6.csv")




## Lifestyle/Diet Covariates & Glycemic traits --------------------------------

## Summarize lm associations of lifestyle/ses vars with glycemic outcomes
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(lifestyle_ses), function(v) { 
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {cov<-models.nofast.l$BMI ; dat<-analysis %>% filter(fast_cat=="0to2hr")} else {
      cov <- models.l$BMI ; dat <- analysis} ; return(
        print_lm(exposure = names(lifestyle_ses)[v], outcome = glycemic_vars_mgdl.l[g], label.outcome=names(glycemic_vars_mgdl.l)[g],
                 covariates = cov, label=paste0("BMI_",lifestyle_ses[v]), data=dat, digits=c(1,1,6)) ) }))  })) #%>% 
#  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_lifestyle_ses_mBMI_mgdl_",tag,".csv"), row.names = F)


## diet traits & 2hr glucose
diet_vars <- c(dietPC3="Diet PC3", dietPC4 = "Diet PC4", dietPC6="Diet PC6", 
               dietPC11="Diet PC11", dietPC12="Diet PC12", dietPC13 = "Diet PC13",
               dietPC15="Diet PC15", dietPC16="Diet PC16", dietPC17 = "Diet PC17")


## Summarize lm associations of lifestyle/ses vars with glycemic outcomes
do.call(rbind.data.frame, lapply(1:3, function(g) {
  do.call(rbind.data.frame, lapply(1:length(diet_vars), function(v) { 
    if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") {cov<-models.nofast.l$BMI ; dat<-analysis %>% filter(fast_cat=="0to2hr")} else {
      cov <- models.l$BMI ; dat <- analysis} ; return(
        print_lm(exposure = names(diet_vars)[v], outcome = glycemic_vars_mgdl.l[g], label.outcome=names(glycemic_vars_mgdl.l)[g],
                 covariates = cov, label=paste0("BMI_",diet_vars[v]), data=dat) ) }))  })) %>% 
  write.csv(paste0("../data/processed/sens_tab_lm_allglycemic_covars_dietPCs_mgdl_",tag,".csv"), row.names = F)


## END_OF_SCRIPT




