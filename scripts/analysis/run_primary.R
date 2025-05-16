# Title: Primary analysis
# Date Updated: 08-27-2024


## load functions
source("../scripts/pantry.R")
library(emmeans)

#command args
args=commandArgs(trailingOnly = T)
ANC="EUR"

## TAG DIRECTORIES
#tag="v7" # winsorized all numeric variables and removed exclusion based on >|5SD| for glucose (since they are now winsorized)
tag="v6" # updated confounder variables for PA & additional variables for smoking & alch added; & glu2hr & restrictions applied in postprocesing
#tag="v5" # differentiated by addign adjustment for assessment center (factor)


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", tag,"_", ANC, ".rda"))

dim(analysis)  # N=241378


# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common))>0.005))
common_diplos_AA <- names(which(prop.table(table(analysis$diplos_common_AA))>0.005))
taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")

glycemic_vars.l <- list("Glucose"="glu", "0-2hr Glucose"="glu2hr", "HbA1c"="hba1c")
glycemic_vars_mgdl.l <- list("Glucose (mg/dL)"="glu_mgdl", "0-2hr Glucose (mg/dL)"="glu2hr_mgdl", "HbA1c (%)"="hba1c")

# create output directories
outDir=paste0("../data/processed/")



################################################
## Descriptive analysis of participant cohort ##
################################################

#Diplotype distributions by assessment center
print_summary_table(data=analysis, vars_to_summarize = c("taste_diplos" = "TAS2R38 Diplotypes"),
                    var_strata = "ac.f", p_types = "descriptive", factor_vars = "taste_diplos") %>%
  write.csv(paste0("../data/processed/descr_tab_tastediplo_ac_",tag,".csv"))
 
## Demographic, behavioral & risk factors
vars_to_summarize <- c(
  age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking Status", cigarettes_per_day="Cigarettes per day (among smokers)", 
  physact_level.lab="Physical Activity Level", physact_met_excess="Physical Activity, excess MET/wk",
  alch_freq.lab="Alcohol Frequency", alch_1to4d.lab="Alcohol Frequency, 1-4 per week",
  alch_drinker_status.lab = "Alcohol Drinker", #alch_heavydrinker = "Heavy alcohol drinker",
  alch_drinks_per_week = "Drinks per week",
  educ_isced.lab="Education Level", income_level.lab = "Income Level",
  sbp_adj="SBP, mmHg (med adjusted)", dbp="DBP, mmHg (med adjusted)",
  tg="Triglyceride, mmol/L", ldl="LDL, mmol/L", hdl="HDL, mmol/L", raw_veg_QT="Raw vegetables",
  coffee_QT="Coffee, cups/day", tea_QT="Tea, cups/day", addsalt_3lvl.lab="Added salt",
  hba1c_max="HbA1C, %" #glu="Glucose, mmol/L",   glu2hr="0-2hr Glucose, mmol/L"
)


# ========================================
## Run over all diplotypes with >0.005 
# ========================================

diplo_sets <- c(diplos_common_AA="diplo_AA_005", taste_diplos="diplo_AA_taste", taster_status="diplo_AA_taste_other")

analysis <- analysis %>% mutate(
  diplos_common_AA = factor(diplos_common_AA, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV", "AVI/AAV", "AAV/PAV")))

## Check for factor coding of relevant variables
lapply(c("smoke_level.lab", "alch_freq.lab", "physact_level.lab", "income_level.lab", "educ_isced.lab"), function(v) {
  str(analysis[[v]]) }) #--> Recode education as factor variable
analysis <- analysis %>% mutate_at("educ_isced.lab", ~factor(., levels=paste0("Level ", 1:5)))
  
for (i in 1:length(diplo_sets)) {
  
  diplo_strata_var = names(diplo_sets)[i]
  diplo_strata_order =  levels(analysis[[diplo_strata_var]])
  
  rbind.data.frame(
    
    # print summary table with digits: 1=means, 0=pct, 4=p-values
    print_summary_table(
      vars_to_summarize = vars_to_summarize, 
      var_strata = diplo_strata_var, var_strata_order = c(diplo_strata_order),
      p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(1,1,4), data=analysis) ,
    
    # Add table for sex-stratified alcohol intake (g/day)
    do.call(rbind.data.frame, lapply(c("Female", "Male"), function(x){
      print_summary_table(vars_to_summarize = c(alch_drinks_per_week=paste0(x, ": Drinks per wweek (drinkers)")),
                          var_strata = diplo_strata_var, var_strata_order = c(diplo_strata_order),
                          p_print = T, p_adjust = c("age", "ac.f"), digits = c(3,0,4), 
                          data=analysis %>% filter(sex==x))
      })),
    
    # Add tables for traits with 0 digits (BP)
    print_summary_table(vars_to_summarize = c(sbp_adj="SBP, mmHg (med adjusted)", dbp_adj="DBP, mmHg (med adjusted)"),
                        var_strata = diplo_strata_var, var_strata_order = c(diplo_strata_order),
                        p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(0,0,4), data=analysis),  
    
    # Add tables for traits with 3 digits (glucose & 2hr glucose)
    print_summary_table(vars_to_summarize = c(glu="Glucose, mmol/L (3 digits)", glu2hr="0-2hr Glucose, mmol/L", alch_drinks_per_week="Drinks per week"),
                        var_strata = diplo_strata_var, var_strata_order = c(diplo_strata_order),
                        p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(3,0,4), data=analysis)
    ) %>% 
    
    write.csv(paste0("../data/processed/descr_tab_",diplo_sets[i], "_", tag,".csv"))
}
  


###############################################
###   Preliminary Analysis of Diet Traits   ###
###############################################

## Descriptive analysis of diet traits across dipllotypes
diet_to_summarise <- c(diet_labels, TCALS="Total Energy, kcal/d", CHO_pct="Carbohydrate, %kcal", 
                       PRO_pct="Protein, %kcal", FAT_pct="Fat, %kcal", MUFA_pct="MUFA, %kcal",
                       PUFA_pct="PUFA, %kcal", SFA_pct="SFA, %kcal")

descr_dietTraits <- print_summary_table(vars_to_summarize = diet_to_summarise, 
  var_strata = "taste_diplos", var_strata_order = taste_diplos,
  p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(3,3,4), data=analysis) 

beta_diplo_diet <- do.call(rbind.data.frame, lapply(1:24, function(d) {
  print_lm(exposure = "taste_diplos.num", outcome=names(diet_labels)[d], covariates = "taste_diplos.num", label=diet_labels[d])
}))

## Pearson correlations of PAV haplotypes (.num) with diet traits
cor.test(analysis$taste_diplos.num, analysis$raw_veg_QT)
summary(lm(raw_veg_QT~taste_diplos.num+age+sex+ac.f, data=analysis))

## Descriptive analysis of diet PCs across diplotypes
dietPCs_labs <- paste0("Diet PC", 1:24) ; names(dietPCs_labs) <- dietPCs
descr_dietPCs <- print_summary_table(
  vars_to_summarize = dietPCs_labs, var_strata = "taste_diplos", 
  var_strata_order = taste_diplos, p_print = T, 
  p_adjust = c("age", "sex", "ac.f"), digits = c(3,3,4), data=analysis) 


# ========================================================================
## Calculate EMMs for Dietary Traits by diplotype, adjusting for age+sex
# ========================================================================
diet_traits <- names(diet_labels)

emm_diplos_diet.l <- lapply(diet_traits, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = "age+sex+ac.f", data=analysis) 
}) ; names(emm_diplos_diet.l) = diet_traits
do.call(rbind.data.frame, lapply(diet_traits, function(d) {
  emm_diplos_diet.l[[d]]$emm} )) %>% mutate(DietLabel=rep(diet_labels, each=3)) %>%
  write.csv(paste0(outDir, "descr_tab_emm_dietTraits_diplos_sexageac_",tag,".csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  Diet=diet_traits,
  F_stat=sapply(1:24, function(i) {emm_diplos_diet.l[[i]]$anv$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_diet.l[[i]]$anv$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(diet_traits, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex+ac.f")), data=analysis)))[2,]
    })) %>% as.data.frame() %>% mutate(Diet=diet_traits, DietLabel=diet_labels) %>% 
      rename(Beta=Estimate, SE='Std. Error', "T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="Diet") %>%
  write.csv(paste0(outDir, "descr_tab_teststat_dietTraits_diplos_sexageac_",tag,".csv"), row.names = F)


cor.test()
cor.test()

# ==============================================================================
## Calculate EMMs for Nutrient intakes by diplotype, adjusting for age+sex+ac.f
# ==============================================================================
nutrients <- c("TCALS", "CHO_pct", "PRO_pct", "FAT_pct", "MUFA_pct", "PUFA_pct", "SFA_pct")

emm_diplos_nutr.l <- lapply(nutrients, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = "age+sex+ac.f", data=analysis) 
}) ; names(emm_diplos_nutr.l) = nutrients
do.call(rbind.data.frame, lapply(nutrients, function(d) {
  emm_diplos_nutr.l[[d]]$emm} )) 

as.data.frame(cbind(
  Diet=nutrients,
  F_stat=sapply(1:7, function(i) {emm_diplos_nutr.l[[i]]$anv$`F value`[1]}),
  P_val_f=sapply(1:7, function(i) {emm_diplos_nutr.l[[i]]$anv$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(nutrients, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex+ac.f")), data=analysis)))[2,3:4]
    })) %>% as.data.frame() %>% mutate(Diet=nutrients) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="Diet")


# ============================================================
## Calculate EMMs for dPCs by diplotype, adjusting for age+sex
# ============================================================
dietPCs <- paste0("dietPC", 1:24)

emm_diplos_dPCs.l <- lapply(dietPCs, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = "age+sex+ac.f", data=analysis) 
}) ; names(emm_diplos_dPCs.l) = dietPCs
do.call(rbind.data.frame, lapply(dietPCs, function(dPC) {
  emm_diplos_dPCs.l[[dPC]]$emm} )) %>%
  write.csv(paste0(outDir, "descr_tab_emm_dietPC_diplos_sexageac_",tag,"_temp.csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  DietPC=1:24,
  F_stat=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(dietPCs, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex+ac.f")), data=analysis)))[2,3:4]
      })) %>% as.data.frame() %>% mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="DietPC") %>%
  write.csv(paste0(outDir, "descr_tab_teststat_dietPC_diplos_sexageac_",tag,"_temp.csv"), row.names = F)



##################################################
###  Preliminary analysis of glycemic traits.  ###
################################################## 

## Add glycemic traits in mg/dL
analysis <- analysis %>% mutate(glu_mgdl = glu*18.018,
                                glu2hr_mgdl = glu2hr*18.018)

glycemic_vars <- c(glu="Glucose", glu2hr="0-2hr Glucose", hba1c="HbA1c") 
glycemic_vars_mgdL <- c(glu_mgdl="Glucose (md/dL)", glu2hr_mgdl="0-2hr Glucose (mg/dL)", hba1c="HbA1c (%)") 

# Describe glucose by diplotype & continuous fasting time: categorical
print_summary_table(data=analysis, vars_to_summarize = glycemic_vars_mgdL, var_strata = "taste_diplos", 
                    var_strata_order = c("AVI/AVI", "AVI/PAV", "PAV/PAV"), digits = c(3,0,6)) %>%
  mutate(P_adjust="Sex+Age", .before="Total") %>%
  write.csv(paste0("../data/processed/descr_tab_allglycemic_tastediplos_noadj_mgdl_",tag,".csv"))

do.call(rbind, lapply(1:3, function(g) {
  print_lm(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], 
           label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = "age+sex", 
           label = "Categorical.Sex+Age", lm_trend = T, data=analysis) }))  %>%
  write.csv(paste0("../data/processed/descr_tab_lm_allglycemic_tastediplos_sexage_mgdl_",tag,".csv"))

do.call(rbind, lapply(1:3, function(g) {
  print_lm(exposure = "taste_diplos.num", outcome = glycemic_vars_mgdl.l[[g]], 
                          label.outcome = names(glycemic_vars_mgdl.l)[g], covariates = "age+sex", 
                          label = "Continuous.Sex+Age", data=analysis, 
                          lm_trend = T) })) %>%
  write.csv(paste0("../data/processed/descr_tab_lm_allglycemic_tastediplos_add_sexage_mgdl_",tag,".csv"))



######################################################
###   Primary Analysis with Canonical Diplotypes   ###
######################################################

## Outcomes
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
glu2hr_var.l <- list(variables = list(glu2h = "0-2hr Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))

glucose_var_mgdl.l <- list(variables = list(glu_mgdl = "Glucose, mg/dL"))
glu2hr_var_mgdl.l <- list(variables = list(glu2h_mgdl = "0-2hr Glucose, mg/dL"))

## Model covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+ac.f+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level.lab+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC3+dietPC4+dietPC6+dietPC11+dietPC12+dietPC13+dietPC15+dietPC17")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)
models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)


# ===============================
## Canonical diplotypes & RG 
# ===============================

## Categorical diplotypes -----
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome="Glucose (md/dL)",
           label=names(models.l)[m], covariates = models.l[[m]], lm_trend = T, 
           data = analysis) } )) %>%
  write.csv(paste0(outDir, "res_tab_lm_rg_tastediplo_mgdl_", tag, ".csv"), row.names = F)

## Continuous diplotypes -----
do.call(rbind.data.frame, lapply(1:length(models.l), function(m){
  print_lm(outcome = "glu_mgdl", exposure = "taste_diplos.num", label.outcome="Glucose (md/dL)",
           covariates = models.l[[m]], label=names(models.l)[m])
})) %>% write.csv(paste0(outDir, "res_tab_lm_rg_tastediplo_add_mgdl_", tag, ".csv"), row.names = F)


# ===============================
## Canonical diplotypes & HbA1c
# ===============================

# categorical diplotypes ----
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label=names(models.l)[m], 
           label.outcome = "HbA1c", covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>% 
  write.csv(paste0(outDir, "res_tab_lm_a1c_tastediplo_", tag, ".csv"), row.names = F)

## Continuous diplotypes -----
do.call(rbind.data.frame, lapply(1:length(models.l), function(m){
  print_lm(outcome = "hba1c", exposure = "taste_diplos.num", label.outcome="Glucose",
           covariates = models.l[[m]], label=names(models.l)[m])
})) %>% write.csv(paste0(outDir, "res_tab_lm_a1c_tastediplo_add_", tag, ".csv"), row.names = F)


# =====================================
## Canonical diplotypes & 2hr Glucose
# =====================================

# Categorical diplotypes -----
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  print_lm(exposure = "taste_diplos_pav", outcome = "glu2hr_mgdl", label.outcome="0-2hr Glucose (md/dL)",
           label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], 
           lm_trend = T, data = analysis %>% filter(fast_cat == "0to2hr")) })) %>% 
  write.csv(paste0(outDir, "res_tab_lm_2hg_tastediplo_mgdl_",tag,".csv"), row.names = F)
  
# Continuous diplotypes -----
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m){
  print_lm(outcome = "glu2hr_mgdl", exposure = "taste_diplos.num", label.outcome="0-2hr Glucose (md/dL)", 
           covariates = models.nofast.l[[m]], label=names(models.nofast.l)[m]) })) %>% 
  write.csv(paste0(outDir, "res_tab_lm_2hg_tastediplo_add_mgdl_", tag, ".csv"), row.names = F)

## Glucose over all fasting windows ----------
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6hr"))
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
    print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = paste0(fast, "-Glucose (mg/dL)"), 
             label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
             data = analysis %>% filter(fast_cat == fast)) })) })) %>% 
    write.csv(paste0(outDir, "res_tab_lm_gluXfast_tastediplo_mgdl_",tag,".csv"), row.names = F)


do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
    print_lm(exposure = "taste_diplos", outcome = "glu_mgdl", label.outcome = paste0(fast, "-Glucose (mg/dL)"), 
             label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], lm_trend = T, 
             data = analysis %>% filter(fast_cat_with_6to12 == fast)) })) })) %>% 
  write.csv(paste0(outDir, "res_tab_lm_gluXfast12to24_tastediplo_mgdl_",tag,".csv"), row.names = F)


## Testing an interaction term in the diet model with Taste x Fasting category
lapply(models.l, function(m) {
  anova(lm(formula(paste0("glu_mgdl~taste_diplos.num*fast_cat+", m)), data=analysis))
})

## ==============================================================
##  Estimated Marginal Means 
## ==============================================================

## All glycemic outcomes ----------
do.call(rbind.data.frame, lapply(1:3, function(g) {
  
  if(glycemic_vars_mgdl.l[[g]]=="glu2hr_mgdl") { 
    models.l.use <- models.nofast.l ; dat.use <- analysis %>% filter(fast_cat=="0to2hr") } else {
      models.l.use <- models.l ; dat.use <- analysis } 
  
  return( do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
    get_emm.fun(exposure = "taste_diplos", outcome = glycemic_vars_mgdl.l[[g]], covars = models.l.use[[m]], 
                label=names(models.l.use)[m], label.outcome = names(glycemic_vars_mgdl.l)[g],
                reference = "AVI/AVI",  data=dat.use)$emm })) ) })) %>%
  write.csv(paste0("../data/processed/res_tab_emm_allglycemic_tastediplos_mgdl_", tag, ".csv"), row.names = F)


## Glucose stratified by fasting time ----------
#fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))
models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)

do.call(rbind.data.frame, lapply(fast_cat.l, function(f) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$Diet.Patterns, 
              label.outcome = paste0(f, "-Glucose (mg/dL)"), label="Diet.Patterns", reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat==f))$emm
  } )) %>% mutate(fast=gsub("[-].*","",outcome),.before=exposure) %>%
  write.csv(paste0(outDir, "res_tab_emm_gluXfast_tastediplo_mDiet_mgdl_",tag,".csv"), row.names = F)


## Glucose stratified by fasting time with 6-12 and 12-24 parsed out --------------

fast_cat_with_6to12 <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6to12hr", "12+hr"))
analysis <- analysis %>% 
  mutate(fast_cat_with_6to12 = factor(ifelse(fasting_hrs >= 6 & fasting_hrs <= 12, "6to12hr", 
                                             ifelse(fasting_hrs >12, "12+hr", fast_cat)),
                                      levels=c("0to2hr", "3hr", "4hr", "5hr", "6to12hr", "12+hr")))

# Diet Model
do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(f) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$Diet.Patterns, 
              label.outcome = paste0(f, "-Glucose (mg/dL)"), label="Diet.Patterns", reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat_with_6to12==f))$emm })) %>% 
  mutate(fast=gsub("[-].*","",outcome),.before=exposure) %>%
  write.csv(paste0(outDir, "res_tab_emm_gluXfast_6to12_tastediplo_mDiet_nofast_mgdl_",tag,".csv"), row.names = F)

# Base Model
do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(f) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$Base, 
              label.outcome = paste0(f, "-Glucose (mg/dL)"), label="Base", reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat_with_6to12==f))$emm })) %>% 
  mutate(fast=gsub("[-].*","",outcome),.before=exposure) %>%
  write.csv(paste0(outDir, "res_tab_emm_gluXfast_6to12_tastediplo_mBase_nofast_mgdl_",tag,".csv"), row.names = F)

# BMI Model
do.call(rbind.data.frame, lapply(fast_cat_with_6to12, function(f) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu_mgdl", covars = models.nofast.l$BMI, 
              label.outcome = paste0(f, "-Glucose (mg/dL)"), label="BMI", reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat_with_6to12==f))$emm })) %>% 
  mutate(fast=gsub("[-].*","",outcome),.before=exposure) %>%
  write.csv(paste0(outDir, "res_tab_emm_gluXfast_6to12_tastediplo_mBMI_nofast_mgdl_",tag,".csv"), row.names = F)


## Compare vs. emmeans with a taste*fasting interaction ----------------
#emmint <- emmeans(lm(formula(paste0("glu ~ taste_diplos*fast_cat+", models.l$Diet.Patterns)), data=analysis),
#                     ~ taste_diplos*fast_cat, rg.limit=35000)
#contrint <- contrast(emmint, "revpairwise", by = "fast_cat", adjust = "bonferroni")
#emmdat <- emmip(emmint, taste_diplos ~ fast_cat, CIs = TRUE, plotit = FALSE) 
#emmdat %>% write.csv(paste0(outDir, "res_tab_emmint_gluXfast_tastediplo_mDiet_",tag,".csv"), row.names = F)
#ggplot(data = emmdat, aes(x = fast_cat, y = yvar-4.82, color = taste_diplos)) +
#  geom_point(position = position_dodge(0.55), size=3) + 
#  geom_errorbar(aes(ymin=LCL-4.82, ymax=UCL-4.82), position = position_dodge(0.55), width=0.35, linewidth=0.65) + 
#  scale_color_manual(values=palettes$nature_diplos_dominant, name="TAS2R38 Diplotype")
  

########################
## Post-hoc analysis  ##
########################

pct_beta_change.fun <- function(b_base, b_adj) {
  b_pct <- ((b_adj - b_base)/b_base)*100
  return(b_pct)
}

## Quantify the % confounding for BMI, lifestyle & food choices on models for glu
m_base <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[1]])), data=analysis))
m_bmi <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[2]])), data=analysis))
m_life <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[3]])), data=analysis))
m_food <- summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[4]])), data=analysis))

# BMI
pct_beta_change.fun(m_base$coef[2,1], m_bmi$coef[2,1]) #1.063%
#Lifestyle
pct_beta_change.fun(m_bmi$coef[2,1], m_life$coef[2,1]) #8.6%
#Food choices
pct_beta_change.fun(m_life$coef[2,1], m_food$coef[2,1]) #3.9


##########################
## TAS2R38 & T2D status ##
##########################

postprocessed <- readRDS("../data/processed/ukb_postprocessed_EUR_20241227.rda")

t2d <- postprocessed %>% 
  filter(n_geno==1 & n_covars)

## Deesctiptive table
t2d %>% group_by(taste_diplos) %>% summarise(T2D=n_pct(t2d_case))

# Associations with T2D (categorical diplotype) -----
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
  print_glm(exposure = "taste_diplos", outcome = "t2d_case", 
           label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], 
           exp = T, data = t2d) })) %>% 
  write.csv(paste0(outDir, "res_tab_glm_t2d_tastediplo_",tag,".csv"), row.names = F)

# Associations with T2D (continuous diplotype) -----
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
  print_glm(exposure = "taste_diplos.num", outcome = "t2d_case", 
            label=names(models.nofast.l)[m], covariates = models.nofast.l[[m]], 
            exp = T, data = t2d) })) %>% 
  write.csv(paste0(outDir, "res_tab_glm_t2d_tastediplo_num_",tag,".csv"), row.names = F)

summary(glm(formula(paste0("t2d_case~taste_diplos+", models.nofast.l$Diet.Patterns)), family=binomial("logit"), data=t2d))


##EOF


