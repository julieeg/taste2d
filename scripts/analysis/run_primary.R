# Title: Primary analysis
# Date Updated: 08-27-2024


## load functions
library(tidyverse) ; library(data.table)


#command args
args=commandArgs(trailingOnly = T)
ANC="EUR"

source("../scripts/pantry/pantry.R")


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(n_geno==1 & n_t2d_contrl==1 & n_complete==1 & n_fast_24hr==1) %>%
  ## Criteria added: v2 to v3
  filter(find_outliers.fun(glu)==0)  

dim(analysis) # N=241378

# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common))>0.005))


# create output directories
outDir_pf=paste0("../data/processed/analysis/", ANC,"_")



###############################################
###   Preliminary Analysis of Diet Traits   ###
###############################################

# ============================================================
## Calculate EMMs for dPCs by diplotype, adjusting for age+sex
# ============================================================

dietPCs <- paste0("dietPC", 1:24)

emm_diplos_dPCs.l <- lapply(dietPCs, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = "age+sex") 
}) ; names(emm_diplos_dPCs.l) = dietPCs
do.call(rbind.data.frame, lapply(dietPCs, function(dPC) {
  emm_diplos_dPCs.l[[dPC]]$emm} )) %>%
  write.csv(paste0(outDir_pf, "emm_dietPC_diplos_sexage_v3.csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  DietPC=1:24,
  F_stat=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv.p$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv.p$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(dietPCs, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex")), data=analysis)))[2,3:4]
      })) %>% as.data.frame() %>% mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="DietPC"
    ) %>%
  write.csv(paste0(outDir_pf, "pstat_dietPC_diplos_sexage_v3.csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype - NO adjustment
 do.call(rbind, lapply(dietPCs, function(d) {
   coef(summary(lm(formula(paste0(d, "~taste_diplos.num")), data=analysis)))[2,3:4]
   })) %>% as.data.frame() %>% 
   mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`) %>%
   write.csv(paste0(outDir_pf, "pstat_dietPC_diplos_raw_v3.csv"), row.names = F)
  

 
######################################################
###   Primary Analysis with Canonical Diplotypes   ###
######################################################

## Outcomes
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))


## Model covariates
#m3_ses = paste0(m2_lifestyle, "+educ_level.lab+income_level.lab")

m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)


# ===============================
## Canonical diplotypes & RG 
# ===============================
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.l)[m], 
           covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>% 
    write.csv(paste0(outDir_pf, "lm_diplo_x_rg_v3.csv"), row.names = F)

## estimated marginal mean glu in each model
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.l[[m]], reference = "AVI/AVI", 
                     data=analysis)
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  return(out) } )) %>%
  mutate(model = rep(names(models.l[[m]]), each=3), .before=emmean) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_glu_v3.csv"), row.names = F)



# ===============================
## Canonical diplotypes & HbA1c
# ===============================
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label=names(models.l)[m], 
           covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>% 
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1c_v3.csv"), row.names = F)


## estimated marginal mean hba1c in each model
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "hba1c", covars = models.l[[m]], reference = "AVI/AVI", 
                     data=analysis)
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  return(out) } )) %>%
  mutate(model = rep(names(models.l[[m]]), each=3), .before=emmean) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_a1c_v3.csv"), row.names = F)



# =====================================
## Canonical diplotypes & 2hr Glucose
# =====================================
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.nofast.l)[m], 
           covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>%
             filter(fast_cat == "0to2hr")) })) %>% 
  write.csv(paste0(outDir_pf, "lm_diplo_x_2hrg_v3.csv"), row.names = F)
  


## ======================================================
##  Canonical diplotypes & Gu over ALL fasting winodows 
## ======================================================

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(fast_cat)

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l))
names(models.nofast.l) <- names(models.l)


## TAS2R38 diplotpes & Glu by fasting time ------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:4, function(m) {
    print_lm(exposure = "taste_diplos", outcome = "glu", label=paste0(fast, ": ", names(models.nofast.l)[m]), 
             covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>%
               filter(fast_cat == fast)) })) })) %>% 
    write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_v3.csv"), row.names = F)


## TAS2R38 diplotype & Glu by fasting time --------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu", 
              covars = models.nofast.l[[4]], reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat == fast))$emm })) %>%
  mutate(fast = rep(fast_cat, each=3), .before=emmean) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_gluXfast_diet_v3.csv"), row.names = F)



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



#######################################################
## Re-Run primary analysis with significant dietPCs  ##
#######################################################

## Model covariates
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")
m_diet_signif = paste0(m3_lifestyle, "+dietPC3+dietPC4+dietPC6+dietPC11+dietPC12+dietPC13+dietPC15+dietPC16+dietPC17")
m_diet_all = paste0(m3_lifestyle, paste0("+dietPC", 1:24, collapse = "+"))
models.l = list("Lifestyle"=m3_lifestyle, "Top10.Diet.Patterns"=m4_diet, "Signig.Diet.Patterns"=m_diet_signif, "All.Diet.Patterns"=m_diet_all)
models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)


## Diplotype & RG ==========
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.l)[m], 
           covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>% 
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_diet_covars.csv"), row.names = F)


## Diplotype & 2hrGlu ==========
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.no.fast.l)[m], 
           covariates = models.no.fast.l[[m]], lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")) } )) #%>% 
#write.csv(paste0(outDir_pf, "lm_diplo_x_rg_v3.csv"), row.names = F)


## estimated marginal mean glu in each model
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.nofast.l[[m]],
                     reference = "AVI/AVI", data=analysis %>% filter(fast_cat == "0to2hr"))
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  return(out) } )) %>%
  mutate(model = rep(names(models.l[[m]]), each=3), .before=emmean) #%>%
  #write.csv(paste0(outDir_pf, "emm_diplo_x_glu_v3.csv"), row.names = F)




##EOF


