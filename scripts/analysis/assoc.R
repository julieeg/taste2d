# load functions
library(tidyverse) ; library(data.table)



#command args
args=commandArgs(trailingOnly = T)
ANC=args[1]

source("../scripts/pantry.R")

## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(N_contrl_compl==1) %>%
  filter(is.na(haplo_0) == F & is.na(haplo_1) == F)

# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common))>0.005))

dim(analysis)

# create output directories
outDir_pf=paste0("../data/processed/analysis/", ANC,"_")



############################
###   Primary Analysis   ###
############################

## exposure: common diplotypes (>0.005)
analysis <- analysis %>%
  mutate(diplos_common = ifelse(diplos_common %in% common_diplos, diplos_common, NA))

# outcomes (formatted) 
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))

# covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat+bmi"
m2_lifestyle = paste0(m1_base, "+smoke_level.lab+physact_level+alch_freq.lab")
m3_ses = paste0(m2_lifestyle, "+educ_level.lab+income_level.lab")
m4_diet = paste0(m3_ses, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")

models.l = list("Base"= m1_base, "Lifestyle"=m2_lifestyle, 
                "SES"=m3_ses, "Diet.Patterns"=m4_diet)



# ============================
## Taster Status & RG 
# ============================

analysis$diplos_AA <- relevel(analysis$diplos_AA, ref = "AVI/AVI")

## Run main effects of taster status on RG
tab_lm_rg <- do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "diplos_AA", outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], 
           data = analysis)  
})) %>% mutate(model=rownames(.), .before=beta) 

tab_lm_rg %>% write.csv(paste0(outDir_pf, "lm_diplo_x_rg.csv"), row.names = F)


# ============================
##Taster Status & A1c ##
# ============================

analysis$diplos_AA <- relevel(analysis$diplos_AA, ref = "AVI/AVI")

## Run main effects of taster status on A1C
tab_lm_a1c <- do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "diplos_AA", outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i],
           data = analysis)  
})) %>% mutate(model=rownames(.), .before=beta) 

tab_lm_a1c %>% write.csv(paste0(outDir_pf, "lm_diplo_x_a1c.csv"), row.names = F)



#############################################
## Taster Status & Glucose by Fasting Time ##
#############################################

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

# covariates without fasting
models.nofast.l <- gsub("[+]fast_cat", "", models.l)
names(models.nofast.l) <- names(models.l)

## Run main effects of taster status on glucose by fasting time
tab_lm_0to2g <- do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = "diplos_AA", outcome = "glu", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                             filter(fast_cat == fast)))  %>%
      mutate(model = rep(names(models.nofast.l)[[i]], 5), .before=beta)
  })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta)
})) 

tab_lm_0to2g %>% write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast.csv"), row.names = F)



################################
###   Sensitivity Analysis   ###
################################

# ============================
## Subsetting to individuals with plausible total energy intakes
## Taster Status & RG
# ============================

analysis_plaus <- analysis %>%
  filter(N_contrl_compl_kcal == 1) %>%
  mutate(diplos_AA = factor(diplos_AA, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV", "AAV/AVI", "AAV/PAV")))

## Run main effects of taster status on RG
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "diplos_AA", outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], data = analysis_plaus)
  })) %>% mutate(model=rownames(.), .before=beta) %>%
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_plauskcal.csv"), row.names = F)


## Run main effects of taster status on A1C
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "diplos_AA", outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i], data = analysis_plaus)
  })) %>% mutate(model=rownames(.), .before=beta) %>%
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1c_plauskcal.csv"), row.names = F)


## Run main effects of taster status on glu by fasting time
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = "diplos_AA", outcome = "glu", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis_plaus %>% 
                             filter(fast_cat == fast)))  %>%
      mutate(model = rep(names(models.nofast.l)[[i]], 4), .before=beta)
  })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>%
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_plauskcal.csv"), row.names = F)



# ================================================================
## Estimated marginal mean glu by fasting time in the diet model
# ================================================================

emm_diplo_gluXfast_diet <- do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  get_emm.fun(exposure = "diplos_AA", outcome = "glu", covars = models.nofast.l[[4]], reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat == fast))$emm })) %>%
  mutate(fast = rep(fast_cat, each=5), .before=emmean) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_gluXfast_diet.csv"), row.names = F)

##EOF



