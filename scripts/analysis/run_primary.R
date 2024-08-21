# load functions
library(tidyverse) ; library(data.table)



#command args
args=commandArgs(trailingOnly = T)
ANC=args[1]

source("../scripts/pantry.R")

## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(N_contrl_compl==1) %>%
  filter(is.na(haplo_0) == F & is.na(haplo_1) == F) %>%
  filter(fasting_hrs <= 24) #%>%
  #filter(glu <= 11.1)

# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common))>0.005))

dim(analysis)

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
  write.csv(paste0(outDir_pf, "emm_dietPC_diplos_sexage.csv"), row.names = F)


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
  write.csv(paste0(outDir_pf, "pstat_dietPC_diplos_sexage.csv"), row.names = F)


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


# ============================
## Taster Status & RG 
# ============================

#analysis$diplos_AA <- relevel(analysis$diplos_AA, ref = "AVI/AVI")

## Run main effects of taster status on RG for canonical diplotypes 
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], 
           data = analysis)  
  })) %>% mutate(model=rownames(.), .before=beta) %>%
  #Add P-values for linear trend
  left_join(
    do.call(rbind, lapply(models.l, function(m){
      coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), data=analysis)))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model=paste0(names(models.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`),
    by="model") %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_v2.csv"), row.names = F)


# ============================
## Taster Status & A1c
# ============================

## Run main effects of taster status on A1C
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = "taste_diplos", outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i],
           data = analysis)  
  })) %>% mutate(model=rownames(.), .before=beta) %>% 
  #Add P-values for linear trend
  left_join(
    do.call(rbind, lapply(models.l, function(m){
      coef(summary(lm(formula(paste0("hba1c~taste_diplos.num+", m)), data=analysis)))[2,3:4]
    })) %>% as.data.frame() %>% mutate(model=paste0(names(models.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`),
    by="model") %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1c_v2.csv"), row.names = F)



## =================================================
##   Taster Status & Glu/A1c by Fasting Time   
## =================================================

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(fast_cat)


# covariates without fasting
models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l))
names(models.nofast.l) <- names(models.l)


## Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                             filter(fast_cat == fast))) %>%
      mutate(model = rep(names(models.nofast.l)[[i]], 3), .before=beta)
  })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
  #Add P-values for linear trend
  left_join(
    do.call(rbind, lapply(fast_cat.l, function(fast){
      do.call(rbind, lapply(models.nofast.l, function(m){
        coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), data=analysis %>% filter(fast_cat == fast))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast.l), "_AVI/AVI")) %>% 
        rename("t" = `t value`, "t_p" = `Pr(>|t|)`) })),
    by="model_fast") %>% 
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_v2.csv"), row.names = F)


## Run main effects of taster status on A1c by fasting time
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = "taste_diplos", outcome = "hba1c_max", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                             filter(fast_cat == fast)))  %>%
      mutate(model = rep(names(models.nofast.l)[[i]], 3), .before=beta)
    })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
  #Add P-values for linear trend
  left_join(
    do.call(rbind, lapply(fast_cat.l, function(fast){
      do.call(rbind, lapply(models.nofast.l, function(m){
        coef(summary(lm(formula(paste0("hba1c~taste_diplos.num+", m)), data=analysis %>% filter(fast_cat == fast))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`) 
      })),
    by="model_fast") %>% write.csv(paste0(outDir_pf, "lm_diplo_x_a1cXfast_v2.csv"), row.names = F) %>%
  #save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1cXfast_v2.csv"), row.names = F)



# ================================================================
## Estimated marginal means glu by fasting time in the diet model
# ================================================================

do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.nofast.l[[4]], reference = "AVI/AVI", 
              data=analysis %>% filter(fast_cat == fast))$emm })) %>%
  mutate(fast = rep(fast_cat, each=3), .before=emmean) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_gluXfast_diet_v2.csv"), row.names = F)




################################
###   Sensitivity Analysis   ###
################################

# ===================================================================
## Subsetting to individuals with plausible total energy intakes
## Plausible: >= 600 & 
##    -for M: <20000 kJ (4780.2 kcal) 
##    -for F: <18000 kj (4302.1 kcal)
# ===================================================================

#Add variables for plausible/missing_implausible kcal 
analysis <- analysis %>%
  mutate(plausible_kcal.f = ifelse(is.na(kcal_plaus)==F, "Plausible", "Implausible_Missing")) 

plausible_kcal.l <- c("Plausible", "Implausible_Missing")

# Run main effects of taster status on RG ------------------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], 
             data = analysis %>% filter(plausible_kcal.f == kcal))
    })) %>% mutate(model=rownames(.), .before=beta) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(models.l, function(m){
        coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), 
                        data=analysis %>% filter(plausible_kcal.f == kcal))))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model=paste0(names(models.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`),
      by="model") 
  })) %>% 
  mutate(plaus_kcal = rep(plausible_kcal.l, each=12), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_plauskcal.csv"), row.names = F)


## Run main effects of taster status on A1c ------------------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = "taste_diplos", outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i], 
             data = analysis %>% filter(plausible_kcal.f == kcal))
    })) %>% mutate(model=rownames(.), .before=beta) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(models.l, function(m){
        coef(summary(lm(formula(paste0("hba1c~taste_diplos.num+", m)), 
                        data=analysis %>% filter(plausible_kcal.f == kcal))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model=paste0(names(models.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`), 
      by="model") 
  })) %>%
  mutate(plaus_kcal = rep(plausible_kcal.l, each=12), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1c_plauskcal.csv"), row.names = F)


## Run main effects of taster status on glucose by fasting time ------------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
      as.data.frame(print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.nofast.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                               filter(fast_cat == fast & plausible_kcal.f == kcal))) %>%
        mutate(model = rep(names(models.nofast.l)[[i]], 3), .before=beta)
    })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(fast_cat.l, function(fast){
        do.call(rbind, lapply(models.nofast.l, function(m){
          coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), 
                          data=analysis %>% filter(fast_cat == fast & plausible_kcal.f == kcal))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast.l), "_AVI/AVI")) %>% 
          rename("t" = `t value`, "t_p" = `Pr(>|t|)`) })),
      by="model_fast")
  })) %>% 
  mutate(plaus_kcal = rep(plausible_kcal.l, each=60), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_plauskcal.csv"), row.names = F)


## Based on observed differences by energy reporting plausibility, test for interaction:
## - Diplotype * SES (income, education)
## - Diplotype * sex
## - Single & multiple interaction terms, in each model

int.l <- c(Plausible_kcal="as.factor(incl_kcal)", Sex="sex", Income="income_level.lab", Educ= "educ_level.lab")

# Interactions with diplo*RG ------------------------
do.call(rbind, lapply(1:length(int.l), function(i) {
  do.call(rbind, lapply(1:length(models.l), function(m){ 
    mod <- lm(formula(paste0("glu~taste_diplos+", models.l[[m]], "+taste_diplos*", int.l[[i]])), data=analysis)
    summary(mod)$coef %>% as.data.frame() %>%
      filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
      mutate(model=names(models.l)[m],
             interaction = rep(names(int.l)[i], nrow(.)),
             term = rownames(.), .before="Estimate")
    })) %>% 
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
  })) %>% write.csv("../data/processed/analysis/int_diplo_x_rg.csv")


# Interactions with diplo.num*RG ------------------------
do.call(rbind, lapply(1:length(int.l), function(i) {
  do.call(rbind, lapply(1:length(models.l), function(m){ 
    summary(lm(formula(paste0("glu~taste_diplos.num+", models.l[[m]], "+taste_diplos.num*", int.l[[i]])), data=analysis))$coef %>% as.data.frame() %>%
      filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
      mutate(model=names(models.l)[m],
             interaction = rep(names(int.l)[i], nrow(.)),
             term = rownames(.),
             .before="Estimate")
  })) %>% 
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_diplo_num_x_rg.csv")



## Interactiond with diplo*glucose by fasting time ------------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind, lapply(1:length(int.l), function(i) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
      summary(lm(formula(paste0("glu~taste_diplos+", models.nofast.l[[m]], "+taste_diplos*", int.l[[i]])), 
                 data=analysis %>% filter(fast_cat == fast)))$coef %>% as.data.frame() %>%
        filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
        mutate(model=names(models.l)[m],
               interaction = rep(names(int.l)[i], nrow(.)),
               term = rownames(.),
               fast=rep(fast, nrow(.)), .before="Estimate")
      }))
    })) %>%
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
  })) %>% write.csv("../data/processed/analysis/int_diplo_x_rg_fasting.csv")
  

## Interactiond with diplo.num*glucose by fasting time ------------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind, lapply(1:length(int.l), function(i) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
      summary(lm(formula(paste0("glu~taste_diplos.num+", models.nofast.l[[m]], "+taste_diplos.num*", int.l[[i]])), 
                 data=analysis %>% filter(fast_cat == fast)))$coef %>% as.data.frame() %>%
        filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
        mutate(model=names(models.l)[m],
               interaction = rep(names(int.l)[i], nrow(.)),
               term = rownames(.),
               fast=rep(fast, nrow(.)), .before="Estimate")
    }))
  })) %>%
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_diplo_num_x_rg_fasting.csv")



## Testing for multiple interactions
do.call(rbind, lapply(2:length(int.l), function(i) {
  do.call(rbind, lapply(1:length(models.l), function(m){ 
    summary(lm(formula(paste0("glu~taste_diplos.num*incl_kcal+", models.l[[m]], "+taste_diplos.num*", int.l[[i]])), data=analysis))$coef %>% as.data.frame() %>%
      filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T | rownames(.) == "incl_kcal") %>% 
      mutate(model=names(models.l)[m],
             interaction = rep(names(int.l)[i], nrow(.)),
             term = rownames(.),
             .before="Estimate")
  })) %>% 
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_kcalplus_diplo_num_x_rg.csv")



##EOF



