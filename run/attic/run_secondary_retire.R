
#############################################################
##  Assoc of Quinine/Caffeine Bitter SNPs with RG & HbA1c  ## 
#############################################################

# ======================
## Set up & load data
# ======================

# load functions & pantry file
library(tidyverse) ; library(data.table)
source("../scripts/pantry/pantry.R")


# load dosage data for bitter snps
bitter_snps <- read.table("../data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(N_contrl_compl==1) %>%
  filter(is.na(haplo_0) == F & is.na(haplo_1) == F) %>%
  filter(fasting_hrs <= 24) 

# create output directories
outDir_pf=paste0("../data/processed/analysis/", ANC,"_")


## Add OTHER bitter SNPs
analysis <- analysis %>% left_join(bitter_snps, by="id")


# ========================
## Prepare bitter SNPs 
# ========================
bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")

analysis <- analysis %>% 
  mutate(rs10772420_A = 2-rs10772420_G)


#####################################################
###   Associations of each Bitter SNP & DietPCs   ###
#####################################################

analysis <- analysis %>% mutate(
  alch_heavydrinker = as.factor(ifelse(alch_freq.lab == "Daily or almost daily" | alch_freq.lab == "3-4 per week", 1, 0)),
  alch_lightdrinker = as.factor(ifelse(alch_freq.lab == "Special occasions only" | alch_freq.lab == "Never", 1, 0))
)

dietPCs <- paste0("dietPC", 1:24)

# Correlations of bitter SNPs with bitter-related diet traits & diet PCs =============
bind_rows(
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    rbind.data.frame(
      print_lm(exposure = snp, outcome="coffee_QT", covariates = "age+sex", label="Coffee"),
      print_lm(exposure = snp, outcome="tea_QT", covariates = "age+sex", label="Tea"),
      print_lm(exposure = snp, outcome="raw_veg", covariates = "age+sex", label="Raw vegetables"),
      print_glm(exposure = snp, outcome="alch_heavydrinker", covariates = "age+sex", label="Heavy drinker"),
      #print_glm(exposure = snp, outcome="alch_lightdrinker", covariates = "age+sex", label="Light drinker"),
      print_lm(exposure = snp, outcome="addsalt_freq_QT", covariates = "age+sex", label="Add Salt (cont)") ) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })),
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    do.call(rbind.data.frame, lapply(dietPCs, function(PC) {
      print_lm(exposure = snp, outcome = PC, covariates = "age+sex", label = PC) })) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })) 
  ) %>% write.csv("../data/processed/analysis/EUR_lm_snps_bitter_diet.csv", row.names=F)

tab<-read.csv("../data/processed/analysis/EUR_lm_snps_bitter_diet.csv")

############################
###   Primary Analysis   ###
############################

# outcomes (formatted) 
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))

# covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat+bmi"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")

models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)


# ============================
## Bitter SNPs & RG 
# ============================


## Run main effects of taster status on RG
tab_lm_bitter_rg <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(exposure = snp, outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], 
           data = analysis)  })) 
}))

tab_lm_bitter_rg %>% write.csv(paste0(outDir_pf, "lm_bittersnps_x_rg.csv"), row.names = F)


# ============================
##Taster Status & A1c ##
# ============================

## Run main effects of taster status on RG
tab_lm_bitter_a1c <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i], 
                  data = analysis)  })) 
}))

tab_lm_bitter_a1c %>% write.csv(paste0(outDir_pf, "lm_bittersnps_x_a1c.csv"), row.names = F)


#######################################################
##   Taster Status & Glucose/HbA1c by Fasting Time   ##
#######################################################

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))


# covariates without fasting
models.nofast.l <- gsub("[+]fast_cat", "", models.l)
names(models.nofast.l) <- names(models.l)


## Run main effects of taster status on glucose by fasting time
tab_lm_bitter_gluXfast <- do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                             filter(fast_cat == fast)))
  })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta)
    })) 
}))

tab_lm_bitter_gluXfast %>% write.csv(paste0(outDir_pf, "lm_bittersnps_x_gluXfast.csv"), row.names = F)



# ================================================================
## Estimated marginal means glu by fasting time in the diet model
# ================================================================

analysis <- analysis %>% mutate(
  rs713598_G.cat = case_when(
    rs713598_G == 0 ~ "00",
    rs713598_G == 1 ~ "01",
    rs713598_G == 2 ~ "11"
  ),
  rs1726866_G.cat = case_when(
    rs1726866_G == 0 ~ "00",
    rs1726866_G == 1 ~ "01",
    rs1726866_G == 2 ~ "11"
  ),
  rs10246939_C.cat = case_when(
    rs10246939_C == 0 ~ "00",
    rs10246939_C == 1 ~ "01",
    rs10246939_C == 2 ~ "11"
  ),
  rs10772420_A.cat = case_when(
    rs10772420_A == 0 ~ "00",
    rs10772420_A == 1 ~ "01",
    rs10772420_A == 2 ~ "11"
  ),
  rs2597979_G.cat = case_when(
    rs2597979_G == 0 ~ "00",
    rs2597979_G == 1 ~ "01",
    rs2597979_G == 2 ~ "11"
  )
)

bitter_snps_cat <- paste0(bitter_snps, ".cat")

do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    get_emm.fun(exposure = snp, outcome = "glu", covars = models.nofast.l[[4]], 
                reference = "00",
                data=analysis %>% filter(fast_cat == fast))$emm })) %>%
    mutate(fast = rep(fast_cat, each=3), .before=emmean) })) %>%
  write.csv(paste0(outDir_pf, "emm_bittersnps_x_gluXfast_diet_v2.csv"), row.names = F)



#######################################################
##  Exploratory: Interactions of TAS2R3 x 24HR data  ## 
#######################################################

# ========================================================
## Explore interactions of TAS2R38 x 24HR on Glu
# ========================================================
# - Diplotype * SES (income, education)
# - Diplotype * sex
# - Single & multiple interaction terms, in each model

int.l <- c(Plausible_kcal="as.factor(incl_kcal)", Sex="sex", Income="income_level.lab", Educ= "educ_level.lab")

## TAS2R38 diplotype x plaus_kcal on RG  ------------------------
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


## TAS2R38 diplotype x plaus_kcal on A1c ------------------------
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


## TAS2R38 diplotype x plaus_kcal on Glu x fasting time  ------------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind, lapply(1:length(int.l), function(i) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
      summary(lm(formula(paste0("glu~taste_diplos+", models.nofast.l[[m]], "+taste_diplos*", int.l[[i]])), 
                 data=analysis %>% filter(fast_cat == fast)))$coef %>% as.data.frame() %>%
        filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
        mutate(model=names(models.l)[m],
               interaction = rep(names(int.l)[i], nrow(.)),
               term = rownames(.),
               fast=rep(fast, nrow(.)), .before="Estimate") }))
  })) %>%
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_diplo_x_rg_fasting.csv")


# # TAS2R38 diplotype x plaus_kcal on A1c x fasting time ------------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind, lapply(1:length(int.l), function(i) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m) {
      summary(lm(formula(paste0("glu~taste_diplos.num+", models.nofast.l[[m]], "+taste_diplos.num*", int.l[[i]])), 
                 data=analysis %>% filter(fast_cat == fast)))$coef %>% as.data.frame() %>%
        filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T) %>% 
        mutate(model=names(models.l)[m],
               interaction = rep(names(int.l)[i], nrow(.)),
               term = rownames(.),
               fast=rep(fast, nrow(.)), .before="Estimate") }))
  })) %>%
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_diplo_num_x_rg_fasting.csv")


## TAS2R38 diplotype x plaus_kcal x [other covariates] on RG --------------
do.call(rbind, lapply(2:length(int.l), function(i) {
  do.call(rbind, lapply(1:length(models.l), function(m){ 
    summary(lm(formula(paste0("glu~taste_diplos.num*incl_kcal+", models.l[[m]], "+taste_diplos.num*", int.l[[i]])), data=analysis))$coef %>% as.data.frame() %>%
      filter(startsWith(rownames(.), "taste") ==T | startsWith(rownames(.), int.l[i]) == T | rownames(.) == "incl_kcal") %>% 
      mutate(model=names(models.l)[m],
             interaction = rep(names(int.l)[i], nrow(.)),
             term = rownames(.), .before="Estimate") })) %>% 
    rename(beta="Estimate", se=`Std. Error`, t_val=`t value`, p=`Pr(>|t|)`)
})) %>% write.csv("../data/processed/analysis/int_kcalplus_diplo_num_x_rg.csv")


## END_OF_SCRIPT









