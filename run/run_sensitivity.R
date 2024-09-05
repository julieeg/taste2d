# Title: Sensitivity analysis - plausible energy intake 
# Date Updated: 08-27-2024


# create output directories
ANC="EUR"
outDir_pf=paste0("../data/processed/analysis/", ANC,"_")

## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(n_geno==1 & n_t2d_contrl==1 & n_complete==1 & n_fast_24hr==1) %>%
  ## Criteria added: v2 to v3
  filter(find_outliers.fun(glu)==0)  

## load analysis dataframe
dim(analysis) # N=241378

# Restrict to those with 
analysis <- analysis %>% filter(glu_keep == 1)



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

## Model covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)

m5_fib2cho = paste0(m4_diet, "+FIB2CHO")
models_sens.l <- c(models.l, FiberCHO=m5_fib2cho)


plausible_kcal.l <- c("Plausible", "Implausible_Missing")


## TAS2R38 diplotype & RG  ------------------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(1:length(models_sens.l), function(i) {
    print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models_sens.l[[i]], label=names(models_sens.l)[i], 
             data = analysis %>% filter(plausible_kcal.f == kcal))
  })) %>% mutate(model=rownames(.), .before=beta) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(models_sens.l, function(m){
        coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), 
                        data=analysis %>% filter(plausible_kcal.f == kcal))))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model=paste0(names(models_sens.l), "_AVI/AVI")) %>% 
        rename("t" = `t value`, "t_p" = `Pr(>|t|)`), by="model") 
  })) %>% 
  mutate(plaus_kcal = rep(plausible_kcal.l, each=15), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_sens_kcal_v3.csv"), row.names = F)


##TAS2R38 diplotype & A1c ------------------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(1:length(models_sens.l), function(i) {
    print_lm(exposure = "taste_diplos", outcome = "hba1c", covariates = models_sens.l[[i]], label=names(models_sens.l)[i], 
             data = analysis %>% filter(plausible_kcal.f == kcal))
  })) %>% mutate(model=rownames(.), .before=beta) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(models_sens.l, function(m){
        coef(summary(lm(formula(paste0("hba1c~taste_diplos.num+", m)), 
                        data=analysis %>% filter(plausible_kcal.f == kcal))))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model=paste0(names(models_sens.l), "_AVI/AVI")) %>% 
        rename("t" = `t value`, "t_p" = `Pr(>|t|)`), by="model") 
  })) %>% 
  mutate(plaus_kcal = rep(plausible_kcal.l, each=15), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_a1c_sens_kcal_v3.csv"), row.names = F)


## TAS2R38 diplotpe & Glu by fasting time ------------------
models.nofast_sens.l <- c(
  models.nofast.l, Fib2CHO=paste0(models.nofast.l[[4]], "+FIB2CHO"))

bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast_sens.l), function(i) {
      as.data.frame(print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.nofast_sens.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast_sens.l)[i]), data=analysis %>% 
                               filter(fast_cat == fast & plausible_kcal.f == kcal))) %>%
        mutate(model = names(models.nofast_sens.l)[[i]], .before=beta)
    })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
    #Add P-values for linear trend
    left_join(
      do.call(rbind, lapply(fast_cat.l, function(fast){
        do.call(rbind, lapply(models.nofast_sens.l, function(m){
          coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), 
                          data=analysis %>% filter(fast_cat == fast & plausible_kcal.f == kcal))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast_sens.l), "_AVI/AVI")) %>% 
          rename("t" = `t value`, "t_p" = `Pr(>|t|)`) })), by="model_fast")
  })) %>% 
  mutate(plaus_kcal = c(rep(plausible_kcal.l[[1]], 75), rep(plausible_kcal.l[[2]], 73)), .before=n) %>%
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_sens_kcal_v3.csv"), row.names = F)

## STOPPED HERE : SEE IF YOU CAN REPLICATE THE DIFFERENCES WHEN RESTRICTING TO NOT 
## EXTREMES ****AND***** MAKE A TABLE OF EXTREMES BY FASTING TIME AND PLAUSIBILITY ********



################################################
###   Alternative measures of blood glucose  ### 
################################################

# load additional glucose data
glu_add <- fread("../data/processed/ukb_glu_add_20240828.csv")

glucose <- analysis %>% left_join(glu_add, by = "id") %>%
  mutate(
    N_gluB.0_lt5sd = ifelse(abs(zscore.fun(gluB.0)) <= 5,1,0),
    N_gluB.1_lt5sd = ifelse(abs(zscore.fun(gluB.1)) <= 5,1,0),
    N_gluM.0_lt5sd = ifelse(abs(zscore.fun(gluM.0)) <= 5,1,0),
    N_gluM.1_lt5sd = ifelse(abs(zscore.fun(gluM.1)) <= 5,1,0)
  ) ; dim(glucose) #N=241378

glucose <- glucose %>% mutate(
  N_gluB.0 = ifelse(is.na(gluB.0),0,1),
  N_gluB.1 = ifelse(is.na(gluB.1),0,1),
  N_gluM.0 = ifelse(is.na(gluM.0),0,1),
  N_gluM.1 = ifelse(is.na(gluM.1),0,1)
)

glucose_vars <- c(
  gluB.0 = "Glucose_Biochemical_T0",
  gluB.1 = "Glucose_Biochemical_T1",
  gluM.0 = "Glucose_NMR_T0",
  gluM.1 = "Glucose_NMR_T1"
  )


# =====================================================
## Replication of primary analysis in total sample
# =====================================================

## glucose measures & RG ---------------
do.call(rbind.data.frame, lapply(1:length(glucose_vars), function(g) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = "taste_diplos", outcome = names(glucose_vars)[[g]], covariates = models.l[[i]], 
             label=paste0(names(models.l)[i], "_", glucose_vars[[g]]), data = glucose)  })) %>% 
    mutate(model=rownames(.), .before=beta) %>%
    left_join(do.call(rbind, lapply(models.l, function(m){
      coef(summary(lm(formula(paste0(names(glucose_vars)[[g]], "~taste_diplos.num+", m)), data=glucose)))[2,3:4] })) %>% 
        as.data.frame() %>% mutate(model=paste0(names(models.l), "_", glucose_vars[[g]], "_AVI/AVI")) %>% 
        rename("t" = `t value`, "t_p" = `Pr(>|t|)`), by="model") })) %>%
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_sens_glu_v3.csv"), row.names = F)

# glucose measures & Glu x fasting -----------
do.call(rbind.data.frame, lapply(1:length(glucose_vars), function(g) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast_sens.l), function(i) {
      as.data.frame(print_lm(exposure = "taste_diplos", outcome = names(glucose_vars)[[g]], covariates = models.nofast_sens.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast_sens.l)[i], "_", glucose_vars[[g]]), data=glucose %>% 
                               filter(fast_cat == fast))) %>%
        mutate(model = names(models.nofast_sens.l)[[i]], .before=beta)
    })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
    #Add P-values for linear trend
    left_join(do.call(rbind, lapply(fast_cat.l, function(fast){
      do.call(rbind, lapply(models.nofast_sens.l, function(m){
        coef(summary(lm(formula(paste0(names(glucose_vars)[[g]], "~taste_diplos.num+", m)), 
                          data=glucose %>% filter(fast_cat == fast))))[2,3:4]
        })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast_sens.l), "_", glucose_vars[[g]], "_AVI/AVI")) %>% 
          rename("t" = `t value`, "t_p" = `Pr(>|t|)`) })), by="model_fast") })) %>% 
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_sens_glu_v3.csv"), row.names = F)


# ========================================================
## Replication of primary analysis by 24HR availability
# ========================================================

## glucose measures & RG ---------------
bind_rows(lapply(plausible_kcal.l, function(kcal) {
  do.call(rbind.data.frame, lapply(1:length(glucose_vars), function(g) {
    do.call(rbind.data.frame, lapply(1:length(models_sens.l), function(i) {
      print_lm(exposure = "taste_diplos", outcome = names(glucose_vars)[[g]], covariates = models_sens.l[[i]], 
               label=paste0(names(models_sens.l)[i], "_", glucose_vars[[g]]), 
               data = glucose %>% filter(plausible_kcal.f == kcal))  })) %>% 
      mutate(model=rownames(.), .before=beta) %>%
      left_join(do.call(rbind, lapply(models_sens.l, function(m){
        coef(summary(lm(formula(paste0(names(glucose_vars)[[g]], "~taste_diplos.num+", m)), 
                        data=glucose %>% filter(plausible_kcal.f == kcal))))[2,3:4] })) %>% 
          as.data.frame() %>% mutate(model=paste0(names(models_sens.l), "_", glucose_vars[[g]], "_AVI/AVI")) %>% 
          rename("t" = `t value`, "t_p" = `Pr(>|t|)`), by="model") }))
  })) %>% mutate(plaus_kcal = rep(plausible_kcal.l, each = 15)) %>%
  write.csv(paste0(outDir_pf, "lm_diplo_x_rg_sens_glu_plaus_v3.csv"), row.names = F)

# glucose measures & Glu x fasting -----------
do.call(rbind.data.frame, lapply(1:length(glucose_vars), function(g) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast_sens.l), function(i) {
      as.data.frame(print_lm(exposure = "taste_diplos", outcome = names(glucose_vars)[[g]], covariates = models.nofast_sens.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast_sens.l)[i], "_", glucose_vars[[g]]), data=glucose %>% 
                               filter(fast_cat == fast))) %>%
        mutate(model = names(models.nofast_sens.l)[[i]], .before=beta)
    })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
    #Add P-values for linear trend
    left_join(do.call(rbind, lapply(fast_cat.l, function(fast){
      do.call(rbind, lapply(models.nofast_sens.l, function(m){
        coef(summary(lm(formula(paste0(names(glucose_vars)[[g]], "~taste_diplos.num+", m)), 
                        data=glucose %>% filter(fast_cat == fast))))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast_sens.l), glucose_vars[[g]], "_AVI/AVI")) %>% 
        rename("t" = `t value`, "t_p" = `Pr(>|t|)`) })), by="model_fast") })) %>% 
  #Save as .csv
  write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_sens_glu_plaus_v3.csv"), row.names = F)



###############################################################
##  Alternative variants for bitter taste: quinine/caffeine  ## 
###############################################################

# ========================
## Prepare bitter SNPs 
# ========================

# load dosage data for bitter snps
bitter_snps <- read.table("../data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 


## Add OTHER bitter SNPs
bitter <- analysis %>% left_join(bitter_snps, by="id")

dietPCs <- paste0("dietPC", 1:24)

bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")

bitter <- bitter %>% 
  mutate(rs10772420_A = 2-rs10772420_G)


# ===========================================================
## Correlation of bitter SNPs with bitter diet traits & PCs
# ===========================================================

bind_rows(
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    rbind.data.frame(
      print_lm(exposure = snp, outcome="coffee_QT", covariates = "age+sex", label="Coffee", data=bitter),
      print_lm(exposure = snp, outcome="tea_QT", covariates = "age+sex", label="Tea", data=bitter),
      print_lm(exposure = snp, outcome="raw_veg", covariates = "age+sex", label="Raw vegetables", data=bitter),
      print_glm(exposure = snp, outcome="alch_heavydrinker", covariates = "age+sex", label="Heavy drinker", data=bitter),
      print_lm(exposure = snp, outcome="addsalt_freq_QT", covariates = "age+sex", label="Add Salt (cont)", data=bitter) ) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })),
  do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
    do.call(rbind.data.frame, lapply(dietPCs, function(PC) {
      print_lm(exposure = snp, outcome = PC, covariates = "age+sex", label = PC, data=bitter) })) %>%
      mutate(SNP = rep(snp, nrow(.)), .before=Model) })) 
) %>% write.csv("../data/processed/analysis/EUR_lm_snps_bitter_diet_v3.csv", row.names=F)



# ===========================================================
## Primary analysis with RG, A1c & Glu x fasting
# ===========================================================

# outcomes (formatted) 
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))

# covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat+bmi"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5+dietPC6+dietPC7+dietPC8+dietPC9+dietPC10")

models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)


## Bitter SNPs & RG -----------------------
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i], 
             data = bitter)  })) 
  })) %>% write.csv(paste0(outDir_pf, "lm_bittersnps_x_rg_v3.csv"), row.names = F)


## Bitter SNPs & A1c -----------------------
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i], 
             data = bitter)  })) })) %>% 
  write.csv(paste0(outDir_pf, "lm_bittersnps_x_a1c_v3.csv"), row.names = F)


## Bitter SNPs & Glu x fasting time -----------------------

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

models.nofast.l <- gsub("[+]fast_cat", "", models.l)
names(models.nofast.l) <- names(models.l)


## Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
      as.data.frame(print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast.l)[i]), data=bitter %>% 
                               filter(fast_cat == fast))) })) %>% 
      mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) })) %>%
  write.csv(paste0(outDir_pf, "lm_bittersnps_x_gluXfast_v3.csv"), row.names = F)


# ==========================================================================
## Estimated marginal mean: bitter snps & Glu x fasting time (diet model)
# =========================================================================

bitter <- bitter %>% mutate(
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

do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    get_emm.fun(exposure = snp, outcome = "glu", covars = models.nofast.l[[4]], 
                reference = "0_0",
                data=bitter %>% filter(fast_cat == fast))$emm })) %>%
    mutate(fast = rep(fast_cat, each=3), .before=emmean) })) %>%
  write.csv(paste0(outDir_pf, "emm_bittersnps_x_gluXfast_diet_v3.csv"), row.names = F)





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



## First, run without restricting to RG o












## END_OF_SCRIPT


