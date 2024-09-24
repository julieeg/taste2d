# Title: Sensitivity analysis - plausible energy intake 
# Date Updated: 09-05-2024

source("../scripts/pantry/pantry.R")

# create output directories
ANC="EUR"
outDir=paste0("../data/processed/")
tag="v5"

## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(n_geno==1 & n_t2d_contrl==1 & n_complete==1 & n_fast_24hr==1) %>%
  ## Criteria added: v2 to v3
  filter(find_outliers.fun(glu)==0)  %>%
  mutate(n_24hr.lab = ifelse(n_24hr==1, "Has 24HR", "No 24HR")) %>%
  mutate(n_24hr_plaus.lab = factor(ifelse(n_24hr_plaus, "Plausible 24HR", "No 24HR"), levels=c("Plausible 24HR", "No 24HR")))

#Set variables for workflow
n24HR.l <- c("Has 24HR", "No 24HR")
n24HR_plaus.l <- c("Plausible 24HR", "No 24HR")

## load analysis dataframe
dim(analysis) # N=241378

# Load technical covariates 
tech <- readRDS("../data/processed/ukb_glucose_assays.rda") %>% select(id, ends_with(".0") | ends_with(".0.lab"))
analysis <- analysis %>% left_join(tech, by="id")# %>% 
  #mutate(glu_2hg = ifelse(fast_cat=="0to2hr", glu, NA))

## Add Assessment Center ----------------------------------------------------------
ac_labs <- list("Barts"=11012, "Birmingham" = 11021, "Bristol" =	11011, "Bury" =	11008, 
                "Cardiff" =	11003, "Cheadle (revisit)" =	11024, "Croydon" =	11020, 
                "Edinburgh" =	11005, "Glasgow" = 11004, "Hounslow" = 11018, "Leeds" = 11010,
                "Liverpool"=11016, "Manchester"=11001, "Middlesborough"=11017, "Newcastle" =11009, 
                "Nottingham"=11013, "Oxford"=11002, "Reading"=11007, "Sheffield"=11014, "Stockport (pilot)"=10003,
                "Stoke"=11006, "Swansea"=	11022,"Wrexham" =11023, "Cheadle (imaging)"=11025,
                "Reading (imaging)"=11026, "Newcastle (imaging)" =11027, "Bristol (imaging)"=11028)

analysis <- analysis %>% mutate(ac.f = descr_label.fun(., "ac", ac_labs))


# ====================
## Build Models
# ====================
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+ac.f+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC3+dietPC4+dietPC6+dietPC11+dietPC12+dietPC13+dietPC15+dietPC16+dietPC17")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)

# sensitivity models
m_diet_all = paste0(m3_lifestyle, "+", paste0("dietPC",1:24,collapse = "+"))
m_f2c = paste0(m4_diet, "+FIB2CHO")
m_ses = paste0(m_f2c, "+income_level.lab+educ_isced.lab")
m_tech = paste0(m_f2c, "+gluB_aliquot.0")
m_diet_ses = paste0(m4_diet, "+income_level.lab+educ_isced.lab")

models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet"=m4_diet, "Fib2Car"=m_f2c, "SES"=m_ses, "AllDietPCs"=m_diet_all, "Diet+SES"=m_diet_ses)

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)



####################################
###  Adjust for all 24 diet PCs  ###
####################################
rbind.data.frame(
  print_lm(exposure = "taste_diplos", outcome = "glu", label="Glucose ~ All Diet Patterns", 
           covariates = models.l$AllDietPCs, lm_trend = T, data = analysis),
  print_lm(exposure = "taste_diplos", outcome = "glu", label="0-2hr.Glucose ~ All Diet Patterns", 
           covariates = models.nofast.l$AllDietPCs, lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")),
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label="HbA1c ~ All Diet Patterns", 
           covariates = models.l$AllDietPCs, lm_trend = T, data = analysis)) %>% 
  write.csv(paste0("../data/processed/lm_sens_diplos_rg_2hg_m_allDietPC_",tag,".csv"), row.names = T)


## % change in Betas
print_lm(exposure = "taste_diplos", outcome = "glu", label="0-2hr.Glucose ~ All Diet Patterns", 
         covariates=models.nofast.l$Diet, lm_trend = T, 
         data = analysis %>% filter(fast_cat=="0to2hr"))

print_lm(exposure = "taste_diplos", outcome = "glu", label="0-2hr.Glucose ~ All Diet Patterns", 
           covariates=models.nofast.l$AllDietPCs, lm_trend = T, 
         data = analysis %>% filter(fast_cat=="0to2hr"))


########################################
###  Adjust forSocioeconomic Status  ###
########################################

rbind.data.frame(
  print_lm(exposure = "taste_diplos", outcome = "glu", label="Glucose ~ SES", 
           covariates = models.l$`Diet+SES`, lm_trend = T, data = analysis),
  print_lm(exposure = "taste_diplos", outcome = "glu", label="0-2hr.Glucose ~ SES", 
           covariates = models.nofast.l$`Diet+SES`, lm_trend = T, data = analysis %>% filter(fast_cat=="0to2hr")),
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label="HbA1c ~ SES", 
           covariates = models.l$`Diet+SES`, lm_trend = T, data = analysis)) %>% 
  write.csv(paste0("../data/processed/lm_sens_diplos_rg_2hg_m_DietSES_",tag,".csv"), row.names = T)



##############################################
###  Run for "fasting" glucose (8-12 hrs)  ###
##############################################

#analysis <- analysis %>% mutate(glu_8to12hr = ifelse(fasting_hrs >= 8  & fasting_hrs <= 12, glu, NA))

#lapply(1:length(models.l), function(m) {
#  print_lm(exposure = "taste_diplos", outcome = "glu_8to12hr", label="Fasting Glucose", 
#           covariates = models.nofast.l[[m]], lm_trend = T, data = analysis)
#  })




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
  educ_isced.lab="Education Level", income_level.lab = "Income Level")

print_summary_table(data=analysis, var_strata = "n_24hr_plaus.lab", vars_to_summarize = decr_vars,
                    var_strata_order = n24HR_plaus.l,
                    factor_vars = c("taste_diplos", "smoke_level.lab", "physact_level",
                                    "alch_freq.lab", "educ_isced.lab", "income_level.lab"), p_print=T) %>% 
  write.csv(paste0("../data/processed/descr_sens_24HRplaus_",tag,".csv"))

print_summary_table(data=analysis, var_strata = "n_24hr.lab", var_strata_order = n24HR.l,
                    vars_to_summarize = c(ac.f="Assessment Center"),
                    factor_vars = c("ac.f"), p_print=T) %>% 
  write.csv(paste0("../data/processed/descr_sens_24HR_center_",tag,".csv"))


## regression
bind_rows(do.call(rbind.data.frame, lapply(4:6, function(m) { 
  print_lm(
    exposure = "taste_diplos", outcome = "glu", label= paste0("Has 24HR: RG~", names(models.l)[m]), 
    covariates = models.l[[m]], lm_trend = T, data = analysis %>% filter(n_24hr_plaus.lab==n24HR_plaus.l[[1]] )) } )),
  do.call(rbind.data.frame, lapply(4:6, function(m) { print_lm(
    exposure = "taste_diplos", outcome = "glu", label= paste0("Has 24HR: 2hr Glu~", names(models.nofast.l)[m]), 
    covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>% 
      filter(fast_cat=="0to2hr" & n_24hr_plaus.lab==n24HR_plaus.l[[1]] )) } )) ) %>% 
  write.csv(paste0("../data/processed/lm_sens_24hr_diplos_rg_2hg_m_dietF2Cses_",tag,".csv"), row.names = T)
  
models.no24hr.l = as.list(gsub("[+]FIB2CHO", "", models.l)) ; names(models.no24hr.l) = names(models.l)
models.no24hr.nofast.l = as.list(gsub("[+]FIB2CHO", "", models.nofast.l)) ; names(models.no24hr.nofast.l) = names(models.nofast.l)

bind_rows(do.call(rbind.data.frame, lapply(c(4,6), function(m) { 
  print_lm(
    exposure = "taste_diplos", outcome = "glu", label= paste0("No 24HR: RG~", names(models.no24hr.l)[m]), 
    covariates = models.no24hr.l[[m]], lm_trend = T, data = analysis %>% filter(n_24hr_plaus.lab==n24HR_plaus.l[[2]]) ) } )),
  do.call(rbind.data.frame, lapply(c(4,6), function(m) { print_lm(
    exposure = "taste_diplos", outcome = "glu", label= paste0("No 24HR: 2hr Glu~", names(models.nofast.l)[m]), 
    covariates = models.no24hr.nofast.l[[m]], lm_trend = T, data = analysis %>% 
      filter(fast_cat=="0to2hr" & n_24hr_plaus.lab==n24HR_plaus.l[[2]] ) ) } )) ) %>% 
  write.csv(paste0("../data/processed/lm_sens_no24hr_diplos_rg_2hg_m_dietF2Cses_",tag,".csv"), row.names = T)
  

###############################################################
##  Alternative variants for bitter taste: quinine/caffeine  ## 
###############################################################

# ========================
## Prepare bitter SNPs 
# ========================

# load dosage data for bitter snps
bitter_snps.df <- read.table("../data/processed/bitter_snps.raw", header = T) %>%
  mutate(rs2597979_G = ifelse(rs2597979_G >= 0 & rs2597979_G <0.5, 0, ifelse(rs2597979_G >=0.5 & rs2597979_G < 1.5, 1, 2))) %>%
  select(id=IID, "rs10772420_G", "rs2597979_G") 


## Add OTHER bitter SNPs
analysis <- analysis %>% left_join(bitter_snps.df, by="id")
bitter_snps <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "rs10772420_A", "rs2597979_G")
analysis <- analysis %>% mutate(rs10772420_A = 2-rs10772420_G)


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
) %>% write.csv(paste0("../data/processed/descr_sens_bitterSNPs_dietBitter_",tag,".csv"), row.names=F)



# ===========================================================
## Primary analysis with RG, A1c & Glu x fasting
# ===========================================================

# outcomes (formatted) 
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))


## Bitter SNPs & RG -----------------------
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu", covariates = models.l[[i]], label=paste0(snp, "_", names(models.l)[i]), 
             data = analysis, lm_trend = T)  })) 
  })) %>% write.csv(paste0("../data/processed/lm_sens_bittersnps_rg_",tag,".csv"), row.names = F)

## Bitter SNPs & 2hg -----------------------
# continuous dosage
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], label=paste0(snp, "_", names(models.nofast.l)[i]), 
             data = analysis %>% filter(fast_cat == "0to2hr"), lm_trend = T)  })) 
})) %>% write.csv(paste0("../data/processed/lm_sens_bittersnps_2hg_",tag,".csv"), row.names = F)

# categorical, per dominant allele 
analysis <- analysis %>% mutate(
  rs713598.a = descr_label_ordered.fun(., "rs713598_G", c("CC"=0, "CG"=1, "GG"=2)),
  rs1726866.a = descr_label_ordered.fun(., "rs1726866_G", c("AA"=0, "AG"=1, "GG"=2)),
  rs10246939.a = descr_label_ordered.fun(., "rs10246939_C", c("TT"=0, "TC"=1, "CC"=2)),
  rs2597979.a = descr_label_ordered.fun(., "rs2597979_G", c("AA"=0, "AG"=1, "GG"=2)),
  rs10772420.a = descr_label_ordered.fun(., "rs10772420_A", c("GG"=0, "GA"=1, "AA"=2))
  ) ; snps.a <- c("rs713598.a", "rs1726866.a", "rs10246939.a", "rs2597979.a","rs10772420.a")

do.call(rbind.data.frame, lapply(snps.a, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], label=paste0(snp, "_", names(models.nofast.l)[i]), 
             data = analysis %>% filter(fast_cat == "0to2hr"), lm_trend = T)  })) 
})) %>% write.csv(paste0("../data/processed/lm_sens_bittersnps_alleles_2hg_",tag,".csv"), row.names = T)


## Bitter SNPs & A1c -----------------------
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "hba1c_max", covariates = models.l[[i]], label=paste0(snp, "_", names(models.l)[i]), 
             data = analysis)  })) })) %>% 
  write.csv(paste0("../data/processed/lm_sens_bittersnps_a1c_",tag,".csv"), row.names = F)


## Bitter SNPs & Glu x fasting time -----------------------
fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l))
names(models.nofast.l) <- names(models.l)


## Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(bitter_snps, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
      as.data.frame(print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], 
                             label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                               filter(fast_cat == fast))) })) %>% 
      mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) })) %>%
  write.csv(paste0("../data/processed/lm_sens_bittersnps_gluXfast_m_diet_",tag,".csv"), row.names = F)


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

do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
    get_emm.fun(exposure = snp, outcome = "glu", covars = models.nofast.l[[4]], 
                reference = "0_0",
                data=analysis %>% filter(fast_cat == fast))$emm })) %>%
    mutate(fast = rep(fast_cat, each=3), .before=emmean) })) %>%
  write.csv(paste0("../data/processed/lm_emm_sens_bittersnps_gluXfast_m_diet_",tag,".csv"), row.names = F)


do.call(rbind.data.frame, lapply(bitter_snps_cat, function(snp) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(exposure = snp, outcome = "glu", covariates = models.nofast.l[[i]], label=paste0(snp, "_", names(models.nofast.l)[i]), 
             data = analysis %>% filter(fast_cat == "0to2hr"), lm_trend = F)  })) 
})) %>%  write.csv(paste0("../data/processed/lm_sens_bittersnps_cat_2hg_",tag,".csv"), row.names = T)


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
    Beta=c(-0.0790, -0.0510, -0.0570), SE=c(0.0022, 0.0021, 0.0021), P_value=c(0.0002571, 0.01228, 0.006042)
    )

## Add Variance, Weights & Inverse Weights
mr_tas2r38_glu <- mr_tas2r38_glu %>% mutate(
  Variance=(SE^2)) %>%
  mutate(W = c(1/Variance)) %>%
  mutate(IVW=Beta*W) %>%
  mutate(weights = 1/(SE^2)) %>% 
  mutate(weights_pct = round(weights/sum(weights)*100, 1)) %>%
  mutate(SNP_alleles=c("rs713598 (C>G)", "rs1726866 (G>A)", "rs10246939 (T>C)")) %>%
  mutate(Bitter_Allele=c("G", "G", "C")) %>%
  mutate(Diplotype=1:3)


## Calculate Pooled Beta & SE  
m.gen <- metagen(TE = Beta, seTE=SE, studlab = SNP_alleles, method.random.ci = "classic",
                 fixed=TRUE, random = FALSE, sm= "SMD", method.tau = "REML", data=mr_tas2r38_glu)

pdf("../output/meta_bittersnps_02hr_magic.pdf", height = 3, width=8)
forest(m.gen, layout = "RevMan5",
       sortvar = Diplotype, 
       #hetstat = F, 
       digits=3, 
       at=seq(-0.1,0.1,0.05),
       colgap.forest=c("10 mm"),
       col.square =  "#96A0B3", col.square.lines = "#96A0B3", col.inside = "black",
       col.diamond = "#435269", col.diamond.lines = "#435269",
       leftcols = c("SNP", "Bitter_Allele",  "w.fixed", "effect.ci"),
       leftlabs = c("TAS2R38\nVariant_EA", "Bitter\nAllele", "Weights\n(%)", "IVW Beta [95% CI]\nper Bitter Allele"),
       smlab = "2-hr Glucose (mmol/L),\nadjusted for BMI", text.fixed = "Pooled Effect Estimate",
       just.addcols.left = "left", just.studlab = "left", just = "left",
       fontsize = 10,
       squaresize = 0.75, lwd.diamond = 5, lwd.square = 1.25, 
       colgap.left = c("5 mm"))
       # fs.axis = 6, fs.heading = 8, fs.smlab = 8,
dev.off()


## Combine summary statitics for summary stats from UKB & **other** study (not significant for now...)



###################################################################
###  Check that covariates are capturing the expected variance  ###
###################################################################

## Suggested by JC (9/17)

analysis <- analysis %>%
  mutate(glu2hr = ifelse(fast_cat=="0to2hr", glu,NA)) %>%
  mutate(educ_isced.lab.f = factor(educ_isced.lab, levels = c(paste0("Level ",1:5))))

lifestyle_ses <- c(smoke_level.lab="Smoking", alch_freq.lab="Alcohol", physact_level="PA Level", 
                   income_level.lab="Income Level", educ_isced.lab="Education Level")

## Summarize lm associations of lifestyle/ses vars 
glu2hr_life_ses <- do.call(rbind.data.frame, lapply(1:length(lifestyle_ses), function(v) {
 make_pretty_lm(print_lm(exposure = names(lifestyle_ses)[v], outcome = "glu2hr", covariates = models.nofast.l$Base,
         data=analysis, label=lifestyle_ses[v], lm_trend = F)) }))

glu_life_ses <- do.call(rbind.data.frame, lapply(1:length(lifestyle_ses), function(v) {
  make_pretty_lm(print_lm(exposure = names(lifestyle_ses)[v], outcome = "glu", covariates = models.l$Base,
                          data=analysis, label=lifestyle_ses[v], lm_trend = F)) }))

a1c_life_ses <- do.call(rbind.data.frame, lapply(1:length(lifestyle_ses), function(v) {
  make_pretty_lm(print_lm(exposure = names(lifestyle_ses)[v], outcome = "hba1c", covariates = models.l$Base,
                          data=analysis, label=lifestyle_ses[v], lm_trend = F)) }))

glu2hr_life_ses %>% write.csv("../data/processed/sensitivity_lm_2hg_lifestyle_ses_base.csv")
glu_life_ses %>% write.csv("../data/processed/sensitivity_lm_glu_lifestyle_ses_base.csv")
a1c_life_ses %>% write.csv("../data/processed/sensitivity_lm_a1c_lifestyle_ses_base.csv")

#Adjusted R-squared
sapply(1:length(lifestyle_ses), function(v) {
  summary(lm(formula(paste0("hba1c~", names(lifestyle_ses)[v])),data=analysis))$adj.r.squared }
)


## Dig more into physical activity -------
pa_id <- fread("../data/processed/pa_id.csv") %>% 
  
  # if walking_frq = -2 (Unable to walk) --> Recode to 0
  mutate_at("walking_frq", ~ifelse(. == -2, 0, .)) %>%
  
  # if activity duration <10 min/day --> Recode to 0
  mutate(across(ends_with("dur"), ~ifelse(.<10,0,.))) %>%
  
  # replace Do not know (-1), Prefer not to answer (-3) or missing (NA) with median
  mutate(across(ends_with("dur") | ends_with("frq") , ~ ifelse(.<0 | is.na(.)==T, 0, .))) %>%
  
  # multipley minutes per day per activity by excess MET score, per activity type
  mutate(met_excess_walk = ((walking_dur * 2.3)/60)*walking_frq,
         met_excess_mod = ((moderate_dur * 3.0)/60)*moderate_frq,
         met_excess_vig = ((vigorous_dur * 7.0)/60)*vigorous_frq) %>%
  
  # calculate excess met-hr/wk
  mutate(physact_met_excess = met_excess_walk + met_excess_mod + met_excess_vig) %>%
  mutate(physact_met_excess = winsorize(physact_met_excess))

physact_met_excess_lvls <- quantile(pa_id$physact_met_excess, probs=seq(0,1,0.33), include.lowest=F)[-1]

pa_id <- pa_id %>% mutate(
  physact_level = case_when(
    physact_met_excess < physact_met_excess_lvls[1] ~ "1",
    physact_met_excess >= physact_met_excess_lvls[1] & 
      physact_met_excess < physact_met_excess_lvls[2] ~ "2",
    physact_met_excess >= physact_met_excess_lvls[2] ~ "3")
) %>%
  mutate(physact_level.lab = factor(physact_level, levels=c("1", "2", "3"), labels=c("Low", "Moderate", "High"))) %>%
  select(id, physact_eMET=physact_met_excess, physact_eMET_level.lab=physact_level.lab)

## Merge into analysis
analysis <- analysis %>% left_join(pa_id, by="id")



###################################
###  Run Sex-Stratified Models  ###
###################################

sex.l <- list("Female", "Male")

## NO significanat sex interactions
summary(lm(formula(paste0("glu2hr~taste_diplos*sex+", models.nofast.l$Base)), data=analysis))
summary(lm(formula(paste0("glu2hr~taste_diplos*sex+", models.nofast.l$Diet)), data=analysis))

View(do.call(rbind.data.frame, lapply(sex.l, function(x) {
  return(make_pretty_lm(print_lm(outcome="glu2hr", exposure="taste_diplos", covariates = gsub("sex[+]", "", models.nofast.l$Diet),
           label=paste0(x, " - Diet"), data=analysis %>% filter(sex==x), lm_trend = T)))
})))






## END_OF_SCRIPT




