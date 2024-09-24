# Title: Primary analysis
# Date Updated: 08-27-2024


## load functions
source("../scripts/pantry/pantry.R")


#command args
args=commandArgs(trailingOnly = T)
ANC="EUR"
tag="v5" # differentiated by addign adjustment for assessment center (factor)


## load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda")) %>%
  filter(n_geno==1 & n_t2d_contrl==1 & n_complete==1 & n_fast_24hr==1) %>%
  ## Criteria added: v2 to v3
  filter(find_outliers.fun(glu)==0) %>%
  mutate(glu2hr=ifelse(fast_cat=="0to2hr",glu,NA)) %>%
  mutate(ALC_drinkers = ifelse(alch_freq.lab == "Never", NA, ALC))

dim(analysis) # N=241378


# "Common" TAS2R38 hapotypes & diplotypes (>0.005)
common_diplos <- names(which(prop.table(table(analysis$diplos_common))>0.005))
common_diplos_AA <- names(which(prop.table(table(analysis$diplos_common_AA))>0.005))
taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")

# create output directories
outDir=paste0("../data/processed/")


## Add Assessment Center ----------------------------------------------------------
ac_labs <- list("Barts"=11012, "Birmingham" = 11021, "Bristol" =	11011, "Bury" =	11008, 
             "Cardiff" =	11003, "Cheadle (revisit)" =	11024, "Croydon" =	11020, 
             "Edinburgh" =	11005, "Glasgow" = 11004, "Hounslow" = 11018, "Leeds" = 11010,
             "Liverpool"=11016, "Manchester"=11001, "Middlesborough"=11017, "Newcastle" =11009, 
             "Nottingham"=11013, "Oxford"=11002, "Reading"=11007, "Sheffield"=11014, "Stockport (pilot)"=10003,
             "Stoke"=11006, "Swansea"=	11022,"Wrexham" =11023, "Cheadle (imaging)"=11025,
             "Reading (imaging)"=11026, "Newcastle (imaging)" =11027, "Bristol (imaging)"=11028)

analysis <- analysis %>% mutate(ac.f = descr_label.fun(., "ac", ac_labs))


################################################
## Descriptive analysis of participant cohort ##
################################################

#Diplotype distributions by assessment center
print_summary_table(data=analysis, vars_to_summarize = c("taste_diplos" = "TAS2R38 Diplotypes"),
                    var_strata = "ac.f", p_types = "descriptive", factor_vars = "taste_diplos") %>%
  write.csv(paste0("../data/processed/descr_diplo_taste_center_",tag,".csv"))
 
## Demographic, behavioral & risk factors
vars_to_summarize <- c(
  age="Age, years", sex="Sex", bmi="BMI, kg/m2",
  smoke_level.lab="Smoking", physact_level="Physical Activity, MET/wk",
  alch_freq.lab="Alcohol Frequency", ALC_drinkers="Alcohol intake (among drinkers)",
  alch_1to4d.lab="Alcohol Frequency, 1-4 per week", alch_heavydrinker = "Heavy alcohol drinker",
  educ_isced.lab="Education Level", income_level.lab = "Income Level",
  #sbp="SBP, mmHg", dbp="DBP, mmHg",
  tg="Triglyceride, mmol/L", ldl="LDL, mmol/L", hdl="HDL, mmol/L",
  coffee_QT="Coffee, cups/day", tea_QT="Tea, cups/day", addsalt_3lvl.lab="Added salt",
  hba1c_max="HbA1C, %" #glu="Glucose, mmol/L",   glu2hr="0-2hr Glucose, mmol/L"
)


# ========================================
## Run over all diplotypes with >0.005 
# ========================================

diplo_sets <- c(diplos_common_AA="diplo_AA_005", taste_diplos="diplo_AA_taste", taster_status="diplo_AA_taste_other")

analysis <- analysis %>% mutate(
  diplos_common_AA = factor(diplos_common_AA, levels=c("AVI/AVI", "AVI/PAV", "PAV/PAV", "AVI/AAV", "AAV/PAV")))
  
for (i in diplo_sets) {
  
  diplo_strata_order =  levels(analysis[[names(i)]])
  
  rbind.data.frame(
    
    # print summary table with digits: 1=means, 0=pct, 4=p-values
    print_summary_table(
      vars_to_summarize = vars_to_summarize, 
      var_strata = names(i), var_strata_order = c(diplo_strata_order),
      p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(1,0,4), data=analysis) ,
    
    # Add table for sex-stratified alcohol intake (g/day)
    do.call(rbind.data.frame, lapply(c("Female", "Male"), function(x){
      print_summary_table(vars_to_summarize = c(ALC_drinkers=paste0(x, ":Alcohol intake (drinkers)")),
                          var_strata = names(i), var_strata_order = c(diplo_strata_order),
                          p_print = T, p_adjust = c("age", "ac.f"), digits = c(3,0,4), 
                          data=analysis %>% filter(sex==x))
      })),
    
    # Add tables for traits with 0 digits (BP)
    print_summary_table(vars_to_summarize = c(sbp="SBP, mmHg", dbp="DBP, mmHg"),
                        var_strata = names(i), var_strata_order = c(diplo_strata_order),
                        p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(0,0,4), data=analysis),  
    
    # Add tables for traits with 3 digits (glucose & 2hr glucose)
    print_summary_table(vars_to_summarize = c(glu="Glucose, mmol/L (3 digits)", glu2hr="0-2hr Glucose, mmol/L"),
                        var_strata = names(i), var_strata_order = c(diplo_strata_order),
                        p_print = T, p_adjust = c("age", "sex", "ac.f"), digits = c(3,0,4), data=analysis)
    ) %>% 
    
    write.csv(paste0("../data/processed/descr_table1_",i, "_", tag,".csv"))
}
  


###############################################
###   Preliminary Analysis of Diet Traits   ###
###############################################

# ============================================================
## Calculate EMMs for dPCs by diplotype, adjusting for age+sex
# ============================================================

dietPCs <- paste0("dietPC", 1:24)

emm_diplos_dPCs.l <- lapply(dietPCs, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = "age+sex+ac.f") 
}) ; names(emm_diplos_dPCs.l) = dietPCs
do.call(rbind.data.frame, lapply(dietPCs, function(dPC) {
  emm_diplos_dPCs.l[[dPC]]$emm} )) %>%
  write.csv(paste0(outDir, "descr_emm_dietPC_diplos_sexagecenter_",tag,".csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  DietPC=1:24,
  F_stat=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv.p$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv.p$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(dietPCs, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex+ac.f")), data=analysis)))[2,3:4]
      })) %>% as.data.frame() %>% mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="DietPC"
    ) %>%
  write.csv(paste0(outDir, "descr_pftstat_dietPC_diplos_sexagecenter_",tag,".csv"), row.names = F)


## For which dPCs do scores differ significantly by TAS2R38 diplotype - NO adjustment
 do.call(rbind, lapply(dietPCs, function(d) {
   coef(summary(lm(formula(paste0(d, "~taste_diplos.num")), data=analysis)))[2,3:4]
   })) %>% as.data.frame() %>% 
   mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`) %>%
   write.csv(paste0(outDir, "descr_pftstat_dietPC_diplos_noadj_",tag,".csv"), row.names = F)
  

######################################################
###   Primary Analysis with Canonical Diplotypes   ###
######################################################

## Outcomes
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))


## Model covariates
#m3_ses = paste0(m2_lifestyle, "+educ_isced.lab+income_level.lab")

m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+ac.f+fast_cat"
m2_bmi = paste0(m1_base, "+bmi")
m3_lifestyle = paste0(m2_bmi, "+smoke_level.lab+physact_level+alch_freq.lab")
m4_diet = paste0(m3_lifestyle, "+dietPC3+dietPC4+dietPC6+dietPC11+dietPC12+dietPC13+dietPC15+dietPC16+dietPC17")
models.l = list("Base"= m1_base, "BMI"=m2_bmi, "Lifestyle"=m3_lifestyle, "Diet.Patterns"=m4_diet)

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)

# ===============================
## Canonical diplotypes & RG 
# ===============================

## Categorical diplotypes
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.l)[m], 
           covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>%
  mutate(mod_exp=rownames(.), .before=n) %>% separate(mod_exp, sep="_", into = c("model", "exposure")) %>%
  write.csv(paste0(outDir, "lm_diplo_rg_", tag, ".csv"), row.names = F)



## Continuous diplotypes
do.call(rbind.data.frame, lapply(1:length(models.l), function(m){
  print_lm(outcome = "glu", exposure = "taste_diplos.num", covariates = models.l[[m]], label=names(models.l)[m])
})) %>% write.csv(paste0(outDir, "lm_diplo_num_rg_", tag, ".csv"), row.names = F)


## estimated marginal mean glu in each model
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.l[[m]], reference = "AVI/AVI", data=analysis)
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  out <- out %>% mutate(model = rep(names(models.l[[m]]), each=3), .before=emmean) 
  return(out) } )) %>%
  write.csv(paste0(outDir, "lm_emm_diplo_glu_", tag, ".csv"), row.names = F)



# ===============================
## Canonical diplotypes & HbA1c
# ===============================
# categorical diplotypes
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label=names(models.l)[m], 
           covariates = models.l[[m]], lm_trend = T, data = analysis) } )) %>% 
  mutate(mod_exp=rownames(.), .before=n) %>% separate(mod_exp, sep="_", into = c("model", "exposure")) %>%
  write.csv(paste0(outDir, "lm_diplo_a1c_", tag, ".csv"), row.names = F)

## estimated marginal mean hba1c in each model
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "hba1c", covars = models.l[[m]], reference = "AVI/AVI", data=analysis)
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  out <- out %>% mutate(model = rep(names(models.l[[m]]), each=3), .before=emmean) 
  return(out) } )) %>%
  write.csv(paste0(outDir, "lm_emm_diplo_a1c_",tag,".csv"), row.names = F)



# =====================================
## Canonical diplotypes & 2hr Glucose
# =====================================
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.nofast.l)[m], 
           covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>%
             filter(fast_cat == "0to2hr")) })) %>% 
  mutate(mod_exp=rownames(.), .before=n) %>% separate(mod_exp, sep="_", into = c("model", "exposure")) %>%
  write.csv(paste0(outDir, "lm_diplo_2hg_",tag,".csv"), row.names = F)
  
# continuous diplotypes
do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(m){
  print_lm(outcome = "glu2hr", exposure = "taste_diplos.num", covariates = models.nofast.l[[m]], label=names(models.nofast.l)[m])
})) %>% write.csv(paste0(outDir, "lm_diplo_num_2hg_", tag, ".csv"), row.names = F)


## estimated marginal mean 2hr-glucose in each model
do.call(rbind.data.frame, lapply(1:length(models.l), function(m) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.nofast.l[[m]], reference = "AVI/AVI", 
                     data=analysis %>% filter(fast_cat=="0to2hr"))
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  out <- out %>% mutate(model = rep(names(models.nofast.l[[m]]), each=3), .before=emmean) 
  return(out) } )) %>%
  write.csv(paste0(outDir, "lm_emm_diplo_2hg_",tag,".csv"), row.names = F)


# Fasting 1hr
do.call(rbind.data.frame, lapply(1:4, function(m) {
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.nofast.l)[m], 
           covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>%
             filter(fasting_hrs == 2)) })) %>% 
  mutate(mod_exp=rownames(.), .before=n) %>% separate(mod_exp, sep="_", into = c("model", "exposure")) 


## ======================================================
##  Canonical diplotypes & Glu over ALL fasting winodows 
## ======================================================

fast_cat <- c("0to2hr", "3hr", "4hr", "5hr", "6+hr")
fast_cat.l <- as.list(fast_cat)

models.nofast.l <- as.list(gsub("[+]fast_cat", "", models.l)) ; names(models.nofast.l) <- names(models.l)


## TAS2R38 diplotpes & Glu by fasting time ------------
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:4, function(m) {
    print_lm(exposure = "taste_diplos", outcome = "glu", label=paste0(fast, "_", names(models.nofast.l)[m]), 
             covariates = models.nofast.l[[m]], lm_trend = T, data = analysis %>%
               filter(fast_cat == fast)) })) })) %>% 
  mutate(mod_exp=rownames(.), .before=n) %>% separate(mod_exp, sep="_", into = c("fast", "model", "exposure")) %>%
    write.csv(paste0(outDir, "lm_diplo_gluXfast_",tag,".csv"), row.names = F)


## estimated marginal mean glucose for all fasting times in DIET model
do.call(rbind.data.frame, lapply(fast_cat.l, function(f) {
  emm <- get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.nofast.l$Diet, reference = "AVI/AVI", 
                     data=analysis %>% filter(fast_cat==f))
  out <- emm$emm %>% mutate(anv.p=emm$anv.p[1,5])
  out <- out %>% mutate(fast = rep(f, each=3), .before=emmean) 
  return(out) } )) %>%
  write.csv(paste0(outDir, "lm_emm_diplo_gluXfast_m_diet_",tag,".csv"), row.names = F)


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



##EOF


