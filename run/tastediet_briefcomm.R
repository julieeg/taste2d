
############
## Set Up ## 
############

lapply(c("tidyverse", "data.table", "paletteer", "ggpubr", "ggtext", #"gplots", "R3port", 
         "xtable", "knitr", "tinytex", "emmeans"),  library, character.only = TRUE)

lapply(list.files("../../pantry/functions/", full.names = T), source)
palettes_pantry = palettes
source("../scripts/pantry.R", echo=F)


###################################
##   Load data for EUR ancestry  ##
###################################

# system arguments
EUR=c("European"="EUR")
tag="v6"

# load postprocessed dataframe
postprocessed <- readRDS("../data/processed/ukb_postprocessed_EUR_20241227.rda")

analysis <- readRDS(paste0("../data/processed/ukb_analysis_EUR.rda"))
range(analysis$glu) #: 2.116 - 8.002

#load dietPCs
dietPCs.dat <- readRDS(paste0("../data/processed/diet/ukb_EUR_dietPCs.rda"))
dietPC_emm <- read.csv(paste0("../data/processed/descr_tab_emm_dietPC_diplos_sexageac_", tag,".csv"))
dietPC_pstat <- read.csv(paste0("../data/processed/descr_tab_teststat_dietPC_diplos_sexageac_",tag,".csv"))

taste_diplos.labs <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")
taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")
taste_diplos_005.labs <- c(taste_diplos.labs, "AVI/AAV", "AAV/PAV")

## Add glucose values with mg/dl
analysis <- analysis %>% mutate(
  glu_mgdl = glu*18.018,
  glu2hr_mgdl = glu2hr*18.018)

## Load preference traits
prefs_id <- fread("../../tasteprefs/data/raw/ukbb_app27892_diet_preference_11062024.csv") %>%
  rename(id=eid, survey_date=p20750, survey_duration=p20751) %>%
  filter(!is.na(survey_date)) %>% # filter out participants without preference data (N=320071)
  mutate(across(starts_with("p"), ~as.numeric(
    case_when(
      . %in% c(1:10) ~ as.character(.),
      . == "extremely dislike" ~ "1",
      . == "neither like nor dislike" ~ "5",
      . == "extremely like" ~ "9",
      . %in% c("never tried", "do not wish to answer") ~ "NA",
      TRUE ~ as.character(NA) ) ))
  ) 

codebook <-readxl::read_xlsx("../../tasteprefs/data/raw/food_pref_vars.xlsx") %>%
  mutate(var_name = gsub("Liking_for", "pref", gsub("/", "", gsub("[(]", "", gsub("[)]", "", gsub("-","", gsub(" ", "_", description)))))))

prefs_id <- prefs_id %>%
  rename_with(
    ~with(codebook, var_name[match(., paste0("p", fieldID))]),
    starts_with("p")
  )

analysis <- left_join(analysis, prefs_id, by="id")

sex.l <- c("Female", "Male")




######################################################
###   Descriptive Analysis of Diet Traits - TOTAL  ###
######################################################

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
baseModel <- paste0("age+sex+ac.f+",paste0("gPC",1:10,collapse = "+"))

emm_diplos_diet.l <- lapply(diet_traits, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = baseModel, data=analysis) 
}) ; names(emm_diplos_diet.l) = diet_traits
do.call(rbind.data.frame, lapply(diet_traits, function(d) {
  emm_diplos_diet.l[[d]]$emm} )) %>% mutate(DietLabel=rep(diet_labels, each=3)) %>%
  mutate(model="Demographic") %>%
  fwrite("../data/processed/briefcomm/emm_means_diplos_diet_mBase_total.csv")


## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  Diet=diet_traits,
  F_stat=sapply(1:24, function(i) {emm_diplos_diet.l[[i]]$anv$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_diet.l[[i]]$anv$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(diet_traits, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+",baseModel)), data=analysis)))[2,]
    })) %>% as.data.frame() %>% mutate(Diet=diet_traits, DietLabel=diet_labels) %>% 
      rename(Beta=Estimate, SE='Std. Error', "T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="Diet") %>%
  fwrite("../data/processed/briefcomm/emm_betas_diplos_diet_mBase_total.csv")


# ============================================================
## Calculate EMMs for dPCs by diplotype, adjusting for age+sex
# ============================================================

dietPCs <- paste0("dietPC", 1:24)

emm_diplos_dPCs.l <- lapply(dietPCs, function(d) {
  get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, covars = baseModel, data=analysis) 
}) ; names(emm_diplos_dPCs.l) = dietPCs
do.call(rbind.data.frame, lapply(dietPCs, function(dPC) {
  emm_diplos_dPCs.l[[dPC]]$emm} )) %>%
  mutate(model="Demographic") %>%
  fwrite("../data/processed/briefcomm/emm_means_diplos_dietPCs_mBase_total.csv", row.names = F)

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
  fwrite("../data/processed/briefcomm/emm_betas_diplos_dietPCs_mBase_total.csv")





###############################################################
###   Descriptive Analysis of Diet Traits - SEX-STRATIFIED  ###
###############################################################

descr_dietTraits_sex <- do.call(rbind.data.frame, lapply(sex.l, function(x) {
  print_summary_table(
    vars_to_summarize = diet_to_summarise, 
    var_strata = "taste_diplos", var_strata_order = taste_diplos,
    p_print = T, p_adjust = c("age", "ac.f"), digits = c(3,3,4), 
    data=analysis %>% filter(sex==x)) %>%
    mutate(Sex=x)
})) 

descr_dietTraits_sex %>% fwrite("../data/processed/briefcomm/descr_dietTraits_sex.csv")


beta_diplo_diet_sex <- do.call(rbind.data.frame, lapply(sex.l, function(x) {
  do.call(rbind.data.frame, lapply(1:24, function(d) {
    print_lm(exposure = "taste_diplos.num", outcome=names(diet_labels)[d], covariates = "taste_diplos.num", label=diet_labels[d],
             data=analysis %>% filter(sex == x))
  })) %>% mutate(sex=x)
}))

beta_diplo_diet_sex %>% fwrite("../data/processed/briefcomm/beta_diplo_diet_sex.csv")


## Descriptive analysis of diet PCs across diplotypes
dietPCs_labs <- paste0("Diet PC", 1:24) ; names(dietPCs_labs) <- dietPCs
descr_dietPCs <- print_summary_table(
  vars_to_summarize = dietPCs_labs, var_strata = "taste_diplos", 
  var_strata_order = taste_diplos, p_print = T, 
  p_adjust = c("age", "sex", "ac.f"), digits = c(3,3,4), data=analysis) 



# ========================================================================
## Calculate EMMs for Dietary Traits by diplotype, adjusting for age+sex
# ========================================================================

baseModel_nosex <- paste0("age+ac.f+",paste0("gPC",1:10,collapse = "+"))

emm_diplos_diet_sex.l <- lapply(sex.l, function(x) {
  emm_diplos_diet_sex <- lapply(diet_traits, function(d) {
    get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = d, label="Demographic",
                covars = baseModel_nosex, data=analysis %>% filter(sex == x)) 
    }) ; names(emm_diplos_diet_sex) = diet_traits ; list(
      emm_diplos_diet_sex, do.call(rbind.data.frame, lapply(diet_traits, function(d) {
        emm_diplos_diet_sex[[d]]$emm %>% mutate(sex=x)} )) %>% mutate(DietLabel=rep(diet_labels, each=3)) 
      ) }) ; names(emm_diplos_diet_sex.l) <- sex.l

do.call(rbind.data.frame, emm_diplos_diet_sex.l) %>% write.csv("../data/processed/briefcomm/emm_means_diplos_dietPCs_mBase_sex.csv")


## For which dPCs do scores differ significantly by TAS2R38 diplotype
names(emm_diplos_diet_sex.l) <- sex.l
emm_diplos_diet_sex <- do.call(rbind.data.frame, lapply(sex.l, function(x) {
  dat.l <- emm_diplos_diet_sex.l[[x]]
    as.data.frame(cbind(
    Diet=diet_traits,
    F_stat=sapply(1:24, function(i) {dat.l[[1]][[i]]$anv$`F value`[1]}),
    P_val_f=sapply(1:24, function(i) {dat.l[[1]][[i]]$anv$`Pr(>F)`[1]}))) %>%
    left_join(
      do.call(rbind, lapply(diet_traits, function(d) {
        coef(summary(lm(formula(paste0(d, "~taste_diplos.num+", baseModel_nosex)), data=analysis %>% filter(sex==x))))[2,]
      })) %>% as.data.frame() %>% mutate(Diet=diet_traits, DietLabel=diet_labels, sex=x) %>% 
        rename(Beta=Estimate, SE='Std. Error', "T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
      by="Diet") 
  })) ; emm_diplos_diet_sex %>%
  write.csv("../data/processed/emm_betas_diplos_dietPCs_mBase_sex.csv")
  


# ============================================================
## Calculate EMMs for dPCs by diplotype, adjusting for age+sex
# ============================================================

dietPCs <- paste0("dietPC", 1:24)

## LEFT OFF HERE ##

emm_diplos_dietPCs_sex.l <- lapply(sex.l, function(x) {
  emm_diplos_dietPCs_sex <- lapply(dietPCs[1:3], function(dPC) {
    get_emm.fun(exposure = "taste_diplos", reference="AVI/AVI", outcome = dPC, label="Demographic",
                covars = baseModel_nosex, data=analysis %>% filter(sex == x)) 
  }) ; names(emm_diplos_dietPCs_sex) = dietPCs[1:3] ; list(
    emm_diplos_dietPCs_sex, do.call(rbind.data.frame, lapply(dietPCs[1:3], function(dPC) {
      emm_diplos_dietPCs_sex[[dPC]]$emm %>% mutate(sex=x)} )) %>% mutate(DietLabel=rep(dPC, each=3)) 
  ) }) ; names(emm_diplos_dietPCs_sex.l) <- sex.l


## For which dPCs do scores differ significantly by TAS2R38 diplotype
names(emm_diplos_dietPCs_sex.l) <- sex.l
emm_diplos_dietPCs_sex <- do.call(rbind.data.frame, lapply(sex.l, function(x) {
  dat.l <- emm_diplos_dietPCs_sex.l[[x]]
  as.data.frame(cbind(
    DietPC=1:24,
    F_stat=sapply(1:24, function(i) {dat.l[[1]][[i]]$anv$`F value`[1]}),
    P_val_f=sapply(1:24, function(i) {dat.l[[1]][[i]]$anv$`Pr(>F)`[1]}))) %>%
    left_join(
      do.call(rbind, lapply(dietPCs, function(d) {
        coef(summary(lm(formula(paste0(d, "~taste_diplos.num+", baseModel_nosex)), data=analysis %>% filter(sex==x))))[2,]
      })) %>% as.data.frame() %>% mutate(DietPC=1:24, DietLabel=diet_labels, sex=x) %>% 
        rename(Beta=Estimate, SE='Std. Error', "T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
      by="Diet") 
})) ; emm_diplos_diet_sex %>%
  write.csv("../data/processed/emm_betas_diplos_dietPCs_mBase_sex.csv")



## For which dPCs do scores differ significantly by TAS2R38 diplotype
as.data.frame(cbind(
  DietPC=1:24,
  F_stat=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv$`F value`[1]}),
  P_val_f=sapply(1:24, function(i) {emm_diplos_dPCs.l[[i]]$anv$`Pr(>F)`[1]}))) %>%
  left_join(
    do.call(rbind, lapply(dietPCs, function(d) {
      coef(summary(lm(formula(paste0(d, "~taste_diplos.num+age+sex+ac.f")), data=analysis)))[2,3:4]
    })) %>% as.data.frame() %>% mutate(DietPC=1:24) %>% rename("T_stat" = `t value`, "P_val_t" = `Pr(>|t|)`),
    by="DietPC") #%>%
#write.csv(paste0(outDir, "descr_tab_teststat_dietPC_diplos_sexageac_",tag,"_temp.csv"), row.names = F)

emm_diplos_dPCs.df <- fread("../data/processed/descr_tab_emm_dietPC_diplos_sexageac_v6.csv")



## EOF





