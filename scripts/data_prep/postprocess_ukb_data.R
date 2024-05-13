## Post-processing

# load packages
library(tidyverse) 
library(data.table)

# load basic functions
source("../scripts/basic_functions.R", echo=F)



#####################
## Select ancestry ##
#####################

# system arguments
args = commandArgs(trailingOnly=TRUE)
ANC = args[1]

phenos_anc <- fread(paste0("../data/processed/ukb_phenos_unrelated_", ANC, ".csv"))
paste("Procecssing data for", ANC, "subset of N =", nrow(phenos_anc), "unrelated participants")



###########################
##  make variable groups ##
###########################

# genetic PCs ---------------
gPC_vars <- c(paste0("gPC", 1:10))


## basic phenotypes  ---------------
pheno_vars <- c("unrelated", "ancestry", "t2d_case.f", "t2d_med_any", "med_mets", "fasting_hrs",
                "ac", "ac_date", "age", "sex", gPC_vars, "smoke.lab", "alch_freq.lab", "ALC", 
                "pa_met_excess_lvl", "pa_met_excess", "income", "educ_level.lab", "educ_isced.lab", "educ_years",
                "glu", "hba1c_max", "bmi", "sbp", "dbp", "crp", "chol", "tg", "ldl" , "hdl")


# nutrient intakes ---------------
nutrient_vars <- c("TCALS", "CHO_pct", "FIBER", "FAT_pct", "MUFA_pct", "PUFA_pct", 
               "SFA_pct", "PRO_pct", "CHO2FIB", "FIB2CHO")

ffq_vars <- c("cooked_veg", "raw_veg", "fresh_fruit", "dried_fruit", 
              "bread_intake", "bread_type_white_vs_brown_or_whole", "bread_type.lab",
              "cereal_intake", "cereal_type_sugar_vs_any_bran", "cereal_type.lab", 
              "cheese", "milk_type_full_vs_low_or_nonfat", "milk_type.lab", 
              "spread_type_butter_vs_any_other", "oily_fish", "nonoily_fish", "procmeat", "poultry", 
              "beef", "lamb", "pork", "tea", "coffee", "coffee_type_decaf_vs_regular",
              "water", "addsalt.lab", "addsalt_freq_QT", "addsalt_always_often_vs_nrs", 
              "hotdrink_temp", "hotdrink_temp_hot_or_vhot_vs_warm")


paste("Done making variable groups")


######################################
##  compile supertaster haplotypes  ##
######################################
paste("Compiling supertaster haplotypes & SNP dosage")

## compile haplotype data -------------------------

haplo_names <- paste0(ANC,"_chunk",10:29,".csv")
haplo_path <- "../data/processed/supertasters"
haplos_id <- do.call(rbind.data.frame, lapply(as.list(paste0(haplo_path, "/", haplo_names)), function(file) fread(file)))

haplos_id <- haplos_id %>%
  mutate(taster_status = case_when(
    haplo_0 == "GGC" & haplo_1 == "GGC" ~ "Supertaster",
    haplo_0 == "GGC" & haplo_1 == "CAT" | haplo_1 == "GGC" & haplo_0 == "CAT" ~ "Taster",
    haplo_0 == "CAT" & haplo_1 == "CAT" ~ "Nontaster",
    TRUE ~ as.character(NA))) %>% 
  mutate(taster_status = ifelse(is.na(taster_status), "Other", taster_status)) %>%
  mutate(taster_status = factor(taster_status, levels = c("Nontaster", "Taster", "Supertaster", "Other")))


## merge haplotype & dosage data -------------------------

genos_id <- haplos_id %>%
  select(c("id", ends_with("_0"), ends_with("_1"), "diplo", "taster_status")) %>%
  left_join(fread("../data/processed/supertasters/supertasters.raw") %>%
              select(id=IID, starts_with("rs")), by = "id") %>%
  # recode to match taster allele 
  mutate("rs713598_G"=2-rs713598_C, 
         "rs10246939_C"=2-rs10246939_T) %>%
  select(-c("rs713598_C", "rs10246939_T"))


paste("Done compiling genotype data")
head(genos_id)


## Load & compile ukb phenotype data  -------------------------
paste("Compiling basic phenotypes")

phenos_id <- phenos_anc %>%   
select("id", all_of(pheno_vars), all_of(gPC_vars), all_of(nutrient_vars), all_of(ffq_vars)) %>%
  mutate(
    sex = case_when(sex == 1 ~ "Male", sex == 0 ~ "Female"),
    income_level.lab=case_when(
      income == "Prefer not to answer" ~ as.character(NA),
      income == "Do not know" ~ as.character(NA),
      income != "Prefer not to answer" & income != "Do not know" ~ as.character(income),
      TRUE ~ as.character(NA)),
    educ_level.lab = case_when(
      educ_level.lab == "Prefer not to answer" ~ as.character(NA),
      educ_level.lab != "Prefer not to answer" ~ as.character(educ_level.lab),
      TRUE ~ as.character(NA)),
    smoke_level.lab = case_when(
      smoke.lab == "No answer" ~ as.character(NA),
      smoke.lab != "No answer" ~ as.character(smoke.lab),
      TRUE ~ as.character(NA)),
    alch_freq.lab = case_when(
      alch_freq.lab == "Prefer not to answer" ~ as.character(NA),
      alch_freq.lab != "Prefer not to answer" ~ as.character(alch_freq.lab),
      TRUE ~ as.character(NA))) %>% 
  mutate(
    income_level.lab = factor(income_level.lab, 
                              levels = c("Less than 18,000", "18,000 to 30,999",
                                         "31,000 to 51,999", "52,000 to 100,000",  
                                         "Greater than 100,000"),
                              labels = c("lt_18000", "from_18000_to_30999",
                                         "from_31000_to_51999", "from_52000_to_100000",  
                                         "gt_100000")),
    #Edu levels & yrs based on: https://www.nature.com/articles/s41380-019-0596-9#MOESM1)
    educ_level.lab = factor(educ_level.lab, 
                            levels = c("College or university degree", # ~20yrs 
                                       "NVQ/HND or equivalent", # 2 of 3 years bachelor's degree ~19yrs
                                       "Other professional qualifications", # e.g., nursing degree, teaching degree ~ 15yrs
                                       "A/AS levels or equivalent", # 1 year bachelor's degree ~13yrs
                                       "O/GCSE levels or equivalent", # HS + Associates degree ~10yrs
                                       "CSEs or equivalent", # completed HS ~10yrs
                                       "None of the above")), #~7yrs
    smoke_level.lab = factor(smoke_level.lab, levels = c("Never", "Previous", "Current")),
    alch_freq.lab = factor(alch_freq.lab, ordered = T,
                           levels = c("Daily or almost daily",
                                      "3-4 per week", "1-2 per week", "1-3 per month",
                                      "Special occasions only", "Never")),
    pa_met_excess_lvl = factor(pa_met_excess_lvl, ordered = T,
                               levels = c("Low", "Moderate", "High")),
    tg_log = log(tg),
    hba1c=hba1c_max,
    fast_cat = case_when(fasting_hrs > 0 & fasting_hrs <= 2 ~ "0to2hr",
                         fasting_hrs == 3 ~ "3hr",
                         fasting_hrs == 4 ~ "4hr",
                         fasting_hrs == 5 ~ "5hr",
                         fasting_hrs >= 6 ~ "6+hr",
                         TRUE ~ as.character(NA)),
    fasting_gt6hr = case_when(fasting_hrs >= 6 ~ 1, fasting_hrs < 6 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt2hr = case_when(fasting_hrs >=2 ~ 1, fasting_hrs < 2 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt3hr = case_when(fasting_hrs >=3 ~ 1, fasting_hrs < 3 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt1hr = case_when(fasting_hrs >= 1 ~ 1, fasting_hrs < 1 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt4hr = case_when(fasting_hrs >= 4 ~ 1, fasting_hrs < 4 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt5hr = case_when(fasting_hrs >= 5 ~ 1, fasting_hrs < 5 ~ 0, TRUE ~ as.numeric(NA))
  )



# Build diet PCs for diet patterns -------------

vars_for_pca <- phenos_anc %>% select(
  id,
  cooked_veg_QT=cooked_veg, raw_veg_QT=raw_veg,
  fresh_fruit_QT=fresh_fruit, dried_fruit_QT=dried_fruit, 
  oily_fish_QT=oily_fish, nonoily_fish_QT=nonoily_fish,
  procmeat_QT=procmeat, poultry_QT=poultry, cheese_QT=cheese,
  beef_QT=beef, lamb_QT=lamb, pork_QT=pork, 
  bread_type_white_vs_brown_or_whole_BIN=bread_type_white_vs_brown_or_whole, 
  bread_intake_QT=bread_intake,
  milk_type_full_vs_low_or_nonfat_BIN=milk_type_full_vs_low_or_nonfat,
  cereal_type_sugar_vs_any_bran_BIN=cereal_type_sugar_vs_any_bran, cereal_intake_QT=cereal_intake,
  spread_type_butter_vs_any_other_BIN=spread_type_butter_vs_any_other,
  coffee_type_decaf_vs_regular_BIN=coffee_type_decaf_vs_regular, coffee_QT=coffee,
  tea_QT=tea, water_QT=water, 
  addsalt_always_often_vs_nrs_BIN=addsalt_always_often_vs_nrs,
  hotdrink_temp_hot_or_vhot_vs_warm_BIN=hotdrink_temp_hot_or_vhot_vs_warm,
) %>%
  mutate(across(where(is.numeric), function(i) remove_outliers.fun(i, SDs=5))) %>%
  #filter_at(vars(-id), any_vars(!is.na(.))) %>%
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x))  # Median impute


## Run PCA
diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA

# compile dietPC scores
dietPCs_id <- as.data.frame(cbind(id=vars_for_pca$id, diet_pcs$x))
colnames(dietPCs_id) <- c("id", paste0("diet", colnames(diet_pcs$x)))


# Save diet PCA results as csv (of factor loadings) & .rda (all outputs)
saveRDS(diet_pcs, file = paste0("../data/processed/diet/ukb_", ANC, "_dietPCs.rda"))

as.data.frame(diet_pcs$rotation) %>% 
  mutate(diet_var = gsub("_BIN", "", gsub("_QT", "", rownames(.)))) %>% 
  fwrite(paste0("../data/processed/diet/ukb_", ANC, "_dietPCloadings.csv"),col.names=T, row.names=F)


print("Done compiling basic phenotypes.")

###################################
## Merge data into ukb_processed ##
###################################

phenos_processed_id <- phenos_id %>%
  left_join(genos_id, by = "id") %>%
  left_join(dietPCs_id, by = "id")


### Exclusion criteria
#* No diabetes (Eastwood algorithm + HbA1c <5.7%)
#* Common TAS2R38 diplotype
#* ANC genetic ancestry
#* Complete data for glucose, genetic ancestry, covariates
#* Plausible total energy intakes (600-4000 kcal/d)

# function to identify missingness
summarise_noNA <- function(x, varname, data=phenos_processed_id) {
  return((data %>% select(c("id", x)) %>% mutate(noNA = ifelse(complete.cases(.), 1, 0)))$noNA) }

# add indicator variables for exclusion criteria
phenos_processed_id <- phenos_processed_id %>% 
  mutate(
    incl_t2d = ifelse(!is.na(t2d_case.f) & t2d_case.f == "Control", 1, 0),
    incl_taste = ifelse(is.na(taster_status) == T | taster_status == "Other", 0, 1),
    incl_compl_glu = ifelse((summarise_noNA("glu") == 1), 1, 0),
    incl_kcal = ifelse(TCALS <600 | TCALS > 4800 | is.na(TCALS)==T, 0, 1),
    incl_compl_fast = summarise_noNA("fasting_hrs"),
    incl_compl_covars = summarise_noNA(c("age", "sex", paste0("gPC", 1:10), "fasting_hrs", "bmi", "smoke_level.lab",
                                         "pa_met_excess_lvl", "alch_freq.lab", "dietPC1", "FIB2CHO"))) %>%
  mutate(incl_compl=ifelse(incl_compl_glu==1 & incl_compl_fast==1 & incl_compl_covars==1, 1, 0)) %>%
  mutate(N_Ancestry = ifelse(ancestry == ANC, 1, 0),
         N_contrl = incl_t2d,
         N_contrl_taste = ifelse(incl_t2d==1 & incl_taste==1, 1, 0),
         N_contrl_taste_kcal = ifelse(incl_t2d==1 & incl_taste==1 & incl_kcal==1, 1,0),
         N_contrl_taste_kcal_compl = ifelse(incl_t2d==1 & incl_taste==1 & incl_kcal==1 & incl_compl==1, 1,0))

## save as rda -------------
saveRDS(phenos_processed_id, file = paste0("../data/processed/ukb_analysis_", ANC, ".rda"))

## Remove continuous outliers with >5SD --------------
phenos_processed_id %>% 
  mutate(across(c(glu, c(paste0("gPC", 1:10)), TCALS, FIB2CHO, dietPC1), function(x) remove_outliers.fun(x, SDs=5))) %>% 
  filter(complete.cases(glu, gPC1, gPC2, gPC3, gPC4, gPC5, gPC6, gPC7, gPC8, gPC9, gPC10, TCALS, FIB2CHO, dietPC1)) %>%
  saveRDS(paste0("../data/processed/ukb_analysis_", ANC, "_rm5sd.rda"))


#EOF
