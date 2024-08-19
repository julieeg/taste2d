## Post-processing

# load packages
library(tidyverse) 
library(data.table)
library(vcfR)

# load basic functions
source("../scripts/pantry.R", echo=F)



#####################
## Select ancestry ##
#####################

# system arguments
args = commandArgs(trailingOnly=TRUE)
ANC = args[1]

# load phenotype data
phenos <- fread(paste0("../data/processed/ukb_phenos_unrelated_", ANC, ".csv"))
paste("Procecssing data for", ANC, "subset of N =", nrow(phenos), "unrelated participants")



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


cat(" \n Done making variable groups. \n ")


######################################
##  compile supertaster haplotypes  ##
######################################

paste("Compiling supertaster haplotypes & SNP dosage")


## compile haplotype data -------------------------

haplos_id <- fread("../data/processed/supertasters/supertaster_phased_haps.csv")

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


## Clean diplotype variables ------------------------

## Recode diplotypes so palindromic heterozygotes are coded the same
## example: CAT/GGC AND GGC/CAT are both coded as CAT/GGC

recode_diplos <- function(diplos) {
  haplo_pairs <- strsplit(diplos, "/")

  # Sort each pair and then combine them back into a string
  recoded <- sapply(haplo_pairs, function(pair) {
    sorted_pair <- sort(pair)
    paste(sorted_pair, collapse = "/")
  }) ; return(recoded) 
}

genos_id <-  genos_id %>% 
  mutate(diplos_pal=recode_diplos(diplo)) 


## Make variable for common diplotypes -------------

diplos_common <- which(prop.table(table(genos_id$diplos_pal))>0.001)
genos_id <- genos_id %>% mutate(
  diplos_common = ifelse(diplos_pal %in% names(diplos_common), diplos_pal, NA)
)


## Recode alleles as amino acids -----------------------

genos_id <- genos_id %>% 
  mutate(diplos_common_AA=gsub("CAT", "AVI", gsub("CAC", "AVV", gsub("CGT", "AAI", 
                               gsub("CGC", "AAV", 
                                     gsub("GGC", "PAV", gsub("GGT", "PAI", 
                                          gsub("GAT", "PVI", gsub("GAC", "PVV", 
                                                diplos_common)))))))))


## Create variable for canonical diplotypes (AVI/AVI, AVI/PAV, PAV/PAV)
taste_diplos <- c("AVI/AVI", "AVI/PAV", "PAV/PAV")
genos_id <- genos_id %>%
  mutate(taste_diplos = ifelse(diplos_common_AA %in% taste_diplos, diplos_common_AA, NA)) 


paste("Done compiling genotype data! ")
head(genos_id)



#########################################
##  Load & compile ukb phenotype data  ## 
#########################################

paste("Compiling basic phenotypes ... ")

phenos_id <- phenos %>%   
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
      TRUE ~ as.character(NA)),
    alch_heavydrinker = as.factor(ifelse(alch_freq.lab == "Daily or almost daily" | alch_freq.lab == "3-4 per week", 1, 0)),
    alch_lightdrinker = as.factor(ifelse(alch_freq.lab == "Special occasions only" | alch_freq.lab == "Never", 1, 0)),
    #Ref: https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=100026
    kcal_plaus = ifelse(sex == "Male" & TCALS >= 600 & TCALS < 4780.1 | 
                          sex == "Female" & TCALS >= 600 & TCALS < 4302.2, TCALS, NA),
    fasting_hrs_lt24 = ifelse(fasting_hrs <= 24, fasting_hrs, NA)) %>%
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
    physact_level = factor(pa_met_excess_lvl, ordered = T,
                               levels = c("Low", "Moderate", "High")),
    tg_log = log(tg),
    hba1c=hba1c_max,
    fast_cat = case_when(fasting_hrs_lt24 >= 0 & fasting_hrs_lt24 <= 2 ~ "0to2hr",
                         fasting_hrs_lt24 == 3 ~ "3hr",
                         fasting_hrs_lt24 == 4 ~ "4hr",
                         fasting_hrs_lt24 == 5 ~ "5hr",
                         fasting_hrs_lt24 >= 6 ~ "6+hr",
                         TRUE ~ as.character(NA)),
    fasting_gt6hr = case_when(fasting_hrs_lt24 >= 6 ~ 1, fasting_hrs_lt24 < 6 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt2hr = case_when(fasting_hrs_lt24 >=2 ~ 1, fasting_hrs_lt24 < 2 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt3hr = case_when(fasting_hrs_lt24 >=3 ~ 1, fasting_hrs_lt24 < 3 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt1hr = case_when(fasting_hrs_lt24 >= 1 ~ 1, fasting_hrs_lt24 < 1 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt4hr = case_when(fasting_hrs_lt24 >= 4 ~ 1, fasting_hrs_lt24 < 4 ~ 0, TRUE ~ as.numeric(NA)),
    fasting_gt5hr = case_when(fasting_hrs_lt24 >= 5 ~ 1, fasting_hrs_lt24 < 5 ~ 0, TRUE ~ as.numeric(NA))
  ) %>%
  mutate(
    taste_diplos.num = case_when(
      taste_diplos == "AVI/AVI" ~ 0,
      taste_diplos == "AVI/PAV" ~ 1,
      taste_diplos == "PAV/PAV" ~ 2,
      TRUE ~ as.numeric(NA)),
    taste_alleles = (rs713598_G+rs1726866_G+rs10246939_C)
    )



# Build diet PCs for diet patterns -------------

vars_for_pca <- phenos %>% select(
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
  hotdrink_temp_hot_or_vhot_vs_warm_BIN=hotdrink_temp_hot_or_vhot_vs_warm) %>%
  
  # Replace missing values with medians
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x))  %>% 
  
  # Winsorize data to 5 SD
  mutate(across(where(is.numeric), function(i) winsorize(i, SDs=5))) %>%
  filter(complete.cases(.)==T)

cat("Dimensions of dietary data \n")
dim(vars_for_pca)


## Run PCA
diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA

# compile dietPC scores
dietPCs_id <- as.data.frame(cbind(id=vars_for_pca$id, diet_pcs$x))
colnames(dietPCs_id) <- c("id", paste0("diet", colnames(diet_pcs$x)))

# Save diet PCA results as csv (of factor loadings) & .rda (all outputs)
saveRDS(diet_pcs, file = paste0("../data/processed/diet/ukb_", ANC, "_dietPCs.rda"))

diet_id <- left_join(vars_for_pca, dietPCs_id, by = "id")

#as.data.frame(diet_pcs$rotation) %>% 
#  mutate(diet_var = gsub("_BIN", "", gsub("_QT", "", rownames(.)))) %>% 
#  fwrite(paste0("../data/processed/diet/ukb_", ANC, "_dietPCloadings.csv"),col.names=T, row.names=F)


print("Done compiling basic phenotypes & ancestry-specific diet patterns!")



###################################
## Merge data into ukb_processed ##
###################################

phenos_processed_id <- phenos_id %>%
  left_join(genos_id, by = "id") %>%
  left_join(diet_id, by = "id")


### Exclusion criteria
#** No diabetes (Eastwood algorithm + HbA1c <5.7%)
#** ANC genetic ancestry
#** Complete data for glucose, genetic ancestry, covariates
#   - sensitivity: Plausible total energy intakes (600-4300 [F] or 4800 [M] kcal/d)
#   - sensitivity: Canonical TAS2R38 diplotypes

# function to identify missingness
summarise_noNA <- function(x, varname, data=phenos_processed_id) {
  return((data %>% select(c("id", x)) %>% mutate(noNA = ifelse(complete.cases(.), 1, 0)))$noNA) }

# add indicator variables for exclusion criteria
phenos_processed_id <- phenos_processed_id %>% 
  mutate(
    incl_geno = ifelse(is.na(haplo_0) == F & is.na(haplo_1) == F, 1, 0),
    incl_t2d = ifelse(!is.na(t2d_case.f) & t2d_case.f == "Control", 1, 0),
    incl_compl_glu = ifelse(is.na(glu)==T, 0, 1),
    incl_taste = ifelse(is.na(taster_status) == T | taster_status == "Other", 0, 1),
    incl_kcal = ifelse(is.na(kcal_plaus)==T, 0, 1)) %>%
  mutate(
    incl_compl_fast = ifelse(is.na(fast_cat)==T, 0, 1),
    incl_compl_covars = summarise_noNA(c("age", "sex", paste0("gPC", 1:10), "fast_cat", "bmi", "smoke_level.lab",
                                 "physact_level", "alch_freq.lab", "dietPC1"))) %>%
  mutate(incl_compl=ifelse(incl_compl_glu==1 & incl_compl_covars==1, 1, 0)) %>%
  mutate(N_Ancestry = ifelse(ancestry == ANC, 1, 0),
         N_geno = ifelse(is.na(haplo_0) == F & is.na(haplo_1) == F,1,0),
         N_contrl = incl_t2d,
         N_contrl_compl = ifelse(incl_t2d==1 & incl_compl == 1, 1, 0),
         N_contrl_compl_kcal = ifelse(incl_t2d==1 & incl_compl==1 & incl_kcal==1, 1, 0),
         N_contrl_compl_kcal_taste = ifelse(incl_t2d==1 & incl_compl == 1 & incl_kcal==1 & incl_taste==1, 1, 0))


# ===============================================
## Compile & save ukb_analysis_ANC.rda
# ===============================================

phenos_processed_id %>% 
  saveRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda"))


cat("Done compiling ukb_analysis_ANC.rda dataset!")

#EOF


