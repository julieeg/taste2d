# Script for calculating diet PCs from UKBB FFFQ data
# Prepared: 2025-01-14

# load packages
library(tidyverse)
library(data.table)



############################################
## Basic functions for preparing FFQ data ## 
############################################

## function to convert frequency values to servings/day (from KEW)
ffq_freq_to_sev <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}


## function to recode negative values as meaninginful
neg_to_num <- function(x) {
  case_when(
    x >= 0 ~ as.numeric(x), # -1 = "Do not know" ; -3 = "Prefer not to answer"
    x == -10 ~ 0.5, # -10 = "Less than 1 serving/day"
    TRUE ~ as.numeric(NA)
  )
}

winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


########################
##  Prepare FFQ data  ## 
########################

### FFQ fields ----------------------------------------

intake_fields <- c(cooked_veg = 1289, raw_veg = 1299,
                   fresh_fruit = 1309, dried_fruit = 1319,
                   bread_intake = 1438, bread_type = 1448, 
                   water = 1528, milk_type = 1418, spread_type = 1428,
                   spread_type_nonbutter=2654,
                   cereal_intake = 1458, cereal_type = 1468,
                   addsalt=1478, tea=1488, coffee = 1498, coffee_type = 1508,
                   hotdrink_temp = 1518) #**non_butter_spread_type

freq_fields <- c(oily_fish = 1329, nonoily_fish = 1339,
                 procmeat = 1349, poultry = 1359, cheese = 1408,
                 beef = 1369, lamb = 1379, pork = 1389)

ffq_fields<-c(intake_fields, freq_fields)

# compile ukb variable names (e.g., f.field.0.0)
ffq_vars <- setNames(paste0("f.", ffq_fields, ".0.0"), names(ffq_fields))


## Gather ids for participants that withdrew consent ---------------------------

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw/w27892_20241217.csv", what=character())


## Build FFQ ukb dataset ------------------------------------

ffq_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE) %>% 
  select(id=as.character(f.eid), ffq_vars) %>% 
  
  # Recode FFQ vars as integers
  mutate(across(names(freq_fields), ffq_freq_to_sev)) %>% # convert frequency vars to sev/day
  mutate(across(names(intake_fields), neg_to_num)) %>% # recode negative values for numeric vars
  
  # Recode vars in units/week --> units/day
  mutate(bread_intake = bread_intake / 7, # bread intake was provided in slices/week
         cereal_intake = cereal_intake / 7 # cereal intake was provided in bowls/week
  ) %>%
  
  # Create binary (from categorical) vars for PCA
  mutate(
    bread_type_white_vs_brown_or_whole_BIN = case_when(      
      bread_type == 1 ~ 1, bread_type == 2 | bread_type == 3 | 
        bread_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_full_vs_low_or_nonfat_BIN = case_when(       
      milk_type == 1 ~ 1, milk_type == 2 | milk_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_rare_never_BIN = case_when(
      milk_type == 6 ~ 1, milk_type != 6 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_butter_vs_any_other_BIN = case_when(          
      spread_type == 1 ~ 1, spread_type == 2 | spread_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_rare_never_BIN = case_when(
      spread_type == 0 ~ 1, spread_type != 0 ~ 0, TRUE ~ as.numeric(NA)),
    cereal_type_sugar_vs_any_bran_BIN = case_when(
      cereal_type == 5 ~ 1, cereal_type != 5 ~ 0, TRUE ~ as.numeric(NA)),
    coffee_type_decaf_vs_regular_BIN = case_when(              
      coffee_type == 1 ~ 1, coffee_type == 2 | coffee_type == 3 | 
        coffee_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    addsalt_freq_QT = addsalt,
    addsalt_always_often_vs_nrs_BIN = case_when(
      addsalt == 3 | addsalt == 4 ~ 1, addsalt == 1 | addsalt == 2 ~ 0, TRUE ~ as.numeric(NA)),
    hotdrink_temp_hot_or_vhot_vs_warm_BIN = case_when(
      hotdrink_temp == 1  ~ 1, hotdrink_temp == 3 | hotdrink_temp == 2 ~ 0, TRUE ~ as.numeric(NA))
  ) %>% 
  
  # Add lablels to FFQ categorical vars (for descriptive purposes)
  mutate(
    bread_type.lab = case_when(
      bread_type == 1 ~ "White", bread_type == 2 ~ "Brown", bread_type == 3 ~ "Wholemeal/Wholegrain",
      bread_type == 4 ~ "Other", bread_type == -1 ~ "Do not know", bread_type == -3 ~ "Prefer not to answer",
      TRUE ~ as.character(NA)),
    
    milk_type.lab = case_when(
      milk_type == 1 ~ "Full cream", milk_type == 2 ~ "Semi-skimmed", milk_type == 3 ~ "Skimmed",
      milk_type == 4 ~ "Soy", milk_type == 5 ~ "Other", milk_type == 6 ~ "Never/rarely have milk",
      milk_type == -1 ~ "Do not know", milk_type == -3 ~ "Prefer not to answer", 
      TRUE ~ as.character(NA)),
    
    spread_type.lab = case_when(
      spread_type == 1 ~ "Butter/spreadable butter", spread_type == 2 ~ "Flora Pro-Active/Benecol",
      spread_type == 3 ~ "Other spread/margarine", spread_type == 0 ~ "Never/rarely use spread",
      spread_type == -1 ~ "Do not know", spread_type == -3 ~ "Prefer not to answer"),
    
    spread_type_nonbutter.lab = case_when(
      spread_type_nonbutter == 4 ~ "Soft (tub) margarine", spread_type_nonbutter == 5 ~	"Hard (block) margarine",
      spread_type_nonbutter == 6 ~ "Olive oil based spread (eg: Bertolli)",
      spread_type_nonbutter == 7 ~ "Polyunsaturated/sunflower oil based spread (eg: Flora)",
      spread_type_nonbutter == 2 ~ "Flora Pro-Active or Benecol",
      spread_type_nonbutter == 8 ~ "Other low or reduced fat spread",
      spread_type_nonbutter == 9 ~ "Other type of spread/margarine", spread_type_nonbutter == -1 ~	"Do not know",
      spread_type_nonbutter == -3 ~ "Prefer not to answer"),
    
    cereal_type.lab = case_when(
      cereal_type == 1 ~ "Bran cereal (e.g. All Bran, Branflakes)", cereal_type == 2 ~ "Biscuit cereal (e.g. Weetabix)",
      cereal_type == 3 ~ "Oat cereal (e.g. Ready Brek, porridge)", cereal_type == 4 ~ "Muesli",
      cereal_type == 5 ~ "Other (e.g. Cornflakes, Frosties)", cereal_type == -1 ~ "Do not know",
      cereal_type == -3 ~ "Prefer not to answer"),
    
    coffee_type.lab = case_when(
      coffee_type == 1 ~ "Decaffeinated coffee (any type)", coffee_type == 2 ~ "Instant coffee",
      coffee_type == 3 ~ "Ground coffee (include espresso, filter etc)", coffee_type == 4 ~ "Other type of coffee",
      coffee_type == -1 ~ "Do not know", coffee_type == -3 ~ "Prefer not to answer"),
    
    addsalt.lab = case_when(
      addsalt == 1 ~ "Never/Rarely", addsalt == 2 ~ "Sometimes", addsalt == 3 ~ "Often",
      addsalt == 4 ~ "Always", TRUE ~ as.character(NA)),
    
    hotdrink_temp.lab = case_when(
      hotdrink_temp == 1 ~ "Very hot", hotdrink_temp == 2 ~ "Hot", hotdrink_temp == 3 ~ "Warm",
      TRUE ~ as.character(NA))
  )


# Build diet PCs for diet patterns -------------

vars_for_pca <- ffq_id %>% 
  mutate(id=as.character(id)) %>% 
  select(id,
  cooked_veg_QT=cooked_veg, raw_veg_QT=raw_veg,
  fresh_fruit_QT=fresh_fruit, dried_fruit_QT=dried_fruit, 
  oily_fish_QT=oily_fish, nonoily_fish_QT=nonoily_fish,
  procmeat_QT=procmeat, poultry_QT=poultry, cheese_QT=cheese,
  beef_QT=beef, lamb_QT=lamb, pork_QT=pork, 
  cereal_intake_QT=cereal_intake,
  bread_intake_QT=bread_intake,
  tea_QT=tea, water_QT=water, 
  bread_type_white_vs_brown_or_whole_BIN, 
  milk_type_full_vs_low_or_nonfat_BIN,
  cereal_type_sugar_vs_any_bran_BIN, 
  spread_type_butter_vs_any_other_BIN,
  coffee_type_decaf_vs_regular_BIN, coffee_QT=coffee,
  addsalt_always_often_vs_nrs_BIN,
  hotdrink_temp_hot_or_vhot_vs_warm_BIN) %>%
  
  # Replace missing values with medians
  mutate_at(vars(-"id"), function(x) ifelse(is.na(x), median(x, na.rm=T), x))  %>% 
  
  # Winsorize data to 5 SD
  mutate(across(where(is.numeric), function(i) winsorize(i, SDs=5))) %>%
  filter(complete.cases(.)==T)


## Run PCA
dietPCs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA

# compile dietPC scores
dietPCs_id <- as.data.frame(cbind(id=vars_for_pca$id, dietPCs$x))
colnames(dietPCs_id) <- c("id", paste0("diet", colnames(dietPCs$x)))

# Save diet PCA results as .rda (all outputs) & csv (of factor loadings) -------
saveRDS(dietPCs, file = paste0("../data/processed/diet/ukb_dietPCdata_MA_forMSG_20250113.rda"))
left_join(
  vars_for_pca %>% left_join(ffq_id %>% select("id", ends_with(".lab")), by = "id"), 
  dietPCs_id, by = "id") %>% 
  write.csv(file = paste0("../data/processed/diet/ukb_dietPCscores_MA_forMSG_20250113.csv"))



###############################################
## WaterFall plot of diet PC Factor Loadings ##
###############################################

# Descriptive diet labels
diet_labels_descriptive <- c(
  raw_veg_QT="Raw vegetable intake",
  cooked_veg_QT="Cooked vegetable intake",
  fresh_fruit_QT="Fresh fruit intake",
  dried_fruit_QT="Dried fruit intake",
  oily_fish_QT="Oily fish intake",
  nonoily_fish_QT="Non-oily fish intake",
  procmeat_QT="Processed meat intake",
  poultry_QT="Poultry intake",
  cheese_QT="Cheese intake",
  beef_QT="Beef intake",
  lamb_QT="Lamb intake",
  pork_QT="Pork intake",
  bread_type_white_vs_brown_or_whole_BIN="Prefer white>brown/whole bread",
  bread_intake_QT="Bread intake",
  milk_type_full_vs_low_or_nonfat_BIN="Prefer full>low/no-fat milk",
  cereal_type_sugar_vs_any_bran_BIN="Prefer sugary>bran cereal",
  coffee_type_decaf_vs_regular_BIN="Prefer decaf>caffeinated coffee",
  cereal_intake_QT="Cerealintake",
  spread_type_butter_vs_any_other_BIN="Prefer butter>other spreads",
  coffee_QT="Coffee intake",
  tea_QT="Tea intake",
  water_QT="Water intake",
  addsalt_always_often_vs_nrs_BIN="Always/often add salt",
  hotdrink_temp_hot_or_vhot_vs_warm_BIN="Prefer hot>warm drink temps"
)

palette_waterfall = c("#888363","#C5C1A5", "#96A0B3", "#435269")

# prepare data for visualization
dietPC.loadings <- dietPCs$rotation %>% 
  as.data.frame() %>%
  mutate(Diet=diet_labels_descriptive[rownames(.)], .before=PC1) %>%
  arrange(PC1)
dietPC.mat <- as.matrix(dietPC.loadings[,-1])
rownames(dietPC.mat) <- dietPC.loadings$Diet

# =========================
## Plot dietPC waterfall
# =========================
plot_dietPC_waterfall.fun <- function(nPCs=10) {
  dietPC.loadings %>%
    select(Diet, c(paste0("PC", 1:nPCs))) %>%
    mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
    pivot_longer(-Diet) %>% 
    mutate(name=factor(name, levels=c(paste0("PC", 1:nPCs)), labels=c(paste0("Diet PC", 1:nPCs)) )) %>%
    mutate(direction = factor(ifelse(value>0.2, "Positive/Major", ifelse(value<0.2 & value>0, "Positive/Minor", 
                                                                         ifelse(value<0 & value>-0.2, "Negative/Minor", "Negative/Major"))),
                              levels=c("Positive/Major", "Positive/Minor", "Negative/Minor", "Negative/Major"))) %>%
    ggplot(aes(x=value, y=Diet, fill=direction)) + 
    facet_wrap(~name, ncol=nPCs/2) + 
    geom_col() + ylab("") + xlab("Rotated Factor Loadings") +
    geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values=palette_waterfall, name = "Factor Loadings") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(color="black"),
          legend.position = "top",
          strip.text = element_text(face="bold", size=16),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          axis.title.x = element_text(size=14),
          legend.key.size = unit(1,"cm"),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16, face="bold"))
}

pdf("../data/processed/diet/ukb_dietPCloadings_MA_waterfall24_forMSG_20250113.pdf", height=18, width=36)
plot_dietPC_waterfall.fun(24)
dev.off()

## END OF SYNTAX


