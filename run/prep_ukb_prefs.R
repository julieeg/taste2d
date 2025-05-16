# load packages
library(tidyverse)
library(data.table)

#######################
##  Basic functions  ##
#######################

## Add descriptive labels 
descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}



##############################################
##  Load demographic & lifestyle variables  ##
##############################################

### Basic phenotypes -----------------------------------------------

all_prepared_phenos <- fread("../../taste2d/data/processed/ukb_phenos.csv")

# participants that withdrew consent
withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw//withdraw27892_459_20240527.txt", what=character())



####################################################
##  Food preference traits & Taste preference PCs ##
####################################################

prefs <- fread("../data/raw/ukbb_app27892_diet_preference_11062024.csv") %>%
  rename(id=eid, survey_date=p20750, survey_duration=p20751) %>%
  filter(!is.na(survey_date)) %>% # filter out participants without preference data (N=320071)
  mutate(across(starts_with("p"), ~as.numeric(
    case_when(
    . %in% c(1:10) ~ as.character(.),
    . == "xtremely dislike" ~ "1",
    . == "neither like nor dislike" ~ "5",
    . == "extremely like" ~ "10",
    . %in% c("never tried", "do not wish to answer") ~ "NA",
    TRUE ~ as.character(NA) ) ))
  )
  

## Taste-based food preference groups
codebook <-readxl::read_xlsx("../data/raw/food_pref_vars.xlsx") 
tastes <- c("bitter", "sweet", "salt", "sour", "umami")

derive_tastepref_pca.fun <- function(taste) {
  
  food_for_pca <- codebook %>% pivot_longer(c(tastes)) %>%
    filter(name == taste & value == 1) %>% select(varID, description)
  
  data_for_pca <- prefs %>% select(id, all_of(food_for_pca$varID)) %>% 
    mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x))
    
  pca_summary <- prcomp(select(data_for_pca, -id), scale.=T)# Run PCA
  pca_scores <- as.data.frame(cbind(id=data_for_pca$id, pca_summary$x)) %>%
    rename_with(., ~ifelse(. == "id", "id", paste0(taste, "_", .)))
  
  # Add meaningful trait labels to PCA summary
  pca_loadings <- pca_summary$rotation %>% as.data.frame() %>%
    mutate(varID=rownames(.)) %>% 
    left_join(food_for_pca %>% rename(Food_Preference_Trait=description), by = "varID") %>%
    select(varID, Food_Preference_Trait, starts_with("PC"))
  
  ## IF PC1 is negatively loading, flip pca_scores & pca_loadings to relfect HIGHER preferences
  if(sum(ifelse(pca_loadings$PC1<0,1,0)) >= nrow(pca_loadings)/2) {
    
    cat(taste, "PCA loadings correspond to lower prefrence --> flip \n")
    pca_loadings <- pca_loadings %>% mutate(across(starts_with("PC"), ~ .*(-1)))
    pca_scores <- pca_scores %>% mutate(across(-id, ~ .*(-1)))
    
    summary = list(scores=pca_scores, loadings=pca_loadings, summary=pca_summary, flipped=T)
    
  } else {
    
    cat(taste, "PCA loadings correspond to higher prefrence --> good to go! \n")
    summary = list(scores=pca_scores, loadings=pca_loadings, summary=pca_summary, flipped=F)
    
  }
  
  return(summary)
  
} 

# Create, name & same 
taste_pref_pca.l <- lapply(tastes, derive_tastepref_pca.fun)
names(taste_pref_pca.l) <- tastes
taste_pref_pca.l %>% saveRDS("../data/processed/taste_pref_pca_20250117.rda")

## Create data.frame with taste preference PC scores for merging 
taste_pref_scores_id <- lapply(taste_pref_pca.l, function(x) x[["scores"]]) %>%
  reduce(inner_join, by = "id")

## Save csv with ALL taste_PC scores
taste_pref_scores_id %>% fwrite("../data/processed/taste_pref_scores_20250117.csv")


###############################################
##  Merge & save dataset for post-processing ##
###############################################

full_join(all_prepared_phenos, taste_pref_scores, by = "id") %>% fwrite("../data/processed/taste_pref_phenos_full_20250117.csv")

left_join(taste_pref_scores, all_prepared_phenos, by = "id") %>% fwrite("../data/processed/taste_pref_phenos_20250117.csv")


###############################################
##  Post-processing ##
###############################################


  





