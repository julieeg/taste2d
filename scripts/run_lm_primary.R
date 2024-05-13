# load functions
library(tidyverse) ; library(data.table)

#command args
args=commandArgs(trailingOnly = T)
ANC=args[1]

# load analysis dataframe
analysis <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, "_rm5sd.rda")) %>% 
  filter(N_contrl_taste_kcal_compl == 1)
dim(analysis)

########################
### Primary Analysis ###
########################

# Build function to run simple lm with covariate adjustment
print_lm <- function(exposure="taster_status", outcome="glu", covariates=m, label=M, data=analysis) {
  mod <- lm(formula(paste0(outcome, "~", exposure, "+", covariates)), data)
  out<-matrix(NA, length(mod$xlevels[[1]]), 5, dimnames = list(paste0(label, "_", mod$xlevels$taster_status), c("beta", "se", "p", "f", "f_p")))
  out[2:nrow(out),c(1:3)] <- summary(mod)$coef[2:3,c(1:2,4)] ; out[1,4:5] <-c(anova(mod)[1,4], anova(mod)[1,5])
  return(out)
}

# outcomes (formatted) 
glucose_var.l <- list(variables = list(glu = "Glucose, mmol/L"))
hba1c_var.l <- list(variables = list(hba1c_max = "HbA1c, %"))

# covariates
m1_base = "age+sex+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+gPC7+gPC8+gPC9+gPC10+fast_cat+bmi"
m2_lifestyle = paste0(m1_base, "+smoke_level.lab+pa_met_excess_lvl+alch_freq.lab")
m3_carbqual = paste0(m2_lifestyle, "+FIB2CHO")
m4_dietpattern = paste0(m3_carbqual, "+dietPC1+dietPC2+dietPC3+dietPC4+dietPC5")

models.l = list("Base"= m1_base, "Lifestyle"=m2_lifestyle, "Carb.Quality"=m3_carbqual, 
                "Diet.Patterns"=m4_dietpattern)

########################
## Taster Status & RG ##
########################

## Run main effects of taster status on RG
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(outcome = "glu", covariates = models.l[[i]], label=names(models.l)[i])  
})) %>% mutate(model=rownames(.), .before=beta) %>%
  write.csv(paste0("../output/main/", ANC, "_lm_rg_taste.csv"), row.names = F)


#########################
## Taster Status & A1c ##
#########################

## Run main effects of taster status on A1C
do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
  print_lm(outcome = "hba1c_max", covariates = models.l[[i]], label=names(models.l)[i])  
})) %>% mutate(model=rownames(.), .before=beta) %>%
    write.csv(paste0("../output/main/", ANC, "_lm_a1c_taste.csv"), row.names = F)


#############################################
## Taster Status & Glucose by Fasting Time ##
#############################################

fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))

# covariates without fasting
models.nofast.l <- gsub("[+]fast_cat", "", models.l)
names(models.nofast.l) <- names(models.l)

## Run main effects of taster status on glucose by fasting time
do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.l), function(i) {
    print_lm(outcome = "glu", covariates = models.nofast.l[[i]], label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% filter(fast_cat == fast))  
  })) %>% mutate(model=rownames(.), .before=beta)
})) %>% write.csv(paste0("../output/main/", ANC, "_lm_gfast_taste.csv"), row.names = F)


##EOF







