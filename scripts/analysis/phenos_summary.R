## Post-processing


lapply(c("tidyverse", "data.table", "paletteer", "ggpubr",
         "R3port", "xtable", "knitr", "tinytex"),  library, character.only = TRUE)


## command args
args <- commandArgs(trailingOnly = TRUE)
ANC=args[1]


source("../scripts/pantry.R", echo=F)


analysis <- readRDS(paste0("../data/processed/ukb_analysis_",ANC,".rda")) %>%
  filter(N_contrl_taste_kcal_compl=="1") %>%
  mutate(taster_status = factor(taster_status, levels=c("Nontaster", "Taster","Supertaster"))) 



###############################
##   Load data by ancestry   ##
###############################

# system arguments
#args = commandArgs(trailingOnly=TRUE)
ANCs = c("European"="EUR", "Afrian"="AFR", "Admixed American"="AMR", "Central/South Asian"="CSA", 
         "East Asian"="EAS", "Middle Eastern"="MID")

#load data by ancestry
analysis_rda.l <-lapply(ANCs, function(ANC) {
  readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda"))
}) ; names(analysis_rda.l)<-ANCs


#load dietPCs by ancestry
dietPC_rot.l <-lapply(ANCs, function(ANC) {
  readRDS(paste0("../data/processed/diet/ukb_", ANC, "_dietPCs.rda"))
}) ; names(dietPC_rot.l)<-ANCs


# N per ANC total sample
sapply(analysis_rda.l, function(dat) dat %>% nrow()) %>% as.data.frame()


# ===================================================================
### Exclusion criteria
# * No diabetes (Eastwood algorithm + HbA1c <5.7%)
# * Common TAS2R38 diplotype
# * [ANC] genetic ancestry
# * Complete data for glucose, genetic ancestry, covariates
# * Plausible total energy intakes (600-4000 kcal/d)
# ===================================================================

# Sample inclusions by ANC
sapply(analysis_rda.l, function(dat) dat %>% select(starts_with("N_")) %>% 
         reframe(N_Ancestry=n_pct(N_Ancestry, level="1"),
                 N_contrl=n_pct(N_contrl, level="1"),
                 N_contrl_taste=n_pct(N_contrl_taste, level="1"),
                 N_contrl_taste_kcal=n_pct(N_contrl_taste_kcal, level="1"),
                 N_contrl_taste_kcal_compl=n_pct(N_contrl_taste_kcal_compl, level="1"))
)


# Taster status by ANC
sapply(analysis_rda.l, function(dat) dat %>% select(taster_status) %>% table()) %>% t()




################################
##   Descriptive statistics   ##
################################

# ======================================
# Taste diplotypes 
# ======================================

plot_Diplo.fun <- function(ANC) {
  
  anc <- analysis_rda.l[[ANC]]
  
  # gather diplotypes with >0.05% frequency
  which(prop.table(table(anc$diplos_pal))*100>0.05)
  
  # Add yscale
  yscale <- max(prop.table(table(anc$diplos_pal))*100)*1.075

  anc %>% 
    filter(incl_t2d==1 & incl_compl==1 & incl_kcal==1) %>%
    reframe(n=table(diplos_pal), pct=prop.table(table(diplos_pal)), diplos_pal=names(table(diplos_pal))) %>%
    arrange(-n) %>% 
    filter(pct>0.0005) %>%
    mutate(
      color = factor(ifelse(diplos_pal == "CAT/CAT", "Nontaster",
                     ifelse(diplos_pal %in% c("CAT/GGC", "GGC/CAT"), "Taster",
                            ifelse(diplos_pal == "GGC/GGC", "Supertaster", "Other"))),
                     levels=c("Nontaster","Taster","Supertaster","Other"))) %>%
    mutate(diplos_order=factor(diplos_pal, levels=diplos_pal[order(n, decreasing=T)])) %>%
    ggplot(aes(x=diplos_order, y=pct*100, fill=color)) + ggtheme +
    geom_bar(stat = "identity") + 
    scale_x_discrete(name = "\n TAS2R38 Diplotypes") +
    ylim(0,yscale) +
    ylab("Frequency (%)") + theme_bw() +
    scale_fill_manual(values = Taster_palette_Labs, name="TAS2R38\nDiplotypes") + 
    theme(axis.text.x = element_text(size=8, angle=25, hjust=1, color = "black"),
          axis.text.y = element_text(size=8, color = "black"),
          axis.title = element_text(color = "black", face = "bold")) +
    geom_text(aes(label=n), vjust=-0.15, size=3) + 
    ggtitle(paste0("TAS2R38 diplotypes among N = (", nrow(analysis_rda.l[[ANC]]),") ", ANC, " individuals"))
}

lapply(ANCs, plot_Diplo.fun)

plot_Diplo.fun("EUR")



# ==================================================================
## Summary tables by taster status - for INCLUDED participants
# ==================================================================

## Subset to included participants
analysis.l <- lapply(analysis_rda.l, function(d) d %>% filter(N_contrl_taste_kcal == 1))

## Demographic, behavioral & risk factors
tabl_descr_basic.fun <- function(ANC) {
  as.data.frame(analysis_rda.l[[ANC]] %>% 
                  filter(incl_t2d==1, incl_kcal==1, incl_compl==1) %>%
  filter(complete.cases(taster_status)) %>%
  group_by(taster_status) %>%
  reframe(
    Age=mean_sd(age,d=1),
    Female=n_pct(sex, level = "Female"),
    BMI=mean_sd(bmi,d=1),
    Smoke_Current=n_pct(smoke_level.lab, level="Current"),
    Smoke_Previous=n_pct(smoke_level.lab, level="Previous"),
    Smoke_Never=n_pct(smoke_level.lab, level="Never"),
    PA_Low=n_pct(pa_met_excess_lvl, level="Low"),
    PA_Moderate=n_pct(pa_met_excess_lvl, level="Moderate"),
    PA_High=n_pct(pa_met_excess_lvl, level="High"),
    Alchohol_Daily=n_pct(alch_freq.lab, level="Daily or almost daily"),
    Alcohol_3to4wk=n_pct(alch_freq.lab, level="3-4 per week"),
    Alcohol_1to2wk=n_pct(alch_freq.lab, level="1-2 per week"), 
    Alcohol_1to3mnth=n_pct(alch_freq.lab, level="1-3 per month"), 
    Alcohol_occsions=n_pct(alch_freq.lab, level="Special occasions only"), 
    Alcohol_never=n_pct(alch_freq.lab, level="Never"),
    Alcohol_g=mean_sd(ALC, d=2),
    Glucose=mean_sd(glu),
    HbA1c=mean_sd(hba1c_max),
    SBP=mean_sd(sbp),
    DBP=mean_sd(dbp),
    Triglyceride=mean_sd(tg),
    LDL=mean_sd(ldl),
    HDL=mean_sd(hdl)) %>% 
  t()) 
}

tabl_descr_basic.fun("EUR") %>% write.csv("../data/processed/EUR_tab_descr.csv")
lapply(ANCs, tabl_descr_basic.fun)


## Additional P-values 
analysis <- analysis %>% mutate(taster_num = as.numeric(factor(taster_status, labels=c(0,1,2))))
table(analysis$taster_status, analysis$taster_num)

# function to calculate p-values across taste diplotypes
p_cont_bytaste_sexageadj.fun <- function(cont_var, data=analysis) {
  d <- data %>% select(age, sex, taster_status, taster_num, TCALS, y=all_of(cont_var))
  anv.P <- c((anova(lm(y~taster_status+age+sex+TCALS, d)))$`Pr(>F)`[1],NA,NA)
  lm.P <- c(NA, summary(lm(y~taster_status+age+sex+TCALS, d))$coef[2:3,4])
  lm.P.trend <- c(summary(lm(y~taster_num+age+sex+TCALS, d))$coef[2,4], NA, NA)
  d %>% group_by(taster_status) %>% 
    reframe(msd=mean_sd(y)) %>%
    mutate(lm.P=lm.P, lm.P.trend=lm.P.trend, Anova.P=anv.P) %>%
    mutate(diet=rep(cont_var,3), .before=msd) %>%
    return()
}


## calculate p-values for all diet traits & PCs
do.call(rbind.data.frame, lapply(names(diet_labels), function(i) p_cont_bytaste_sexageadj.fun(i)))
lapply(dietPCs, function(i) p_cont_bytaste_sexageadj.fun(i))


## Glucose/A1c by fasting time 
tabl_descr_gluA1c_fast.fun <- function(ANC) {
  bind_rows(
  analysis.l[[ANC]] %>% filter(complete.cases(taster_status)) %>%
    group_by(fast_cat, taster_status) %>%
    reframe(Glucose=mean_sd(glu)) %>%
    pivot_wider(names_from="taster_status", values_from="Glucose") %>%
    mutate(Biomarker="Glucose", .before = fast_cat),
  analysis.l[[ANC]] %>% filter(complete.cases(taster_status)) %>%
    group_by(fast_cat, taster_status) %>%
    reframe(HbA1c=mean_sd(hba1c_max)) %>%
    pivot_wider(names_from="taster_status", values_from="HbA1c") %>%
    mutate(Biomarker="HbA1c", .before = fast_cat))
}

tabl_descr_gluA1c_fast.fun("EUR")

lapply(ANCs, tabl_descr_gluA1c_fast.fun)



## waterfall plot of Factor loading
plot_dietPC_waterfall.fun <- function(ANC, nPCs=10) {
  as.data.frame(dietPC_rot.l[[ANC]]$rotation) %>%
    mutate(Diet=diet_labels[gsub("_QT", "", gsub("_BIN", "", rownames(.)))]) %>%
    select(Diet, c(paste0("PC", 1:nPCs))) %>%
    mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
    pivot_longer(-Diet) %>% 
    mutate(name=factor(name, levels=c(paste0("PC", 1:10)) )) %>%
    mutate(name=gsub("PC", "Diet Pattern PC", name)) %>%
    mutate(direction = factor(ifelse(value>0.2, "Positive/Major", ifelse(value<0.2 & value>0, "Positive/Minor", 
                                     ifelse(value<0 & value>-0.2, "Negative/Minor", "Negative/Major"))),
                              levels=c("Positive/Major", "Positive/Minor", "Negative/Minor", "Negative/Major"))) %>%
    ggplot(aes(x=value, y=Diet, fill=direction)) + 
    facet_wrap(~name, ncol=5) + 
    geom_col() + ylab("") + xlab("Rotated Factor Loadings") +
    geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values=c(palettes$NatComms[c(1,5,6,4)]), name = "Factor Loadings") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(color="black"),
          legend.position = "top",
          strip.text = element_text(face="bold", size=6),
          axis.text.y=element_text(size=8),
          axis.text.x=element_text(size=8),
          axis.title.x = element_text(size=8),
          legend.key.size = unit(0.25,"cm")) #+
      #ggtitle(paste0("Rotated factor loadings for dietPCs derived among ", ANC, " participants"))
}

ggsave(plot_dietPC_waterfall.fun("EUR",5), filename="../output/EUR_plot_dietPCs_top5.pdf", height=4, width=8)


## Summary table by taster status
tabl_descr_diet.fun <- function(ANC) {
  as.data.frame(
  analysis.l[[ANC]] %>% filter(complete.cases(taster_status)) %>% 
    group_by(taster_status) %>%
    reframe(
      Energy_kcal=mean_sd(TCALS),
      Carbohydrate_pct = mean_sd(CHO_pct),
      Fat_pct = mean_sd(FAT_pct),
      Protein_pct = mean_sd(PRO_pct), 
      Alcohol_g=mean_sd(ALC),
      Diet_Pattern_PC1 = mean_sd(dietPC1),
      Diet_Pattern_PC2 = mean_sd(dietPC2),
      Diet_Pattern_PC3 = mean_sd(dietPC3),
      Diet_Pattern_PC4 = mean_sd(dietPC4),
      Diet_Pattern_PC5 = mean_sd(dietPC5),
      Diet_Pattern_PC6 = mean_sd(dietPC6),
      Diet_Pattern_PC7 = mean_sd(dietPC7),
      Diet_Pattern_PC8 = mean_sd(dietPC8),
      Diet_Pattern_PC9 = mean_sd(dietPC9),
      Diet_Pattern_PC10 = mean_sd(dietPC10)) %>% 
    t()) %>% 
  mutate(P_Ftest=c(NA, sapply(dietVars, function(y) {
    anova(lm(formula(paste0(y, "~taster_status+age+sex")), 
             data=analysis %>% filter(taster_status != "Other")))$`P`[1]} )))
}

#EOF

# ==========================================================================================
#HOLD
# ==========================================================================================

analysis_all.l$ALL %>%
  group_by(diplos_common) %>%
  reframe(
    n=n(),
    raw_veg=mean_sd(raw_veg_QT, d=1),
    cooked_veg=mean_sd(cooked_veg_QT, d=1),
    fresh_fruit=mean_sd(fresh_fruit_QT, d=1),
    dried_fruit=mean_sd(dried_fruit_QT, d=1),
    oily_fish=mean_sd(oily_fish_QT, d=1),
    nonoily_fish=mean_sd(nonoily_fish_QT,d=1),
    procmeat=mean_sd(procmeat_QT,d=1),
    poultry=mean_sd(poultry_QT, d=1),
    cheese=mean_sd(cheese_QT, d=1),
    beef=mean_sd(beef_QT, d=1),
    lamb=mean_sd(lamb_QT, d=1),
    pork=mean_sd(pork_QT, d=1),
    bread_intake=mean_sd(bread_intake_QT, d=1),
    coffee=mean_sd(coffee),
    tea=mean_sd(tea_QT, d=1),
    water=mean_sd(water_QT, d=1),
    choose_white_bread=n_pct(bread_type_white_vs_brown_or_whole_BIN, level="1"),
    choose_fullfat_milk=n_pct(milk_type_full_vs_low_or_nonfat_BIN, level="1"),
    choose_sugary_cereal=n_pct(cereal_type_sugar_vs_any_bran_BIN, level="1"),
    choose_butter_spread=n_pct(spread_type_butter_vs_any_other_BIN, level="1"),
    choose_decaf_coffee=n_pct(coffee_type_decaf_vs_regular_BIN, level="1"),
    frequently_add_salt=n_pct(addsalt_always_often_vs_nrs_BIN, level="1"),
    prefer_hot_drinks=n_pct(hotdrink_temp_hot_or_vhot_vs_warm_BIN, level="1")
  ) %>%
  t() %>%
  mutate(P_Ftest=c(NA, sapply(dietVars, function(y) {
    anova(lm(formula(paste0(y, "~taster_status+age+sex")), 
             data=analysis %>% filter(taster_status != "Other")))$`P`[1]} )))


# ==========================================================================================
## ARCHIVE 2
# ==========================================================================================

msd_dietByanc <- do.call(rbind.data.frame, lapply(ANCs, function(ANC) {
  tab <- analysis.l[[ANC]] %>% select(c(id, raw_veg_QT, cooked_veg_QT, fresh_fruit_QT, dried_fruit_QT, coffee_QT, tea_QT, addsalt_freq_QT)) %>%
    reframe(
      raw_veg=c(mean(raw_veg_QT), sd(raw_veg_QT)),
      cooked_veg=c(mean(cooked_veg_QT), sd(cooked_veg_QT)),
      fresh_fruit=c(mean(fresh_fruit_QT), sd(fresh_fruit_QT)),
      dried_fruit=c(mean(dried_fruit_QT), sd(dried_fruit_QT)),
      coffee=c(mean(coffee_QT), sd(coffee_QT)),
      tea=c(mean(tea_QT), sd(tea_QT)),
      addsalt=c(mean(addsalt_freq_QT, na.rm=T), sd(addsalt_freq_QT, na.rm=T))
    ) %>% t() %>%
    as.data.frame() %>%
    mutate(ancestry=rep(ANC, nrow(.)),
           food=gsub(".*[.]", "", rownames(.)))
  names(tab)[1:2] <- c("mean", "sd")
  return(tab)
}))

msd_dietByanc %>%
  mutate(food=factor(food, levels=c("raw_veg", "cooked_veg", "fresh_fruit", "dried_fruit", "coffee", "tea", "addsalt"))) %>%
  mutate(ancestry=factor(ancestry, levels=c(ANCs))) %>%
  ggplot(aes(x=food, y=mean, ymin=mean, ymax=mean+sd, fill=food, color=food)) +
  facet_grid(~ancestry) +
  geom_bar(stat="identity", position = position_dodge(0.93)) + 
  geom_errorbar(width=0.35, position = position_dodge(0.93)) + 
  scale_fill_manual(values=c(palettes$greens[c(2,4)], palettes$oranges[c(2,4)], palettes$purples[c(2,4)], "brown")) +
  scale_color_manual(values=c(palettes$greens[c(2,4)], palettes$oranges[c(2,4)], palettes$purples[c(2,4)], "brown"))


# ==========================================================================================
## Table of bitter food intakes by common diplotypes within each ancestry group
# ==========================================================================================

tab_bfood_byAnc.l <- lapply(ANCs, function(ANC){
  analysis.l[[ANC]] %>%
    group_by(diplos_common) %>%
    reframe(
      addsalt=mean_sd(addsalt_freq_QT, d=2),
      coffee=mean_sd(coffee_QT, d=2),
      tea=mean_sd(tea_QT, d=2),
      raw_veg=mean_sd(raw_veg_QT, d=2)) %>% 
    as.data.frame() %>%
    pivot_longer(c(addsalt, coffee, tea, raw_veg)) %>%
    pivot_wider(values_from=value, names_from=diplos_common) %>%
    mutate(Ancestry = rep(ANC, nrow(.))) }
)

do.call(cbind.data.frame, lapply(ANCs, function(ANC) analysis.l[[ANC]] %>%
                                   reframe(
                                     addsalt=mean_sd(addsalt_freq_QT, d=1),
                                     coffee=mean_sd(coffee_QT, d=1),
                                     tea=mean_sd(tea_QT, d=1),
                                     raw_veg=mean_sd(raw_veg_QT, d=1)) %>%
                                   t()
))



