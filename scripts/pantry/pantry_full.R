# Basic functions


# load packages
library(tidyverse) ; library(table1)


#########################################
## ~~~~  Variable Lists & Labels  ~~~~ ##
#########################################

# ===================================
## Diet variables
# ===================================

dietVars <- c("TCALS", "CHO_pct", "FIBER", "CHO2FIB", "FAT_pct", "PRO_pct", "ALC", paste0("dietPC", 1:10))

diet_labels <- c(
  raw_veg_QT="Raw vegetables",
  cooked_veg_QT="Cooked vegetables",
  fresh_fruit_QT="Fresh fruit",
  dried_fruit_QT="Dried fruit",
  oily_fish_QT="Oily fish",
  nonoily_fish_QT="Non-oily fish",
  procmeat_QT="Processed meat",
  poultry_QT="Poultry",
  cheese_QT="Cheese",
  beef_QT="Beef",
  lamb_QT="Lamb",
  pork_QT="Pork",
  bread_type_white_vs_brown_or_whole_BIN="Prefer white bread",
  bread_intake_QT="Bread",
  milk_type_full_vs_low_or_nonfat_BIN="Prefer full-fat milk",
  cereal_type_sugar_vs_any_bran_BIN="Prefer sugary cereal",
  coffee_type_decaf_vs_regular_BIN="Prefer decaf coffee",
  cereal_intake_QT="Cereal",
  spread_type_butter_vs_any_other_BIN="Prefer butter>oil",
  coffee_QT="Coffee",
  tea_QT="Tea",
  water_QT="Water",
  addsalt_always_often_vs_nrs_BIN="Often add salt",
  hotdrink_temp_hot_or_vhot_vs_warm_BIN="Prefer hot drinks"
)

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

dietVars <- c("raw_veg_QT", "cooked_veg_QT", "fresh_fruit_QT", "dried_fruit_QT", "oily_fish_QT", "nonoily_fish_QT", 
              "procmeat_QT", "poultry_QT", "cheese_QT", "beef_QT", "lamb_QT",  "pork_QT", "bread_type_white_vs_brown_or_whole_BIN",
              "bread_intake_QT", "milk_type_full_vs_low_or_nonfat_BIN", "cereal_type_sugar_vs_any_bran_BIN", "cereal_intake_QT", 
              "spread_type_butter_vs_any_other_BIN", "coffee_type_decaf_vs_regular_BIN",  "coffee_QT", "tea_QT", "water_QT",
              "addsalt_always_often_vs_nrs_BIN",  "hotdrink_temp_hot_or_vhot_vs_warm_BIN")


################################################################################
## ~~~~~~~~ Customized functions to clean, analyse & visualize data  ~~~~~~~~ ##
################################################################################


##########################################
##  Data display for descriptive tables ##
##########################################

# ======================================
## Print continuous vars as mean +- SD 
# ======================================
mean_sd<-function(x, d=2) {
  sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
          round(sd(x, na.rm=T), digits = d))
}


# ==================================
## Print categorical vars as n (%)
# ==================================
n_pct <- function(x, level=F, d=1) {
  if(level==F) {
    sapply(as.list(names(table(x))), function(lvl) {
      paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100, d), "%)") }) } 
  else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100, d), "%)")}
}


# =====================================================
## Print continuous vars as median [25th, 75th %tile]
# =====================================================
median_25to75<-function(x, d=2) {
  qs<-round(quantile(x, breaks=seq(0,1,0.25), na.rm=T), d)
  sprintf("%s [%s, %s]", round(median(x, na.rm=T), digits = d), 
          qs[[2]], qs[[4]])
}



################################################
##  ~~~~ Data wrangling & transformation ~~~~ ##
################################################

# ========================
## Remove outliers by SD
# ========================
remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}

find_outliers.fun <- function(x, SDs=5, recode_outliers=1) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- as.factor(ifelse(x>bounds[1] & x<bounds[2], 0, recode_outliers))
  x
}

find_in_range.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- as.factor(ifelse(x>bounds[1] & x<bounds[2],1,0))
  x
}



# =============================
## Flag missing values as 1/0
# =============================

recode_na_as.fun <- function(x, recode_na=0) {
  x <- ifelse(is.na(x) == T, recode_na, 1)
  x
}

recode_na_as_factor.fun <- function(x, recode_na=0) {
  x <- as.factor(ifelse(is.na(x) == T, recode_na, 1))
  x
}

# ================================================
## Median impute for negative or missing values
# ================================================
median_imp.fun <- function(x) {
  x.new <- ifelse(x == -1 | x == -3 | x == -9 | is.na(x) == T, median(x, na.rm=T), x)
  return(x.new)
}


# ===================
## Calculate zscore
# ===================
zscore.fun <- function(x) {
  z<-((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  return(z)
}


# ========================
## Winsorize data by SD
# ========================
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


# =========================
## Add descriptive labels 
# =========================

#e.g., labs_vals = female_labs <- list("Female" = 0, "Male" = 1)
descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}


descr_label_ordered.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; temp <- factor(temp, levels=names(labs_vals)) 
  return(temp)
}



# ===================================================
## Function to convert genotypes to amino acids
# ===================================================
geno_to_aa <- c("CAT"="AVI", "CAC"="AVV", "CGT"="AAI", 
                "CGC"="AAV", "CAV"="AVI", "GGC"="PAV", 
                "GGT"="PAI", "GAT"="PVI", "GAC"="PVV"
)

geno_to_aa.fun <- function(ANC.df) {
  haplo_aa <- ANC.df %>% 
    mutate(across(c("haplo_0", "haplo_1"), ~factor(
      gsub("CAT", "AVI", gsub("CAC", "AVV", gsub("CGT", "AAI", gsub("CGC", "AAV", gsub("CAV", "AVI",gsub("GGC", "PAV", gsub("GGT", "PAI", gsub("GAT", "PVI", gsub("GAC", "PVV", haplo_0)))))))))))) %>% 
    select(id, haplo_0, haplo_1)
  names(haplo_aa) <- c("id", "haplo_0_aa", "haplo_1_aa")
  return(haplo_aa)
}



##########################################################
##  ~~~~~ Statistical analysis/summary functions ~~~~~  ##
##########################################################

# =======================================
## Print descriptive table 1 by strata
# =======================================
source("../scripts/pantry/print_summary_table_fun.R", echo=F)


# ==============================================
## Function to get estimated marginal means
# ==============================================
get_emm.fun <- function(exposure, outcome, covars, reference, label=F, label.outcome=F, data=analysis) {
  exp.dat <- data %>% select(exp=exposure)
  dat <- data %>% mutate(exp=exp.dat$exp) ; dat$exp <- relevel(as.factor(dat$exp), ref=reference)
  mod <- lm(formula(paste0(outcome, "~", "exp", "+", covars)), dat)
  emm <- as.data.frame(emmeans(mod, ~ exp, rg.limit=180000)) %>% mutate(outcome=outcome, .after=exp) %>% 
    mutate(n = as.vector(table(mod$model$exp)), .before=emmean)
  anv <- anova(mod) #$'Pr(>F)'[1]
  
  out<-matrix(NA, nrow(emm), 10, dimnames = list(NULL, c("outcome", "model", "exposure", "level", "n", "emmean", "SE", "df", "lowCI", "upCI")))
  out[,1]=rep(ifelse(label.outcome==F,outcome,label.outcome),nrow(out))
  out[,2]<-rep(ifelse(label==F,NA,label),nrow(out)) ; out[,3]<-rep(exposure,nrow(out))
  out[,4:10]<-as.matrix(emm[,c(1,3:8)]) ; out <- as.data.frame(out) %>% mutate(anv.p=anv[1,5])
  
  return(list(emm=as.data.frame(out), anv=anv.p))
}


# =================
## summarize lm
# =================
print_lm <- function(exposure="taste_diplos", outcome="glu", covariates=m, label, 
                     label.outcome=F, round=F, digits=c(3,3), lm_trend=T, data=analysis) {
  
  mod <- lm(formula(paste0(outcome, "~", exposure, "+", covariates)), data)
  if(label==F) {lab <- outcome} else {lab <- label}
  if(label.outcome==F) {lab.outcome <- outcome} else{lab.outcome <- label.outcome}
  
  # For categorical exposure variable 
  if(is.numeric(data[[exposure]]) == F) {
    nX <- length(mod$xlevels[[exposure]])
    out<-matrix(NA, nX, 9, dimnames = list(NULL, c("outcome", "model", "exposure", "n", "beta", "se", "p", "f", "f_p")))
    out[,1] <- rep(lab.outcome, nrow(out)) ; out[,2] <- rep(lab, nrow(out)) ; out[,3] <- mod$xlevels[[exposure]]
    out[2:nrow(out),c(5:7)] <- summary(mod)$coef[2:nX, c(1:2,4)] ; out[1,8:9] <-c(anova(mod)[1,4], anova(mod)[1,5])
    out[,4] <- as.vector(table(mod$model[[exposure]]))
    # If lm_trend should be printed
    if(lm_trend == T) {
      exposure_num <- as.numeric(data[[exposure]])
      data <- data %>% mutate(exposure.num=exposure_num)
      lm.trend <- as.vector(summary(lm(formula(paste0(outcome, "~exposure.num+", covariates)), data))$coef[2,])
      out <- as.data.frame(out) %>% mutate(trend_p=c(lm.trend[4], rep(" ", nX-1) )) } 
  } 
  
  # For continuous exposure variable
  if(is.numeric(data[[exposure]]) == T) {
    out<-matrix(NA, 1, 7, dimnames = list(NULL, c("outcome", "model", "exposure", "n",  "beta", "se", "p")))
    out[5:7] <- summary(mod)$coef[2,c(1:2,4)]
    out[,4] <-length(mod$fitted.values) ; out[,c(1:3)] <- c(lab.outcome, lab, exposure) 
    out <- as.data.frame(out) %>% mutate(across(c("beta", "se", "p"), ~as.numeric(.)))
  }
  
  # Convert values to numeric
  out <- out %>% mutate(across(c("beta", "se"), ~as.numeric(.)))
  
  # If values should be rounded
  if(round == T) {
    out <- data.frame(out) %>%
      mutate(across(c("beta", "se"), ~round(as.numeric(.), digits[1]))) %>%
      mutate_at("p", ~round(as.numeric(.), digits[2]) ) 
  }
  return(as.data.frame(out))
}


print_glm <- function(exposure="taste_diplos", outcome, covariates=m, label, data=analysis) {
  
  mod <- glm(formula(paste0(outcome, "~", exposure, "+", covariates)), family=binomial("logit"), data) 
  
  # For categorical exposure variable 
  if(is.numeric(data[, ..exposure][[1]]) == F) {
    nX <- length(mod$xlevels[[exposure]])
    out<-matrix(NA, nX, 6, dimnames = list(paste0(label, "_", mod$xlevels[[exposure]]), c("n", "beta", "se", "p", "f", "f_p")))
    out[2:nrow(out),c(2:4)] <- summary(mod)$coef[2:nX, c(1:2,4)] ; out[1,5:6] <-c(anova(mod)[1,4], anova(mod)[1,5])
    out[,1] <- as.vector(table(mod$model[[exposure]])) 
  }
  
  # For continuous exposure variable
  if(is.numeric(data[, ..exposure][[1]]) == T) {
    out<-matrix(NA, 1, 5, dimnames = list(exposure, c("n", "Model", "beta", "se", "p")))
    out[3:5] <- summary(mod)$coef[2,c(1:2,4)]
    out[,1] <-length(mod$fitted.values) ; out[,2] <- NA ; out[,2] <- label
    out <- as.data.frame(out) %>% mutate(across(c("beta", "se", "p"), ~as.numeric(.)))
  }
  
  return(as.data.frame(out))
  
}


# ==================================
## summarize lm with interaction
# ==================================
print_lm_interaction <- function(exposure="taste_diplos", interaction, outcome="glu", 
                                 model=m, label_model, label_interaction, print_interaction_term, 
                                 lm_trend=T, round=T, digits=c(3,6), data) {
  
  mod <- lm(formula(paste0(outcome, "~", exposure, "*", interaction, "+", model)), data)
  mod.anova <- anova(mod)
  if(label_model==F) {model_label <- outcome} ; if(label_interaction==F) {
    label_interaction <- paste0(exposure, "x", outcome) }
  
  # For categorical exposure variable 
  if(is.numeric(data[, ..exposure][[1]]) == F) {
    nExp <- length(mod$xlevels[[exposure]]) ; nInt <- length(mod$xlevels[[interaction]])
    out<-matrix(NA, nExp+nInt, 6, dimnames = list(c(mod$xlevels[[exposure]], mod$xlevels[[interaction]]), c("n", "beta", "se", "p", "f", "f_p")))
    out[2:nExp,c(2:4)] <- summary(mod)$coef[2:nExp, c(1:2,4)] ; out[1,5:6] <- c(anova(mod)[1,4], anova(mod)[1,5])
    out[(nExp+2):(nExp+nInt),c(2:4)] <- summary(mod)$coef[3:(3+nInt-2), c(1:2,4)] ; out[(nExp+1),5:6] <- c(mod.anova[2,4], mod.anova[2,5])
    out[,1] <- c(as.vector(table(mod$model[[exposure]])), as.vector(table(mod$model[[interaction]])))
    
    out <- rbind(c(rep(NA, 4), mod.anova[(nrow(mod.anova)-1),4], mod.anova[(nrow(mod.anova)-1),5]), out)
    rownames(out)[1] <- label_interaction
   
     # If lm_trend should be printed
    #if(lm_trend == T) {
     # exposure_num <- as.numeric(data[, ..exposure][[1]])
     # data <- data %>% mutate(exposure.num=exposure_num)
     # lm.trend <- as.vector(summary(lm(formula(paste0(outcome, "~exposure.num+", covariates)), data))$coef[2,])
    # out <- as.data.frame(out) %>% mutate(t_p=c(lm.trend[4], rep(" ", nX-1) )) } 
  } ; 

  # If values should be rouded
  if(round == T) {
    out <- data.frame(out) %>%
      mutate(across(c("beta", "se"), ~round(as.numeric(.), digits[1]))) %>%
      mutate_at("p", ~round(as.numeric(.), digits[2]) ) 
  }
  
  return(as.data.frame(out))
}


# ===================================
## Make print_lm tables "pretty"
# ===================================

make_pretty_lm <- function(lm_table, digits=c(3,3),show_SE=F, scientific_P=T) {
  
  # Set parameters
  d_est<-digits[1];d_pval<-digits[2]
  se<- as.numeric(lm_table$se)
  
  # Format P-value function
  format_p.fun <- function(p) {
    P<-as.numeric(ifelse(p <0.001, format(p, scientific = T, digits=d_pval), 
           round(p, digits=d_pval))) 
    return(P)
  }
  
  #Rename P-values if only p detected
  pretty_table <- lm_table %>%
    mutate(across(ends_with("_p"), ~as.numeric(., d_pval))) %>%
    mutate(across(ends_with("_p"), ~format_p.fun(.))) %>%
    mutate(beta_se = sprintf("%s (%s, %s)", round(beta, d_est), round(beta-1.96*se, d_est), round(beta+1.96*se, d_est)))
  if("f" %in% names(pretty_table) & nrow(pretty_table) > 1) {
    pretty_table <- pretty_table %>% 
      rename(Outcome=outcome, Model=model, Exposure=exposure) %>% 
      mutate_at("p", ~format(., scientific=T, digits=d_pval)) %>%
      rename(P_t.test=p) %>%
      rename_with(~paste0("P_", gsub("[_].*", "", .), ".test"), ends_with("_p"))
  } else {
    pretty_table <- pretty_table %>% 
      mutate(Outcome=outcome, Model=model, Exposure=exposure) %>% 
      mutate_at("p", ~format_p.fun(.)) %>%
      rename(P_t.test=p) %>% mutate(N=n)
  } ; pretty_table <- pretty_table %>%
    select(Outcome, Model, Exposure, N=n, Beta_95CI=beta_se, starts_with("P")) %>%
    mutate_all(~ifelse(is.na(.), "-", .)) %>%
    mutate_at("Beta_95CI", ~ifelse(. == "NA (NA, NA)", "-", .))
  
  if(show_SE==T) { pretty_table <- pretty_table %>% mutate(SE=ifelse(is.na(se), "-", round(se, d_pval)), .after=Beta_95CI) } else{
    pretty_table <- pretty_table
  } ; rownames(pretty_table) <- NULL
  
  return(pretty_table)
}
  


######################################
##  ~~~~~ Plotting functions ~~~~~  ##
######################################

## ==============================
## Plot diplotypes as barplot
## ==============================

plot_Diplo.fun <- function(ANC) {
  
  anc <- analysis.l[[ANC]]
  
  # gather diplotypes with >0.05% frequency
  which(prop.table(table(anc$diplos_pal))*100>0.05)
  
  # Add yscale
  yscale <- max(prop.table(table(anc$diplos_pal))*100)*1.15
  
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
    scale_fill_manual(values = palettes$taster_labs, name="TAS2R38\nDiplotypes") + 
    theme(axis.text.x = element_text(size=8, angle=25, hjust=1, color = "black"),
          axis.text.y = element_text(size=8, color = "black"),
          axis.title = element_text(color = "black", face = "bold")) +
    geom_text(aes(label=n), vjust=-0.15, size=3) + 
    ggtitle(paste0(ANC, " (N = ", nrow(analysis_rda.l[[ANC]]), ")"))
}



# =========================
## Plot dietPC waterfall
# =========================
plot_dietPC_waterfall.fun <- function(dPC.df, nPCs=10) {
  dietPC.rot %>%
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
    scale_fill_manual(values=palettes$nature_main4, name = "Factor Loadings") +
    ggtheme_waterfall
}



###############################################
##  ~~~~  Color palettes & plot themes  ~~~~ ##
###############################################

# ============================================
## pre-built color palettes & ggplot themes
# ============================================

ggtheme_basic <- theme(
  axis.text = element_text(color="black", size=10), 
  axis.title.y = element_text(color="black", size=12, vjust=-0.5),
  axis.title.x = element_text(color="black", size=12, vjust=-0.5),
        legend.position = "bottom", 
        legend.box.background = element_rect(color = "black"),
        legend.key.size = unit(0.5, 'line'),
        legend.margin = margin(0.5,0.5,0.5,0.5, unit="pt"),
        legend.text = element_text(size=8), 
        legend.title = element_text(face="bold", size=8),
        plot.title=element_text(size=10),
        strip.text = element_text(face="bold", size=10))

# set ggplot theme
ggtheme_fancy <- theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.text.x = element_text(color="black", size=8), 
        axis.text.y = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=10),
        legend.position = "right", 
        legend.box.background = element_rect(color = "black"),
       # legend.key.size = unit(0.35, 'line'),
        legend.key.size = unit(0.5, 'line'),
        legend.margin = margin(0.5,0.5,0.5,0.5, unit="pt"),
        legend.text = element_text(size=8), 
        legend.title = element_text(face="bold", size=8),
        plot.title=element_text(size=8),
        strip.text = element_text(face="bold", size=8)
  )

ggtheme_figures <- theme(axis.text.x = element_text(color="black", size=8), 
        axis.text.y = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=9),
        panel.background = element_rect(color="#00000010"),
        legend.position = "right", 
        #legend.box.background = element_rect(color = "black"),
        #legend.key.size = unit(0.5, 'line'),
        legend.margin = margin(0.5,0.5,0.5,0.5, unit="pt"),
        legend.text = element_text(size=8), 
        legend.title = element_text(face="bold", size=8),
        plot.title=element_text(size=8),
        strip.text = element_text(face="bold", size=8)
  )

ggtheme_waterfall <- 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "top",
        strip.text = element_text(face="bold", size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x = element_text(size=6),
        legend.key.size = unit(0.25,"cm"))


ggtheme_emm <- theme(axis.title.y = element_text(size=11),
      axis.title.x = element_text(size=10, vjust=0.5),
      axis.text.y = element_text(size=11), axis.text.x = element_text(size=10),
      panel.background = element_rect(fill="#00000008"),
      panel.border = element_rect(color="black", fill=NA),
      legend.position = "right", legend.text = element_text(size=9),
      legend.title = element_text(face="bold", size=10))

# =================================
## Build pre-set color palettes 
# =================================
library("paletteer") ; library("RColorBrewer")
palettes <- list(NatComms= paletteer_d("ggsci::nrc_npg", n=10),
                 greens5 = paletteer_dynamic("cartography::green.pal", 5),
                 greens = rev(paletteer_dynamic("cartography::green.pal", 10)),
                 oranges = rev(paletteer_dynamic("cartography::orange.pal", 10)),
                 blues = rev(paletteer_dynamic("cartography::blue.pal", 10)),
                 purples=rev(paletteer_dynamic("cartography::purple.pal", 10))[3:10],
                 nature_main4=c("#888363","#C5C1A5", "#96A0B3", "#435269"),
                 blured_bin = c("#5496CE", "#DC6464"),
                 
                 #Categorical/assigned 
                 nature_haplos_dominant=c("AVI"="#C5C1A5", "AAV"="grey80", "PVI"="grey80", "PAV"="#96A0B3"),
                 nature_diplos_dominant=c("AVI/AVI"="#C5C1A5", "AVI/AAV"="grey80", "AVI/PAV"="#96A0B3", "PAV/PAV"="#435269", "AAV/PAV"="grey80"),
                 
                 # Diverging heatmaps
                 #hm_palette = colorRampPalette(c("#888363", "white", "#435269"))(n=100),
                 hm_palette_rwb = colorRampPalette(c("#00488D", "white", "#8F2F24"))(n=100),
                 hm_corr_blue = colorRampPalette(c("white", "#CAE5EF", "#96CED4", "#48BDBC", "#0096A0"))(n=100),
                 
                 #Binary codings
                 plaus_24hr=c("forestgreen", "grey38"), has_24hr=c("forestgreen", "grey38"),
                 
                 #3-level continuous
                 alleles_greyblue = c("#C6CAD6", "#96A0B3", "#6D788D"),
                 alleles_yellows = c("#F3DD92", "#E9C550", "#CB9B23"),
                 alleles_browns=c("#DDBCA1", "#BC9778", "#744F3D")
)



#A5A083

# ============================================
## ggplot themes
# ============================================

# set ggplot theme standards
ggplot_theme_standard_continuous <- theme_bw() + theme(
  axis.text.x = element_text(size=10, vjust=0.65, color = "black"),
  axis.text.y = element_text(size=10, color="black"), 
  strip.text=element_text(size=8, face="bold"),
  axis.title = element_text(size=10, color = "black"))

ggplot_theme_standard_categorical <- theme_bw() + theme(
  axis.text.x = element_text(size=10, color = "black", angle=30, hjust=0.9),
  axis.text.y = element_text(size=10, color="black"), 
  strip.text=element_text(size=8, face="bold"),
  axis.title = element_text(size=10, color = "black"))



# ============================================
## basic descriptive plots 
# ============================================

plot_continuous <- function(cont_var, data=analysis) {
  data %>% select(var=all_of(cont_var)) %>% filter(!is.na(var)) %>%
    ggplot(aes(x=var)) + geom_histogram(bins=30) +
    labs(title=cont_var, x=cont_var, y="frequency") + 
    theme(plot.title = element_text(face="bold", size=8)) #+
  #ggplot_theme_standard_continuous
}

plot_categorical <- function(cat_var, data=analysis) {
  data %>% 
    select(var=all_of(cat_var)) %>% filter(!is.na(var)) %>%
    ggplot(aes(x=factor(var))) + geom_bar(stat="count") +
    labs(title=cat_var, x=cat_var) + xlab("") +
    theme(axis.text.x = element_text(angle=35, hjust=0.75, size=7)) +
    theme(plot.title = element_text(face="bold", size=8))
  #ggplot_theme_standard_categorical
}

## EOF

