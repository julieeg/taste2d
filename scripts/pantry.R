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
descr_label.fun <- function(data, base_var, labs_vals) {
  
  base <- data %>% select(base_var) 
  temp <- rep(NA, length(base))
  
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}


# ===================================================
## Function to convert genotypes to amino acids
# ===================================================
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
source("../scripts/functions/print_summary_table_fun.R", echo=F)

# ==============================================
## Function to get estimated marginal means
# ==============================================
get_emm.fun <- function(exposure, outcome, covars, reference, data=analysis) {
  exp.dat <- data %>% select(exp=exposure)
  dat <- data %>% mutate(exp=exp.dat$exp) ; dat$exp <- relevel(as.factor(dat$exp), ref=reference)
  mod <- lm(formula(paste0(outcome, "~", "exp", "+", covars)), dat)
  emm <- as.data.frame(emmeans(mod, ~ exp, rg.limit=152000)) %>% mutate(outcome=outcome, .after=exp) %>% 
    mutate(n = as.vector(table(mod$model$exp)), .before=emmean)
  anv.p <- anova(mod) #$'Pr(>F)'[1]
  return(list(emm=as.data.frame(emm), anv.p=anv.p))
}


# =================
## summarize lm
# =================
print_lm <- function(exposure="taste_diplos", outcome="glu", covariates=m, label, data=analysis, interaction_term=F) {
  
  mod <- lm(formula(paste0(outcome, "~", exposure, "+", covariates)), data)
  
  # For categorical exposure variable 
  if(is.numeric(data[, ..exposure][[1]]) == F) {
    nX <- length(mod$xlevels[[exposure]])
    out<-matrix(NA, nX, 6, dimnames = list(paste0(label, "_", mod$xlevels[[exposure]]), c("n", "beta", "se", "p", "f", "f_p")))
    out[2:nrow(out),c(2:4)] <- summary(mod)$coef[2:nX, c(1:2,4)] ; out[1,5:6] <-c(anova(mod)[1,4], anova(mod)[1,5])
    out[,1] <- as.vector(table(mod$model[[exposure]]))
    
   # # If interaction_term == T 
  #  if(interaction_term != F) {
  #    nInt <-  length(mod$xlevels[[interaction_term]])
  ##    int <- summary(mod)$coef %>% as.data.frame() %>% filter(startsWith(rownames(.), exposure)==T)
  #    out.int<-matrix(NA, nX, 6, dimnames = list(paste0(label, "_", mod$xlevels[[exposure]]), c("n", "beta", "se", "p", "f", "f_p")))
  #    out[2:nrow(out),c(2:4)] <- [(nrow(mod.coefs)-nX):nrow(mod.coefs), c(1:2,4)] ; out[1,5:6] <-c(anova(mod)[1,4], anova(mod)[1,5])
    
  #}
    
  }
  
  # For continuous exposure variable
  if(is.numeric(data[, ..exposure][[1]]) == T) {
    out<-matrix(NA, 1, 4, dimnames = list(exposure, c("Model", "beta", "se", "p")))
    out[2:4] <- summary(mod)$coef[2,c(1:2,4)]
    out[,1] <- NA ; out[,1] <- label
  }
  
  
  # Return summary table
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
    out<-matrix(NA, 1, 4, dimnames = list(exposure, c("Model", "beta", "se", "p")))
    out[2:4] <- summary(mod)$coef[2,c(1:2,4)]
    out[,1] <- NA ; out[,1] <- label
  }
  
  return(as.data.frame(out))
  
}


######################################
##  ~~~~~ Plotting functions ~~~~~  ##
######################################

## ===========================
## Plot diplotypes as barplot
## ===========================

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

# set ggplot theme
ggtheme <- theme_bw() + 
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


library("paletteer") ; library("RColorBrewer")
palettes <- list(NatComms= paletteer_d("ggsci::nrc_npg", n=10),
                 greens5 = paletteer_dynamic("cartography::green.pal", 5),
                 greens = rev(paletteer_dynamic("cartography::green.pal", 10)),
                 pugr = c("#9C50CE", "#9C50CE95","#30900195", "#309001"),
                 #oranges5 = rev(paletteer_dynamic("cartography::orange.pal", 5)),
                 oranges = rev(paletteer_dynamic("cartography::orange.pal", 10)),
                 #blues5 = rev(paletteer_dynamic("cartography::blue.pal", 5)),
                 blues = rev(paletteer_dynamic("cartography::blue.pal", 10)),
                 purples=rev(paletteer_dynamic("cartography::purple.pal", 10))[3:10],
                 CovarGroups=c(brewer.pal(9,"Oranges")[4:6], brewer.pal(9, "Greens")[5:4], 
                               brewer.pal(9, "Blues")[5:6],  brewer.pal(9, "Purples")[c(8:4)], 
                               brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6]),
                 nature_main4=c("#888363","#C5C1A5", "#96A0B3", "#435269"),
                 nature_diplosAA_canonical=c("AVI/AVI"="#CB9B23", "AVI/AAV"="grey50", "AVI/PAV"="#E96900",
                                   "PAV/PAV"="#0096A0", "AAV/PAV"="grey50"),
                 nature_haplos=c("AVI"="#CB9B23","PAV"="#0096A0", "AAV"="grey50","PVI"="grey50"),
                 
                 nature_diplosAA_canonical_yrb=c("AVI/AVI"="#957628", "AVI/AAV"="grey50", "AVI/PAV"="#00488D",
                                             "PAV/PAV"="#8F2F24", "AAV/PAV"="grey50"),
                 nature_haplos_yrb=c("AVI"="#957628","PAV"="#8F2F24", "AAV"="grey50","PVI"="grey50"),
                 hm_palette = colorRampPalette(c("#888363", "white", "#435269"))(n=100)
)



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







# ================================= HOLD =============================== ##

########################################
##  Save tables/figures as csvs/pdfs  ##
########################################

#write wrapper function to generate latex summary output
#save_latex <- function(output, type, title, order, plwidth=0.9, show=F, ...) {
#  
#  if(type == "plot") {
#    ltx_plot(output, title=title, lwidth=paste0(plwidth,"\\linewidth"),
#             fontsize = 1, show=show, orientation = "portrait",
#             out = tempfile(paste0(ANC, "_", order), fileext = ".tex"), ...)
#  } ; if (type == "table") {
#    ltx_list(output, group=0, title=title, show=show,
#             footnote = "none", size=8, orientation = "portrait",
#             out=tempfile(paste0(ANC, "_", order), fileext = ".tex"))
#  }
#}

#outputs.l <- list(
#  out1 = list(output=tabl_sample, type="table", title=paste0("Sample Inclusions: ",ANC), order=1),
#  out2 = list(output=plot_descr_diplo, type="plot", title="Common TAS2R38 diplotypes", order=2, pheight=4, pwidth=8),
#  out3 = list(output=plots_descr_base_bm, type="plot", title="Demographic lifestyle and health characteristics by taster status", order=3, pheight=8, pwidth=6),
#  out4 = list(output=tabl_descr_gluA1c_fast, type="table", title="Descriptive statistics: Glucose and hbA1c by taster stauts", order=4),
#  out5 = list(output=plot_dietPC_waterfall, type="plot", title="Waterfall plot of food group factor loadings for diet PCs", pheight=7, pwidth=10, order=5),
#  out6 = list(output=tabl_descr_diet, type="table", title="Diet traits by taster status", order=6)
#)

#lapply(outputs.l, function(out) do.call(save_latex, out))
#ltx_combine(combine = tempdir(),  out=paste0("../output/", ANC, "_phenos_summary.tex"),
#            template=paste0(system.file(package="R3port"), "/default.tex"), orientation="portrait") 

# unlink tempdir to delete temporary files
#unlink(tempdir(), recursive = T)









