## Post-processing

# load packages
lapply(c("tidyverse", "data.table", "paletteer", "ggpubr", "R3port", "xtable", "knitr", "tinytex"),  library, character.only = TRUE)
#tinytex::install_tinytex()

# load basic functions
source("../scripts/basic_functions.R", echo=F)


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


# Write wrapper functions to save outputs for LaTex


# build basic color palettes
palettes <- list(NatComms=paletteer_d("ggsci::nrc_npg", n=10),
                HeatMapYlRd=rev(paletteer::paletteer_c("grDevices::heat.colors", n=100)))

Taster_palette_Labs <- c("Nontaster" = palettes$NatComms[1], "Taster"= palettes$NatComms[4], 
                         "Supertaster"= palettes$NatComms[3], "Other" = "grey50")


# output path
out_prefix=paste0("../output/", ANC, "_")


###############################
##   Load data by ancestry   ##
###############################

# system arguments
#args = commandArgs(trailingOnly=TRUE)
ANC = "EUR"

#load data by ancestry
analysis_rda <- readRDS(paste0("../data/processed/ukb_analysis_", ANC, ".rda"))

#system(paste0("mkdir -p ", ))



# ===================================================================
### Exclusion criteria
# * No diabetes (Eastwood algorithm + HbA1c <5.7%)
# * Common TAS2R38 diplotype
# * [ANC] genetic ancestry
# * Complete data for glucose, genetic ancestry, covariates
# * Plausible total energy intakes (600-4000 kcal/d)
# ===================================================================

tabl_sample <- analysis_rda %>%
  select(id, starts_with("N_")) %>%
  pivot_longer(-id, names_to="criterion") %>%
  group_by(criterion) %>%
  summarise(Yes = as.integer(sum(value)),
            No = sum(!value)) %>%
  setNames(c("Inclusion criterion", "Yes", "No"))

analysis <- analysis_rda %>%
  filter(N_contrl_taste_kcal_compl==1)


# save as csv
#write.csv(tab_sample, paste0(out_prefix, "_tab_sample_N.csv"))

#paste0("../output/", ANC, "/tabl_sample.tex")
#ltx_list(tab_sample, vargroup = c("", "Included", "Included"), group=1,
 #        title="Sample Inclusion/Exclusions", footnote = "none", show=F,
 #        out=paste0(outputDir, "/tab1_sample.tex"), size=12, orientation="portrait")



################################
##   Descriptive statistics   ##
################################

plot_continuous <- function(cont_var) {
  analysis %>% select(var=all_of(cont_var)) %>% filter(!is.na(var)) %>%
    ggplot(aes(x=var)) + geom_histogram(bins=30) +
    labs(title=cont_var, x=cont_var, y="frequency") + 
    theme(plot.title = element_text(face="bold", size=8)) #+
  #ggplot_theme_standard_continuous
}

plot_categorical <- function(cat_var) {
  analysis %>% 
    select(var=all_of(cat_var)) %>% filter(!is.na(var)) %>%
    ggplot(aes(x=factor(var))) + geom_bar(stat="count") +
    labs(title=cat_var, x=cat_var) + xlab("") +
    theme(axis.text.x = element_text(angle=35, hjust=0.75, size=7)) +
    theme(plot.title = element_text(face="bold", size=8))
  #ggplot_theme_standard_categorical
}


# ====================================
## Basic phenotypes
# ====================================

plot_descr_basic <- ggarrange(
  ggarrange(plot_continuous("age"), plot_categorical("sex"), plot_continuous("bmi"), ncol=3, align="hv"),
  ggarrange(plot_categorical("smoke_level.lab"), plot_categorical("pa_met_excess_lvl"), 
            plot_categorical("educ_level.lab"), ncol=3, widths=c(1,1,1.33), align="hv"),
  ggarrange(plot_continuous("ALC"), plot_categorical("alch_freq.lab"), ncol=2, align="hv", 
            widths=c(1,1.5)),
  nrow=3, heights=c(1,1.75, 1.55))

plot_descr_bm <- ggarrange(
  ggarrange(plot_continuous("glu"), plot_continuous("hba1c"), nrow=1, align="hv"),
  ggarrange(plot_continuous("sbp"), plot_continuous("dbp"), plot_continuous("tg_log"),
  plot_continuous("ldl"), plot_continuous("hdl"), ncol=3, nrow=2, align="hv"), 
  nrow=2, heights = c(1,1.75))


plots_descr_base_bm <- list(plot_descr_basic, plot_descr_bm)


#save for latex
pdf(paste0(out_prefix, "_plot_descr_basic.pdf"))
plot_descr_basic
plot_descr_bm
dev.off()

#paste0("../output/", ANC, "/descr_basic.tex")
#ltx_plot(list(pl_descr_basic, pl_descr_bm), titlepr="Figure 1", orientation="portrait",
#         title = paste0("Descripive plots of Basic phenotypes for ", ANC, "participants"), 
#         lwidth="0.85\\linewidth", pheight = 7, fontsize = 1,
#         out=paste0(outputDir, "/p1_descr-basic.tex"), show=T)


# ======================================
# Taster diplotypes 
# ======================================

diplo_gt05=names(which(prop.table(table(
  (analysis_rda %>% filter(incl_t2d==1 & incl_kcal==1 & incl_compl==1))$diplo))>0.0005))

plot_descr_diplo <- analysis_rda %>%
  filter(diplo %in% diplo_gt05, incl_t2d==1, incl_kcal==1, incl_compl==1) %>%
  reframe(n=table(diplo), pct=prop.table(table(diplo)), diplo=names(table(diplo))) %>%
  arrange(-n) %>%
  mutate(
    color = ifelse(diplo == "CAT/CAT", "Nontaster",
                   ifelse(diplo %in% c("CAT/GGC", "GGC/CAT"), "Taster",
                          ifelse(diplo == "GGC/GGC", "Supertaster", "Other")))) %>%
  mutate(diplo_order=factor(diplo, levels=diplo[order(n, decreasing=T)])) %>%
  ggplot(aes(x=diplo_order, y=pct*100, fill=color)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(name = "\n TAS2R38 Diplotypes") +
  ylab("Frequency (%)") + #ylim(0, 0.5) +
  scale_fill_manual(values = Taster_palette_Labs, name="Taster Status") + 
  theme(axis.text.x = element_text(size=10, angle=35, vjust=0.65, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        axis.title = element_text(color = "black", face = "bold")) +
  geom_text(aes(label=n), vjust=-0.15, size=3) #+ 
  #ggtitle("TAS2R38 Diplotype Frequency Table") 

#paste0("../output/", ANC, "/descr_basic.tex")
#ltx_plot(pl_descr_diplo, title = paste0("Distribution of taster diplotypes"), 
 #        lwidth="0.90\\linewidth", pheight = 7, fontsize = 1, orientation="portrait",
  #        out=paste0(outputDir, "/p2_descr_diplo.tex"), show=T)


# ========================================
## Summary tables by taster status
# ========================================

## Demographic, behavioral & risk factors
tabl_descr_basic <- as.data.frame(analysis %>% 
  filter(complete.cases(taster_status)) %>%
  group_by(taster_status) %>%
  reframe(
    Age=mean_sd(age),
    Female=n_pct(sex, level = "Female"),
    BMI=mean_sd(bmi),
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
  t()) %>% 
  filter(row.names(.) != "taster_status") %>% 
  rename("Nontaster"=V1, "Taster"=V2, "Supertaster"=V3) %>%
  mutate(Variable=row.names(.), .before="Nontaster") 

# save for latex
#ltx_list(tabl_descr_basic, vargroup=c("", c(rep("Taster Status", 3))), group = 1, 
#         title="Participant characteritics by taster status", size=12,
#         fill="-",footnote="none", mancol="p{3cm}llll", show=F,
#         out=paste0(outputDir, "/tab2_descr_basic.tex"))


## Glucose/A1c by fasting time 
tabl_descr_gluA1c_fast <- bind_rows(
  analysis %>% filter(complete.cases(taster_status)) %>%
    group_by(fast_cat, taster_status) %>%
    reframe(Glucose=mean_sd(glu)) %>%
    pivot_wider(names_from="taster_status", values_from="Glucose") %>%
    mutate(Biomarker="Glucose", .before = fast_cat),
  analysis %>% filter(complete.cases(taster_status)) %>%
    group_by(fast_cat, taster_status) %>%
    reframe(HbA1c=mean_sd(hba1c_max)) %>%
    pivot_wider(names_from="taster_status", values_from="HbA1c") %>%
    mutate(Biomarker="HbA1c", .before = fast_cat)
) %>% tibble()


#save for latex
#ltx_list(tabl_descr_gluA1c_fast, vargroup=c("", "", c(rep("Taster Status", 3))), group = 1, 
 #        title="Participant characteritics by taster status", size=12,
  #       fill="-",footnote="none", mancol="p{3cm}llll", show=T,
   #      out=paste0(outputDir, "/tab3_descr_gluA1c_fast.tex"))

# ===================================
## Diet patterns 
# ===================================

dietPC.dat <- readRDS(paste0("../data/processed/diet/ukb_", ANC, "_dietPCs.rda"))
dietVars <- c("TCALS", "CHO_pct", "FIBER", "CHO2FIB", "FAT_pct", "PRO_pct", "ALC", paste0("dietPC", 1:10))

## waterfall plot of Factor loading
plot_dietPC_waterfall <- as.data.frame(dietPC.dat$rotation) %>%
  mutate(diet=gsub("_QT", "", gsub("_Bin", "", rownames(.)))) %>%
  select(diet, c(paste0("PC", 1:10))) %>%
  mutate(diet=factor(diet, levels=diet[order(PC1, decreasing=T)])) %>%
  pivot_longer(-diet) %>% 
  mutate(name=factor(name, levels=c(paste0("PC", 1:10)) )) %>%
  mutate(direction = ifelse(value <=0, "neg", "pos")) %>%
  ggplot(aes(x=value, y=diet, fill=direction)) + 
  facet_wrap(~name, ncol=5) + 
  #geom_vline(xintercept = 0, color = "black") +
  geom_col() + ylab("Diet traits derived from UKB FFQ") + xlab("Rotated Factor Loadings") +
  geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values=c(palettes$NatComms[7], palettes$NatComms[5]), 
                    name = "Factor Loading", labels = c("Negative", "Positive")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        legend.position = "top")


## Summary table by taster status
tabl_descr_diet <- as.data.frame(
  analysis %>% filter(complete.cases(taster_status)) %>% 
    group_by(taster_status) %>%
    reframe(
      Energy_kcal=mean_sd(TCALS),
      Carbohydrate_pct = mean_sd(CHO_pct),
      Fiber_g = mean_sd(FIBER),
      Carb_to_Fiber=mean_sd(CHO2FIB),
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
             data=analysis %>% filter(taster_status != "Other")))$`P`[1]} ))) %>%
  rename("Nontasters"=V1, "Tasters"=V2, "Supertasters"=V3) %>% 
  filter(rownames(.) != "taster_status") %>%
  mutate("Variable" = row.names(.), .before="Nontasters")


# save for latex
#ltx_plot(pl_dietPC_waterfall, out = paste0(outputDir, "/p3_dietPCs.tex"), 
 #        lwidth="0.85\\linewidth", pheight = 7, fontsize = 1, show=F)

#ltx_list(tabl_descr_diet, vargroup = c("", c(rep("Taster Status",3))), group=1,
 #        title="Diet characteritics by taster status", size=12,
  #       fill="-",footnote="none", show=T, 
   #      out=paste0(outputDir, "/tab4_descr_diet.tex"))



########################################
##  Save tables/figures as csvs/pdfs  ##
########################################

#write wrapper function to generate latex summary output
save_latex <- function(output, type, title, order, plwidth=0.9, show=F, ...) {
  
  if(type == "plot") {
    ltx_plot(output, title=title, lwidth=paste0(plwidth,"\\linewidth"),
             fontsize = 1, show=show, orientation = "portrait",
             out = tempfile(paste0(ANC, "_", order), fileext = ".tex"), ...)
  } ; if (type == "table") {
    ltx_list(output, group=0, title=title, show=show,
             footnote = "none", size=8, orientation = "portrait",
             out=tempfile(paste0(ANC, "_", order), fileext = ".tex"))
  }
}

outputs.l <- list(
  out1 = list(output=tabl_sample, type="table", title=paste0("Sample Inclusions: ",ANC), order=1),
  out2 = list(output=plot_descr_diplo, type="plot", title="Common TAS2R38 diplotypes", order=2, pheight=4, pwidth=8),
  out3 = list(output=plots_descr_base_bm, type="plot", title="Demographic lifestyle and health characteristics by taster status", order=3, pheight=8, pwidth=6),
  out4 = list(output=tabl_descr_gluA1c_fast, type="table", title="Descriptive statistics: Glucose and hbA1c by taster stauts", order=4),
  out5 = list(output=plot_dietPC_waterfall, type="plot", title="Waterfall plot of food group factor loadings for diet PCs", pheight=7, pwidth=10, order=5),
  out6 = list(output=tabl_descr_diet, type="table", title="Diet traits by taster status", order=6)
  )

lapply(outputs.l, function(out) do.call(save_latex, out))
ltx_combine(combine = tempdir(),  out=paste0("../output/", ANC, "_phenos_summary.tex"),
            template=paste0(system.file(package="R3port"), "/default.tex"), orientation="portrait") 

# unlink tempdir to delete temporary files
unlink(tempdir(), recursive = T)

#EOF



