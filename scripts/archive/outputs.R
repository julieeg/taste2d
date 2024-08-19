## Manuscript Tables & Figures *Formatted*


# command args
ANC="EUR"

analysis <- readRDS(paste0("../data/processed/ukb_analysis_",ANC,".rda")) %>%
  filter(N_contrl_taste_kcal_compl=="1") %>%
  mutate(taster_status = factor(taster_status, levels=c("Nontaster", "Taster","Supertaster"))) 

dim(analysis)


# set variable names, palettes
palettes <- list(
  taster_labs=c("Nontaster" = "#F98400", "Taster"= "#F4B02E", "Supertaster"= "#00A08A",  "Other" = "grey50"),
  taster_cols=c("#F98400", "#F4B02E", "#00A08A", "grey50"),
  NatComms=paletteer_d("ggsci::nrc_npg", n=10)
)

dietPCs <- c(paste0("dietPC",1:10))


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
        legend.key.size = unit(0.25, 'line'),
        legend.margin = margin(0.2,0.2,0.2,0.2, unit="pt"),
        legend.text = element_text(size=8), 
        legend.title = element_text(face="bold", size=8),
        plot.title=element_text(size=8),
        strip.text = element_text(face="bold", size=8)
)


##############
##  TABLES  ##
##############


# =============================================
## Diplotypes with amino acid codings
# =============================================

recode_diplo_to_aa <- c(
  "CAT/GGC"="AVI / PAV", "CAT/CAT"="AVI / AVI", "GGC/GGC"="PAV / PAV",
  "CAT/CGC"="AVI / AAV", "CGC/GGC"="AAV / PAV", "CGC/CGC"="AAV/AAV", 
  "CAT/GAT"="AVI/PVI", "GAT/GGC"="PVI/PAV"
)

#Filters: filter(incl_t2d==1 & incl_compl==1 & incl_kcal==1) %>%
plot_Diplo.fun("EUR") + 
  xlab("") +
  theme(legend.position = "right",
        legend.title = element_text(face="bold",size=10),
        legend.text = element_text(size=9),
        legend.background = element_rect(color="black", linewidth=0.25),
        legend.key.size = unit(0.4,"cm"),
        axis.text.x = element_text(size=8.5, color = "black")) + 
  ggtitle("") +
  scale_x_discrete(labels=recode_diplo_to_aa)

ggsave(paste0("../output/EUR_diplosAA_legend.pdf"), height=3.5, width=5.5)

# =============================================
## Table 1 - Descriptives by taster status
# =============================================



# ===================================================
## Table 2 - GLM of taster status & glycemic traits
# ===================================================



# ===================================================
## Table 2 - GLM of taster status & 0-2hr glucose
# ===================================================





###############
##  FIGURES  ##
###############

# =====================================================
## Mean DP scores by taster status
# =====================================================

plot_mdp <- analysis_rda %>% 
  filter(incl_t2d==1 & incl_taste==1 & incl_kcal==1) %>% 
  select(taster_status, c(paste0("dietPC",1:5))) %>%
  filter(complete.cases(taster_status)) %>%
  group_by(taster_status) %>%
  reframe(dietPC1_mean = mean(dietPC1, na.rm=T),
          dietPC2_mean = mean(dietPC2, na.rm=T),
          dietPC3_mean = mean(dietPC3, na.rm=T),
          dietPC4_mean = mean(dietPC4, na.rm=T),
          dietPC5_mean = mean(dietPC5, na.rm=T),
          dietPC1_se = sd(dietPC1, na.rm=T) / sqrt(n()),
          dietPC2_se = sd(dietPC2, na.rm=T) / sqrt(n()),
          dietPC3_se = sd(dietPC3, na.rm=T) / sqrt(n()),
          dietPC4_se = sd(dietPC4, na.rm=T) / sqrt(n()),
          dietPC5_se = sd(dietPC5, na.rm=T) / sqrt(n())) %>%
  pivot_longer(cols = !taster_status, names_sep="_", names_to=c("diet", "msr")) %>%
  pivot_wider(names_from = msr, values_from = value) %>% 
  mutate(Diet=as.factor(gsub("dietPC", "Diet Pattern PC", diet))) %>%
  ggplot(aes(x = taster_status, y = mean, ymin=mean-1.96*se, ymax=mean+1.96*se, fill = taster_status)) +
  facet_grid(~Diet) + ggtheme +
  geom_hline(yintercept = c(seq(-0.2,0.1,0.05)), linewidth=0.15, color="grey") + 
  geom_bar(stat="identity", width=0.85) +
  geom_errorbar(width=0.15, linewidth=0.35) + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous(limits=c(-0.2,0.1), breaks=seq(-0.2,0.2,0.05)) +
  ylab("Mean (SE) PC scores") + xlab("") +
  scale_fill_manual(values=palettes$taster_labs) + 
  theme(axis.text.x = element_text(size=8, angle=25, hjust=1, color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"),
        legend.position = "none")
plot_mdp

ggsave(plot_mdp, filename="../output/EUR_plotb_dietPCs_meanse.pdf", height=3.25, width=8)
ggsave(plot_dietPC_waterfall.fun("EUR",5), filename="../output/EUR_plot_dietPCs_top5.pdf", height=4, width=8)



## Add pvalues
# function to calculate p-values across taste diplotypes
p_cont_bytaste.fun <- function(cont_var, data=analysis) {
  d <- data %>% select(taster_status,  y=all_of(cont_var)) #taster_num,
  anv.P <- c((anova(lm(y~taster_status, d)))$`Pr(>F)`[1],NA,NA)
  lm.P <- c(NA, summary(lm(y~taster_status, d))$coef[2:3,4])
 # lm.P.trend <- c(summary(lm(y~taster_num, d))$coef[2,4], NA, NA)
  d %>% group_by(taster_status) %>% 
    reframe(msd=mean_sd(y)) %>%
    mutate(lm.P=lm.P, 
           #lm.P.trend=lm.P.trend, 
           Anova.P=anv.P) %>%
    mutate(diet=rep(cont_var,nrow(.)), .before=msd) %>%
    return()
}

lapply(dietPCs, function(i) p_cont_bytaste.fun(i, data=analysis_rda %>% filter(incl_t2d==1 & incl_taste==1)))
lapply(dietPCs, function(i) p_cont_bytaste.fun(i, data=analysis))



# =====================================
## diet PCs for significnant oods
#=====================================

as.data.frame(dietPC_rot.l[[ANC]]$rotation) %>%
  mutate(Diet=diet_labels[gsub("_QT", "", gsub("_BIN", "", rownames(.)))]) %>%
  select(Diet, c(paste0("PC", 1:nPCs))) %>%
  mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
  pivot_longer(-Diet) %>% 
  mutate(name=factor(name, levels=c(paste0("PC", 1:10)) )) %>%
  mutate(name=gsub("PC", "Diet Pattern PC", name)) %>%
  mutate(direction = factor(ifelse(value>0.3, "Positive/Major", ifelse(value<0.3 & value>0, "Positive/Minor", 
                                                                       ifelse(value<0 & value>-0.3, "Negative/Minor", "Negative/Major"))),
                            levels=c("Positive/Major", "Positive/Minor", "Negative/Minor", "Negative/Major"))) %>%
  mutate(keep=ifelse(direction=="Negative/Major" | direction == "Positive/Major",1,0)) %>%
  filter(keep == 1) %>%
  ggplot(aes(x=value, y=Diet, fill=direction)) + 
  facet_wrap(~name, ncol=5) + 
  geom_col() + ylab("") + xlab("Rotated Factor Loadings") +
  geom_vline(xintercept = c(-0.3, 0.3), color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values=c(palettes$NatComms[c(1,4)]), name = "Factor Loadings") +
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




# =====================================================
## Beta (95% CI) RG by taster status / all models 
# =====================================================

tab_lm_rg %>%
  mutate(Model= gsub("[.]", " ", gsub("_.*"," Model", model))) %>%
  mutate(Model= factor(Model, levels=c("Base Model", "Lifestyle Model", "Carb Quality Model", "Diet Patterns Model"))) %>%
  mutate(taster_status=factor(gsub(".*_", "", rownames(.)), levels=c("Nontaster", "Taster", "Supertaster"))) %>%
  mutate(across(c(beta, se), ~ifelse(is.na(.), 0, .))) %>%
  ggplot(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se, x=Model, color=rev(taster_status))) +
  coord_flip() +
  ggtheme + 
  geom_hline(yintercept = 0, size=0.15) +
  geom_point(size=2, position=position_dodge(0.65)) + 
  geom_errorbar(width=0.15, position=position_dodge(0.65), linewidth=0.35) + 
  scale_y_continuous(limits=c(-0.04, 0.04)) +
  scale_color_manual(values=c(rev(palettes$taster_cols[1:3])), name="TAS2R38 Diplotypes") +
  ylab("Beta (95 %CI)") + xlab("")
  


# =============================================================
## Estimated marginal mean RG by taster status / all models 
# =============================================================

## Function to get estimated marginal means
get_emm.fun <- function(covars, outcome="glu", data=analysis) {
  library(emmeans)
  mod <- lm(formula(paste0(outcome, "~", "taster_status+", covars)), data)
  emm.glu <- emmeans(mod, ~ taster_status)
  emm.glu <- as.data.frame(emm.glu)
  anv.p <- anova(mod)$'Pr(>F)'[1]
  return(list(emm.glu=as.data.frame(emm.glu), anv.p=anv.p))
}


## Random Glucose in all models ---------------------------
emm_rg <- do.call(rbind.data.frame, lapply(1:4, function(m) {
  as.data.frame((get_emm.fun(models.l[[m]]))$emm.glu) %>% 
    mutate(model=rep(names(models.l)[[m]], nrow(.))) }))

emm_rg_p <- do.call(rbind.data.frame, lapply(1:4, function(m) {
  as.data.frame((get_emm.fun(models.l[[m]]))$anv.p) }))
p_labs = (tab_lm_rg %>% mutate(p=ifelse(is.na(p),"-",round(p,4))))$p


plot_emm_rg <- as.data.frame(emm_rg) %>% 
  mutate(model=factor(model, levels=c("Diet.Patterns", "Carb.Quality","Lifestyle", "Base"),
                      labels=c("Diet Patterns\nModel", "Carb Quality\nModel",
                               "Lifestyle\nModel", "Base\nModel") )) %>%
  mutate(taster_status = factor(taster_status, levels=c("Supertaster", "Taster", "Nontaster"))) %>%
  ggplot(aes(y=emmean, ymin=lower.CL, ymax=upper.CL, x=rev(model), 
             fill=taster_status, color=taster_status, shape=model)) +
  coord_flip() +
  scale_y_continuous(limits = c(4.89, 5.0), breaks=seq(4.9,4.96,0.02)) +
  scale_x_discrete(expand=expansion(add=c(0.5, 1.15))) +
  geom_hline(yintercept = seq(4.9,4.96,0.02), linewidth=0.15, color="grey") +
  ggtheme + 
  geom_point(size=2.5, position=position_dodge(0.65)) + 
  geom_errorbar(width=0.25, position=position_dodge(0.65), linewidth=0.55) + 
  ylab("Estimated marginal mean (95 %CI)") + xlab("") +
  scale_shape_manual(values=c(21,22,24,23)) +
  scale_fill_manual(values=c(palettes$taster_labs), name="TAS2R38 Diplotypes") + 
  scale_color_manual(values=c(palettes$taster_labs), name="TAS2R38 Diplotypes") +
  theme(legend.position = "none",
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(face="bold", size=8, hjust=0.325)) 

plot_emm_rg

ggsave(plot_emm_rg, filename="../output/EUR_plot_emmRG_v2.pdf", height=3, width=5)
     
             
## 0-2hr Glucose in the diet pattern model ---------------------------

fast_cat.l <- as.list(c("0to2hr", "3hr", "4hr", "5hr", "6+hr"))
emm_0to2hrg <- do.call(rbind.data.frame, lapply(fast_cat.l, function(f) {
  as.data.frame((get_emm.fun(models.nofast.l[[4]], data=analysis %>%
                               filter(fast_cat==f)))$emm.glu) %>% 
    mutate(fast=rep(f, nrow(.))) }))

#P-values
emm_2hrg_p <- do.call(rbind.data.frame, lapply(fast_cat.l, function(f) {
  as.data.frame((get_emm.fun(models.nofast.l[[4]], data=analysis %>% 
                               filter(fast_cat==f)))$anv.p) }))

tab_lm_0to2g %>% filter(model=="Diet.Patterns") %>% select(model_fast, p)

## Dot plot of EMM 
plot_emm_2hrg <- as.data.frame(emm_0to2hrg) %>% 
  mutate(fast=factor(fast, levels=c("6+hr", "5hr", "4hr", "3hr", "0to2hr"),
         labels=c("Fasting\n6+hr", "Fasting\n5hr", "Fasting\n4hr", "Fasting\n3hr", "Fasting\n0-2hr"))) %>%
  mutate(taster_status = factor(taster_status, levels=c("Supertaster", "Taster", "Nontaster"))) %>%
  ggplot(aes(y=emmean, ymin=lower.CL, ymax=upper.CL, x=fast, fill=taster_status, color=taster_status)) +
  coord_flip() + ggtheme + 
  scale_y_continuous(limits = c(4.84, 5.095), breaks=seq(4.85,5.0,0.05)) +
  scale_x_discrete(expand=expansion(add=c(0.5, 1.15))) +
  geom_hline(yintercept = seq(4.8,5.02,0.05), linewidth=0.15, color="grey") +
  geom_point(size=2, position=position_dodge(0.65), shape=23) + 
  geom_errorbar(width=0.25, position=position_dodge(0.65), linewidth=0.55) + 
  ylab("Estimated marginal mean (95 %CI)") + xlab("") +
  scale_fill_manual(values=c(palettes$taster_labs), name="TAS2R38 Diplotypes") + 
  scale_color_manual(values=c(palettes$taster_labs), name="TAS2R38 Diplotypes") + 
  theme(legend.position = "none",
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(face="bold", size=8, hjust=0.325)) 

plot_emm_2hrg

ggsave(plot_emm_2hrg, filename="../output/EUR_plot_emm2hrg_v2.pdf", height=2.95, width=5)

