# packages
library(tidyverse) ; library(data.table) ; library(paletteer)

palettes <- list(NatComms=paletteer_d("ggsci::nrc_npg", n=10))

# system arguments
commandArgs(trailingOnly = T)
ANC=args[1]
ANCs <- c("ALL", ANC)

genotypes <- c("rs713598_0", "rs713598_1", "rs1726866_0", "rs1726866_1", "rs10246939_0", "rs10246939_1", "haplo_0", "haplo_1", "diplo", "taster_status")
#genotypes <- c("rs713598_G", "rs1726866_G", "rs10246939_C", "haplo_0", "haplo_1", "diplo", "taster_status")


# Load ancestry-specific analysis file & subset to genotypes haplotypes by ancestry
df <- lapply(ANCs, function(x) {readRDS(paste0("../data/processed/ukb_analysis_", x, ".rda")) %>%
                select(id, all_of(genotypes)) %>% rename_with(~str_c(., "_", x), .cols=-"id")} ) %>%
  Reduce(function(df1, df2) right_join(df1, df2, by = "id"), .)

# Compare genotypes between ALL and ANC
match_sum <- function(var) {
  df %>% select(starts_with(var)) %>% rename_with(~gsub(".*ALL", "ALL", gsub(paste0(".*",ANC), "ANC", .))) %>%
    mutate(across(c("ALL", "ANC"), ~ifelse(is.na(.),  "xNA", .))) %>%
    mutate(match=ifelse(ANC == ALL, 1, 0)) %>% group_by(ANC) %>%
    summarise(Agree=sum(match==1, na.rm=T),
              Disagree=sum(match==0, na.rm=T))
}

# function to generate barplot of genotype matching
plot_match_sum <- function(var) {
  yscale <- max(match_sum(var)[,2:3])*1.10
  match_sum(var) %>% pivot_longer(-ANC) %>% 
    ggplot(aes(x=ANC, y=value, group=name, fill=name)) + ylim(0, yscale) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    geom_text(aes(label=value), position=position_dodge(0.8), size=3, angle=35) +
    scale_fill_manual(values=c(pals$NatComms[1:2]),
                      name=paste("Genotype match: ALL &", ANC)) +
    ylab("Frequency") + xlab("Genotype") +
    ggtitle(paste(var)) +
    theme(axis.text.x = element_text(color = "black", size=10, angle = 35)) +
    theme(axis.text.y = element_text(color = "black", size=10),
          legend.text = element_text(color = "black", size = 10),
          legend.title = element_text(color = "black", size = 12, face = "bold"),
          legend.position = "top") 
}

# save tables & plots
library(gridExtra) ; library(ggpubr)
pdf(paste0("../output/", ANC, "_genotype_comparison.pdf"), height = 8, width = 11)
par(mar=c(1.2, 2.1, 2.1, 1.2))
ggarrange(
  plot_match_sum(genotypes[1]), plot_match_sum(genotypes[2]),
  plot_match_sum(genotypes[3]), plot_match_sum(genotypes[4]),
  plot_match_sum(genotypes[5]), plot_match_sum(genotypes[6]),
  ncol=2, nrow=3, common.legend=T, legend = "top")
ggarrange(ggarrange(
  plot_match_sum(genotypes[7]), plot_match_sum(genotypes[8]), 
  ncol=2, common.legend=T, legend = "top"), plot_match_sum(genotypes[9]),
  nrow=2, ncol=1)
plot_match_sum(genotypes[10])
dev.off()


  plot_match_sum(genotypes[9])+theme(axis.text.x = element_text(angle=35))


length(genotypes)


# THEN MOVE ON TO BUILD OUT PHENOS_DESCR SYNTAX (CAN BE USED FOR REDI MR TOO) THEN SIMPLE ANALYSIS FRAMEWORK - CAN BE DONE TOFAY !!

# ASK TXEMA ABOUT ANCESTRY-SPECIFIC...once you have some more details on people with discrepant genotypes 

eur



haplos_id <- haplos %>%
  mutate(taster_status = case_when(
    haplo_0 == "GGC" & haplo_1 == "GGC" ~ "supertaster",
    haplo_0 == "GGC" & haplo_1 == "CAT" | haplo_1 == "GGC" & haplo_0 == "CAT" ~ "taster",
    haplo_0 == "CAT" & haplo_1 == "CAT" ~ "nontaster",
    TRUE ~ as.character(NA))) %>% 
  mutate(taster_status = ifelse(is.na(taster_status), "other", taster_status))


#-->look at combinationsof 2 not 3 snps
# use a positive control and comparison: qwhich phasing (if any) still uphold known associations
# tasters and other with salty outcome --> compare across


left_join(
  analysis %>% select(id, rs713598_G, rs1726866_G, rs10246939_C, haplo_0, haplo_1, diplo, taster_status), 
  all %>% select(id, rs713598_G, rs1726866_G, rs10246939_C, haplo_0, haplo_1, diplo, taster_status), by = "id", 
  suffix = c("_eur", "_all")) %>% 
  pivot_longer(-id, names_to="genotypes") %>%
  group_by(taster_status_eur) %>%
  summarise(Nontasters_EUR = sm()
            kable(caption = "Comparison of taste diplotypes from phasing on EUR vs ALL participants")
            
            ```