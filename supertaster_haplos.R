
# load packages
library(tidyverse) ; library(data.table) ; library(vcfR)

# system arguments
args = commandArgs(trailingOnly=TRUE)
i = args[1]

vcf <- read.vcfR(paste0("../data/processed/genos/supertasters/chunk", i, ".vcf"))

haps <- extract.haps(vcf, mask=F, unphased_as_NA=TRUE, verbose=T)
haps_id <- as.data.frame(t(cbind.data.frame(haps))) %>% mutate(tag=rownames(.)) %>% 
  mutate(parent=gsub(".*_", "", rownames(.)), id=gsub("_.*", "", rownames(.)))

left_join(
  haps_id %>% filter(parent==0) %>% rename_with(., ~paste0(., "_0"), starts_with("rs")) %>% select(-"parent"),
  haps_id %>% filter(parent==1) %>% rename_with(., ~paste0(., "_1"), starts_with("rs")) %>% select(-"parent"), by = "id") %>% 
  select("id", ends_with("0"), ends_with("1")) %>%
  
  mutate(haplo_0 = paste0(rs713598_0, rs1726866_0, rs10246939_0), 
         haplo_1 = paste0(rs713598_1, rs1726866_1, rs10246939_1)) %>%
  
  mutate(diplo = paste0(haplo_0, "/", haplo_1)) %>%
  
  mutate(taster_status = case_when(
    haplo_0 == "GGC" & haplo_1 == "GGC" ~ "supertaster",
    haplo_0 == "GGC" & haplo_1 == "CAT" | haplo_1 == "GGC" & haplo_0 == "CAT" ~ "taster",
    haplo_0 == "CAT" & haplo_1 == "CAT" ~ "nontaster",
    TRUE ~ as.character(NA))) %>% 
  mutate(taster_status = ifelse(is.na(taster_status), "other", taster_status)) %>%
      
  fwrite(paste0("../data/processed/genos/supertasters/chunk", i, "_haplos_id.csv"), row.names=F, col.names=T, quote=F)

#EOF


