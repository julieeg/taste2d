# Title: Sensitivity analysis - ≠ fasting categories
# Date Updated: 08-27-2024





# ==========================================================================
### Explore what's underlying the 6-24 hr differences in RG by diplotype
# ==========================================================================
## LEFT OFF HERE ## 
# was goign to run all models in each time framt to comment on which time it was most seen in .. then do sex- or other stratified analyses to understand the 6+ phenom

# Look at more granular units of time
prop.table(table(analysis$fasting_hrs, analysis$sex),1) # 0-1; 2; 6; 7-8; >8
analysis <- analysis %>% mutate(
  fast_cat_v2 = case_when(
    fasting_hrs == 0 | fasting_hrs == 1 ~ "0to1hr",
    fasting_hrs == 2 ~ "2hr",
    fasting_hrs == 3 ~ "3hr",
    fasting_hrs == 4 ~ "4hr",
    fasting_hrs == 5 ~ "5hr",
    fasting_hrs == 6 ~ "6hr",
    fasting_hrs == 7 | fasting_hrs == 8 ~ "7to8hr",
    fasting_hrs > 8 & fasting_hrs <= 24 ~ "8+hr",
  )
)

fast_cat_v2.l <- c("0to1hr", "2hr", "3hr", "4hr", "5hr", "6hr", "7to8hr", "8+hr")



do.call(rbind.data.frame, lapply(fast_cat_v2.l, function(fast) {
  do.call(rbind.data.frame, lapply(1:length(models.nofast.l), function(i) {
    as.data.frame(print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.nofast.l[[i]], 
                           label=paste0(fast, "_", names(models.nofast.l)[i]), data=analysis %>% 
                             filter(fast_cat_v2 == fast)))  %>%
      mutate(model = rep(names(models.nofast.l)[[i]], 3), .before=beta)
  })) %>% mutate(fast=rep(fast, nrow(.)), model_fast=rownames(.), .before=beta) })) %>% 
  #Add P-values for linear trend
  left_join(
    do.call(rbind, lapply(fast_cat_v2.l, function(fast){
      do.call(rbind, lapply(models.nofast.l, function(m){
        coef(summary(lm(formula(paste0("glu~taste_diplos.num+", m)), data=analysis %>% filter(fast_cat_v2 == fast))))[2,3:4]
      })) %>% as.data.frame() %>% mutate(model_fast=paste0(fast, "_", names(models.nofast.l), "_AVI/AVI")) %>% rename("t" = `t value`, "t_p" = `Pr(>|t|)`) 
    })),
    by="model_fast") %>% write.csv(paste0(outDir_pf, "lm_diplo_x_gluXfast_sens_fastcatv2.csv"), row.names = F)
#save as .csv

