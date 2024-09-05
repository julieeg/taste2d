






## Sex-stratified
sex.l <- c("Female", "Male")
models.nofast.nosex.l <- lapply(models.nofast.l, function(m) gsub("[+]sex", "", m))

do.call(rbind.data.frame,lapply(sex.l, function(mf){
  do.call(rbind.data.frame, lapply(1:length(models.nofast.nosex.l), function(m) { 
    do.call(rbind.data.frame, lapply(fast_cat.l, function(fast) {
      get_emm.fun(exposure = "taste_diplos", outcome = "glu", covars = models.nofast.nosex.l[[m]], reference = "AVI/AVI", 
                  data=analysis %>% filter(fast_cat == fast & sex == mf))$emm %>%
        mutate(fast_cat = fast, .after=n) })) %>%
      mutate(model = names(models.nofast.nosex.l)[[m]], .after=outcome) })) %>%
    mutate(sex = mf, .before=model)  })) %>%
  write.csv(paste0(outDir_pf, "emm_diplo_x_gluXfast_v2_allmodels_bysex.csv"), row.names = F)


