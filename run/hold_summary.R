


### Repeat analysis at visit 1 


# PAUSE: WILL NEED TO DOWNLOAD ALL OTHER COVARIATES FROM TIME1 TOO.....
sensitivity <- sensitivity %>% 
  mutate(fasting_hrs_keep.1 = ifelse(fasting_hrs.1 >= 25, NA, fasting_hrs.1)) %>%
  mutate(
    fast_cat_v1 = case_when(
      fasting_hrs_keep.1 >= 0 & fasting_hrs_keep.1 <=2 ~ "0to2hr",
      fasting_hrs_keep.1 ==3 ~ "3hr", fasting_hrs_keep.1 == 4 ~ "4hr",
      fasting_hrs_keep.1 == 5 ~ "5hr", fasting_hrs_keep.1 >=6 & fasting_hrs_keep.1 <=24 ~ "6+hr",
      TRUE ~ as.character(NA))
  ) ; sensitivity %>% ggplot(aes(x=as.factor(fasting_hrs), y=gluB_keep.0, fill=n_24hr.lab)) +
  geom_boxplot(position=position_dodge(0.9)) + scale_fill_manual(values=palettes$has_24hr)
sensitivity %>% ggplot(aes(x=as.factor(fast_cat_v1), y=gluB_keep.0, fill=n_24hr.lab)) +
  geom_boxplot(position=position_dodge(0.9)) + scale_fill_manual(values=palettes$has_24hr)

models.v1.l <- gsub("fast_cat", "fast_cat_v1", c(models.l[1:4]))
models.v1.l <- c(models.v1.l, paste0(m4_diet, "+dilution_keep.1+gluB_aliquot_keep.1.lab"))
names(models.v1.l) <- names(models.l)

## TAS2R38 diplotype & RG in diet patterns model ------------------------
sensitivity_diplo_v1_rg <- rbind(
  print_lm(exposure = "taste_diplos", outcome = "gluB_keep.1", label=names(models.v1.l[4]), 
           covariates = models.v1.l$Diet.Patterns, lm_trend = T, data = sensitivity),
  print_lm(exposure = "taste_diplos", outcome = "gluB_keep.1", covariates = models.v1.l$Assay.Performance, 
           lm_trend = T, label=names(models.v1.l)[5], data=sensitivity) )
make_pretty_lm(sensitivity_diplo_v1_rg) %>% kable(caption="TAS2R38 diplotype & RG with assay perf adjustment")

## TAS2R38 diplotype & 2hrGlu  ------------------------
sensitivity_assay_diplo_2hg <- rbind(
  print_lm(exposure = "taste_diplos", outcome = "glu", label=names(models.nofast.l[4]), 
           covariates = models.nofast.l$Diet.Patterns, lm_trend = T, data = sensitivity %>%
             filter(fast_cat == "0to2hr")),
  print_lm(exposure = "taste_diplos", outcome = "glu", covariates = models.nofast.l$Assay.Performance, 
           lm_trend = T, label=names(models.nofast.l)[5], data=sensitivity %>% 
             filter(fast_cat == "0to2hr")) )
sensitivity_assay_diplo_2hg %>% kable(caption="TAS2R38 diplotype & 2hr glucose with assay perf adjustment")

## TAS2R38 diplotype & HbA1c in diet patterns model ------------------------
sensitivity_assay_diplo_a1c <- rbind(
  print_lm(exposure = "taste_diplos", outcome = "hba1c", label=names(models.l[4]), 
           covariates = models.l$Diet.Patterns, lm_trend = T, data = sensitivity),
  print_lm(exposure = "taste_diplos", outcome = "hba1c", covariates = models.l$Assay.Performance, 
           lm_trend = T, label=names(models.l)[5], data=sensitivity) )
sensitivity_assay_diplo_a1c %>% kable(caption="TAS2R38 diplotype & HbA1c glucose with assay perf adjustment")
```





###### View markdown file for results tabulated in the whole cohort, before removing outliers


```{r, include=F}
### Q: How does GLUCOSE differ based on assay performance?
glucose %>% select(id, visit.0=gluB_keep.0, visit.1=gluB_keep.1) %>%
  pivot_longer(-id) %>%
  ggplot(aes(x=value)) + 
  geom_histogram() + 
  facet_wrap(~name, scales="free_y", nrow=1) + 
  theme(legend.position = "none")


# plotting function -----
plot_glu_by_assay <- function(var, Label, y_correction=4) {
  glucose %>%
    select(starts_with("gluB."), contains(var) & ends_with("lab"), "n_24hr.lab") %>%
    rename_all(., ~gsub("[.]lab", "", gsub(paste0("gluB_", var), "visit", .))) %>%
    pivot_longer(c("gluB.0", "gluB.1"), values_to="glu", names_to="gluB.visit") %>% 
    pivot_longer(c("visit.0", "visit.1"), names_to="visit", values_to="var") %>%
    filter(complete.cases(glu)) %>% group_by(var, n_24hr, visit) %>%
    reframe(mean=mean(glu, na.rm=T), sd=sd(glu, na.rm=T)) %>%
    filter(complete.cases(var)) %>%
    ggplot(aes(x = var, y=mean-y_correction, ymin=mean-y_correction-sd, ymax=mean-y_correction+sd, 
               fill=n_24hr, group=n_24hr)) +
    facet_wrap(~visit, scales="free_y") + 
    geom_bar(stat="identity", position=position_dodge(0.9)) + 
    geom_errorbar(position=position_dodge(0.9), linewidth=0.45, width=0.35) +
    scale_y_continuous(breaks=seq(0,2.5,0.5), labels=seq(0,2.5,.5)+y_correction) +
    scale_fill_manual(values=palettes$has_24hr) + 
    ylab("Glucose, mmol/L") + xlab(Label) + 
    theme(legend.position = "none",
          axis.title.x = element_text(face="bold"), axis.text.x = element_text(angle=35, hjust=1))
} ; n_24hr.l <- unique(glucose$n_24hr.lab)


# aliquot -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_aliquot_keep.0.lab", outcome="gluB_keep.0", label=paste0("Aliquot, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by aliquot, 24HR & visit") ; plot_glu_by_assay("aliquot", "Aliquot")

# correction type -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_correction_type_keep.0.lab", outcome="gluB_keep.0", label=paste0("Correction Level, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by correction type, 24HR & visit") ; plot_glu_by_assay("correction_type", "Correction Level")

# correction reason -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_correction_reason_keep.0.lab", outcome="gluB_keep.0", label=paste0("Correction Reason, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by correction reason, 24HR & visit")
```



```{r, out.width="25%", include=F}

# sample with glucose data -----
n_sensitivity_glu <- c(n_gluB.0="Blood glucose at v0, n(%)", n_gluB.1="Blood glucose at v1, n(%)",
                       n_gluM.0="NMR glucose at v0, n(%)", n_gluM.1="NMR glucose at v1, n(%)") 
print_summary_table(data=glucose, var_strata = F, vars_to_summarize = n_sensitivity_glu,
                    factor_vars = c(names(n_sensitivity_glu)), p_print=F) %>% 
  kable(caption="Sensitivity analysis sample after basic exclusions")

# n (%) and mean/sd blood & NMR glucose by 24HR ------
glucose_vars <- c(gluB.0="Blood glucose at v0, mmol/L", gluB.1="Blood glucose at v1, mmol/L",
                  gluM.0="Blood glucose at v0, mmol/L",gluM.1="Blood glucose at v1, mmol/L") 
print_summary_table(data=glucose, vars_to_summarize = glucose_vars, digits = c(4,1,6),
                    var_strata = "n_24hr.lab", p_types = "descriptive") %>%
  kable(caption = "Description of blood assays by 24HR data")


### Q: Are there differences in sample dilution factor by 24HR?

# overall dilution -----
print_summary_table(vars_to_summarize = c(dilution_pct.0 = "Sample dilution at v0, %",
                                          dilution_pct.1 = "Sample dilution at v1, %"), digits = c(4,1,6),
                    var_strata = "n_24hr.lab", p_types = "descriptive", data=glucose) %>%
  kable(caption = "Description of blood assays by 24HR data") ; glucose %>% 
  select(n_24hr.lab, starts_with("dilution.")) %>%
  pivot_longer(-n_24hr.lab) %>%
  mutate_at("name", ~gsub("dilution.", "Visit ", .)) %>%
  ggplot(aes(x=n_24hr.lab, y=value, fill=n_24hr.lab)) +
  facet_wrap(~name, scales="free_y") + 
  geom_hline(yintercept = 1) +
  geom_boxplot() + xlab("24HR data") + ylab("Sample dilution fraction, %") +
  scale_fill_manual(values=palettes$has_24hr, name="24HR data") + 
  theme(legend.position = "bottom") + 
  ggtitle(paste0("Sample dilution by 24HR data")) ; print_lm(
    exposure = "n_24hr.lab", outcome = "dilution.0", covariates = "dilution.0", 
    label = "Sample dilution x 24HR data", data=glucose) %>% as.data.frame() %>%
  kable(caption = "Sex/age adjusted P-value")

# dilution x fasting -----
glucose %>% 
  group_by(fasting_hrs.0, n_24hr.lab) %>%
  reframe(mean=mean(dilution.0, na.rm=T),
          sd=sd(dilution.0, na.rm=T)) %>%
  ggplot(aes(x=fasting_hrs.0, y=mean-0.985, ymin=mean-0.985-sd, ymax=mean-0.985+sd, fill=n_24hr.lab)) +
  geom_hline(yintercept = 1-0.985) +
  scale_y_continuous(limits=c(0,0.025), breaks=seq(0,0.025,0.005), labels=seq(0,0.025,0.005)+0.985) +
  scale_x_continuous(breaks=seq(0,23,1), labels=seq(1,24,1)) +
  geom_bar(stat="identity", position=position_dodge(0.9)) + 
  geom_errorbar(position=position_dodge(0.9), linewidth=0.45, width=0.35) +
  xlab("Fasting time at visit 0") +
  ylab("Sample Dilution Fraction, %") +
  scale_fill_manual(values=palettes$has_24hr, name="24HR data") + 
  theme(legend.position = "bottom") + 
  ggtitle("Sample dilution (fraction) by fasting time & 24HR data")


### Q: Are there differences in assay performance? 

# aliquot & correction type -----
print_summary_table(data=glucose, vars_to_summarize = c(
  gluB_aliquot.0.lab = "Blood glucose at v0, aliquot",
  gluB_aliquot.1.lab = "Blood glucose at v1, aliquot",
  gluB_correction_type.0.lab = "Blood glucose at v0, correction type",
  gluB_correction_type.1.lab = "Blood glucose at v1, correction type"), 
  digits = c(4,1,6), var_strata = "n_24hr.lab", p_types = "descriptive") %>%
  kable(caption = "Blood assay performance by 24HR data") ; ggarrange(
    plot_assay_perf_by_24hr(var="aliquot.", Label="Aliquot"),
    plot_assay_perf_by_24hr(var="correction_type.", Label="Correction Type"), 
    ncol=2, nrow = 1, common.legend = T)

#correction reason -------
print_summary_table(data=glucose, vars_to_summarize = c(
  gluB_correction_reason.0.lab = "Blood glucose at v0, correction reason",
  gluB_correction_reason.1.lab = "Blood glucose at v1, correction reason"), 
  digits = c(4,1,6), var_strata = "n_24hr.lab", p_types = "descriptive") %>%
  kable(caption = "Description of blood assays by 24HR data")


### Q: How does GLUCOSE differ based on assay performance?

# glucose histograms -----
glucose %>% select(id, visit.0=gluB.0, visit.1=gluB.1) %>%
  pivot_longer(-id) %>%
  ggplot(aes(x=value)) + 
  geom_histogram() + 
  facet_wrap(~name, scales="free_y", nrow=1) + 
  theme(legend.position = "none")

# aliquot -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_aliquot.0.lab", outcome="gluB.0", label=paste0("Aliquot, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by aliquot, 24HR & visit") ; plot_glu_by_assay("aliquot", "Aliquot")

# correction type -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_correction_type.0.lab", outcome="gluB.0", label=paste0("Correction Level, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by correction type, 24HR & visit") ; plot_glu_by_assay("correction_type", "Correction Level")

# correction reason -----
do.call(rbind.data.frame, lapply(n_24hr.l, function(k) {
  print_lm(exposure="gluB_correction_reason.0.lab", outcome="gluB.0", label=paste0("Correction Reason, v0 - ",k),
           covariates="age+sex", data=glucose %>% filter(n_24hr.lab == k )) })) %>% 
  kable(caption="Glucose by correction reason, 24HR & visit")

```




```{r}

# HOLD

## dilution with both timepoints
p_dilute_24hr <- glucose %>% 
  select(n_24hr.lab, starts_with("dilution_keep.")) %>%
  pivot_longer(-n_24hr.lab) %>%
  mutate_at("name", ~gsub("dilution_keep.", "Visit ", .)) %>%
  filter(complete.cases(value)) %>%
  ggplot(aes(x=n_24hr.lab, y=value, fill=n_24hr.lab)) +
  facet_wrap(~name, scales="free_y") + 
  geom_hline(yintercept = 1) +
  geom_boxplot() + xlab("24HR data") + ylab("Sample dilution fraction, %") +
  scale_fill_manual(values=palettes$has_24hr, name="24HR data") + 
  theme(legend.position = "bottom") + 
  ggtitle(paste0("Sample dilution by 24HR data")) ; print_lm_pretty(print_lm(
    exposure = "n_24hr.lab", outcome = "dilution.0", covariates = "dilution.0", 
    label = "Sample dilution x 24HR data", data=glucose)) %>% as.data.frame() %>%
  kable(caption = "Sex/age adjusted P-value - without outliers (>5SD")


