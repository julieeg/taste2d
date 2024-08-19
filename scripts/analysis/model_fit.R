#model_building

# INCORPORATE FROM: taste_t2d_fib2cho_fast3hr_sens.Rmd // 
# LINK: file:///Users/jg1671/Documents/taste2d/Archive/workflow/run/taste_t2d_fib2cho_A.html



## Considerations on fasting covariates

#*Comparison of a) sample restricted to ≥3 hrs + continuous fasting hrs covariate vs. b) sample restricted to ≥2 and 3hours with categorical covariats*

get_modelfit.fun <- function(model, data=analysis) {
  mod <- model
  rmse <- sqrt(mean(mod$residuals^2))
  adj.R2 <- summary(mod)$adj.r.squared
  df<-summary(mod)$df[2]
  return(list(model_sum=summary(model), rmse=rmse, adj.R2=adj.R2, df=df))
} 


## Comparing alternative coding for Fasting Time variable

lm_fast_v1 <- lm(formula(paste0("glu~taste_diplos+", models.l$Diet.Patterns)), data=analysis)

get_modelfit.fun(lm_fast_v1)


## Fasting variables
bar_glu_x_fast <- analysis %>% 
  mutate(fasting_hrs_lt24 = ifelse(fasting_hrs >= 24, 24, fasting_hrs)) %>% 
  group_by(fasting_hrs_lt24) %>% 
  filter(complete.cases(fasting_hrs_lt24)) %>% 
  summarise(N=n(),
            mean=mean(glu, na.rm=T),
            sd=sd(glu, na.rm=T)) 

bar_glu_x_fast %>%
  ggplot(aes(x=fasting_hrs_lt24, y=mean-3.5, ymin=mean-3.5-sd, ymax=mean-3.5+sd)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_y_continuous(breaks=seq(0,3,0.5), limits = c(0,3), name="Blood glucose (mmol/L)", labels = seq(0,3,0.5)+3.5) +
  scale_x_continuous(limits=c(0,25), expand=c(0,0), breaks=seq(1,24.5,1)) +
  xlab("Time since last meal (hours)") +
  geom_errorbar(width=0.2) + 
  annotate("text", x=1:24, y=(bar_glu_x_fast$mean-3.5)+bar_glu_x_fast$sd+0.1, label=bar_glu_x_fast$N, size=2.5)




















## Lists of models, ALL restricting to ≥3 hr fasting
models.nofast.l <- list("Base"=base_covars, "BMI"=bmi_covars, "Lifestyle"=life_covars,
                        "Carb\nQuality"= chqIa_covars, "Gen.\nDiet"=diet_covars)
models.fastcont.l <- lapply(models.nofast.l[1:5], function(x) paste0(x, "+fasting_hrs"))
models.fastcat.l <- lapply(models.nofast.l[1:5], function(x) paste0(x, "+fast_cat")) 

## Among ≥3hr fasting
fit.fast.type.gt3 = as.data.frame(
  rbind(
    do.call(rbind.data.frame,lapply(models.nofast.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(exclude_fasting3hr==0)) } )),
    do.call(rbind.data.frame, lapply(models.fastcont.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(exclude_fasting3hr==0)) } )),
    do.call(rbind.data.frame, lapply(models.fastcat.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(exclude_fasting3hr==0)) } ))
  )) %>%
  mutate("type" = c(rep("A. No Fasting", 5), rep("B. Fasting (cont)", 5), rep("C. Fasting (cat)", 5)), .before=rmse) %>%
  mutate("model.order"=rep(LETTERS[1:5], 3), .before=rmse)

## Among ≥2hr fasting
fit.fast.type.gt2 = as.data.frame(
  rbind(
    do.call(rbind.data.frame,lapply(models.nofast.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt2hr==1)) } )),
    do.call(rbind.data.frame, lapply(models.fastcont.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt2hr==1)) } )),
    do.call(rbind.data.frame, lapply(models.fastcat.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt2hr==1)) } ))
  )) %>%
  mutate("type" = c(rep("A. No Fasting", 5), rep("B. Fasting (cont)", 5), rep("C. Fasting (cat)", 5)), .before=rmse) %>%
  mutate("model.order"=rep(LETTERS[1:5], 3), .before=rmse)

## Among ≥1hr fasting
fit.fast.type.gt1 = as.data.frame(
  rbind(
    do.call(rbind.data.frame,lapply(models.nofast.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt1hr==1)) } )),
    do.call(rbind.data.frame, lapply(models.fastcont.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt1hr==1)) } )),
    do.call(rbind.data.frame, lapply(models.fastcat.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(fasting_gt1hr==1)) } ))
  )) %>%
  mutate("type" = c(rep("A. No Fasting", 5), rep("B. Fasting (cont)", 5), rep("C. Fasting (cat)", 5)), .before=rmse) %>%
  mutate("model.order"=rep(LETTERS[1:5], 3), .before=rmse)

as.data.frame(rbind(fit.fast.type.gt1, fit.fast.type.gt2, fit.fast.type.gt3)) %>%
  mutate(fast=c(rep("≥1hr",15), rep("≥2hr", 15), rep("≥3hr", 15) )) %>% 
  ggplot(aes(x=model.order, y=adj.R2-0.02, group=type, fill=type)) + 
  facet_wrap(~fast, scales = "free") +
  geom_bar(stat = "identity", position=position_dodge(0.6), width = 0.5) +
  scale_x_discrete(labels = names(models.fastcat.l)) +
  xlab("") + ylab("Model Adjusted R2 (higher = better fit)")+
  scale_y_continuous(limits=c(0, 0.02), labels=seq(0,0.02,0.005)+0.02)




## Considerations on BMI inclusion 


get_modelfit.fun <- function(model_covars, data=ukb_sens) {
  mod <- lm(formula(paste0("glu~taster_status+", model_covars)), data=data)
  rmse <- sqrt(mean(mod$residuals^2))
  adj.R2 <- summary(mod)$adj.r.squared
  df<-summary(mod)$df[2]
  return(list(rmse=rmse, adj.R2=adj.R2, df=df))
} 

#!# Run using the ukb_sens, which is only CONTROLS but NOT restricted for fasting time --> add by hand.

## Main models, restricting to ≥3 hr fasting
models.l <- list("Base"=base_covars, "BMI"=bmi_covars, "Lifestyle"=life_covars,
                 "Carb\nQuality"= chqIa_covars, "Gen.\nDiet"=diet_covars, "Fasting"=fast_covars)

fit.fastgt3 = as.data.frame(
  rbind(
    do.call(rbind.data.frame,lapply(models.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(exclude_fasting3hr==0)) } )),
    do.call(rbind.data.frame, lapply(models.int.l, function(x) {
      get_modelfit.fun(x, data=ukb_sens %>% filter(exclude_fasting3hr==0)) } ))
  )) %>%
  mutate("type" = c(rep("Main Effect", 6), rep("With Interaction",6)), .before=rmse) %>%
  mutate("model.order"=rep(LETTERS[1:6],2), .before=rmse)

# plot of adjusted R2
fit.fastgt3 %>% 
  ggplot(aes(x=model.order, y=adj.R2)) + 
  facet_wrap(~type, scales = "free") +
  geom_bar(stat = "identity", position=position_dodge())+
  scale_x_discrete(labels = names(models.l)) +
  xlab("") + ylab("Model Adjusted R2 (higher = better fit") +
  scale_y_continuous(limits=c(0, 0.035))

# plot of RMSE
fit.fastgt3 %>% 
  ggplot(aes(x=model.order, y=rmse-0.53)) + 
  facet_wrap(~type, scales = "free") +
  geom_bar(stat = "identity", position=position_dodge())+
  scale_x_discrete(labels = names(models.l)) + xlab("") + ylab("Model RMSE (lower = better fit)") +
  scale_y_continuous(limits=c(0,0.01), breaks=seq(0, 0.01, 0.0025), labels=seq(0,0.01,0.0025)+0.53)

#kable(fit.main.fastgt3, caption = "Model fit statistics for main effects models")
#kable(fit.main.int.fastgt3, caption = "Model fit statistics for interaction effects models")
