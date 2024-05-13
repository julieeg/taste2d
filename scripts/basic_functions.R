## Basic functions


# load packages
library(tidyverse) ; library(table1)


# Functions

## Print continuous vars as mean +- SD 
mean_sd<-function(x, d=2) {
  sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
          round(sd(x, na.rm=T), digits = d))
}


## Print categorical vars as n (%)
n_pct <- function(x, level=F) {
  if(level==F) {
    sapply(as.list(names(table(x))), function(lvl) {
      paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100,1), "%)") }) } 
  else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100,1), "%)")}
}


## Remove outliers by SD
remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}


## Calculate zscore
zscore.fun <- function(x) {
  z<-((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  return(z)
}


## Median impute for negatie or missing values
median_imp.fun <- function(x) {
  x.new <- ifelse(x == -1 | x == -3 | x == -9 | is.na(x) == T, median(x, na.rm=T), x)
  return(x.new)
}

## Add descriptive labels 
descr_label.fun <- function(data, base_var, labs_vals) {
  
  base <- data %>% select(base_var) 
  temp <- rep(NA, length(base))
  
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}

## ESS F-test 

## For continuous oucomes - Extra Sum of Squares F-Test: function: ess_test.fun()
F_test<-function(out, exp, covars, d=2, data=analysis){
  # Full & reduced models
  full<-lm(formula(paste0(out,"~", exp,"+", covars)), data)
  red<-lm(formula(paste0(out,"~", covars)), data)
  
  # Identify the relevant r/cs in anova output to calcuate F-test
  N_vars<-dim(anova(red))[1]-1 #N variables
  
  #Calculate F-Statistic = ((RSSres-RSSfull)/(Kfull-Kred))/Mean RSSfull
  f<-((deviance(red)-deviance(full))/anova(full)$Df[1])/anova(full)$'Mean Sq'[rownames(anova(full)) == "Residuals"]
  
  #Caluclate P-Values from F-stat (df1=#predictors testing, in total; df2=residual DF in full model)
  p<-pf(q=f, df1=sum(anova(full)$Df[1:N_vars]), df2=4, lower.tail = F) 
  
  ess<-list("F.stat"=f,"P.value"=round(p, 5)); return(ess)
}




