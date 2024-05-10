## Basic functions


# load packages
library(tidyverse) ; library(table1)


# Functions

# 1. Print continuous vars as mean +- SD 
mean_sd.fun<-function(x, d=2) {
  sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
          round(sd(x, na.rm=T), digits = d))
}



# 2. Print categorical vars as n (%)
n_pct.fun <- function(x) {
  require(table1)
  (sapply(stats.default(x), function(y) with(y, sprintf("%s (%s)", round(FREQ, d=2), round(PCT, 0)))))
}

n_pct <- function(x, level=F) {
  if(level==F) {
    sapply(as.list(names(table(x))), function(lvl) {
      paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100,1), "%)") }) } 
  else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100,1), "%)")}
}


# 3. Remove outliers by SD
remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}


# 4. Calculate zscore
zscore.fun <- function(x) {
  z<-((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  return(z)
}


# 5. Median impute for negatie or missing values
median_imp.fun <- function(x) {
  x.new <- ifelse(x == -1 | x == -3 | x == -9 | is.na(x) == T, median(x, na.rm=T), x)
  return(x.new)
}

# 5. Function to create variable with descriptive labels 
descr_label.fun <- function(data, base_var, labs_vals) {
  
  base <- data %>% select(base_var) 
  temp <- rep(NA, length(base))
  
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}




