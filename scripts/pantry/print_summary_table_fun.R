## Print Summary Table


print_summary_table <- function(
    data, # dataframe
    vars_to_summarize, # vector of variables and labels; c(var1 = "Label1, unit", var2 = "Label2, unit")
    var_strata = F, # name of strata variable; "taste_diplos_AA"
    var_strata_order=F, # F = based on variable coding in data
    factor_vars=F,
    digits=c(2,1,4), # d_cont; d_pct; d_pval
    p_print=T, # If T, print P-values in summmary table; default = F
    p_types=c("descriptive", "trend"), # P-value types; descriptive(t.test-P/ANOVA-F for cont; X2-P for cat); trend(continuous exposure levels using lm)
    p_adjust="agesex", # P-value adjustments; default = c(agesex); otherwise, list variables for adjustment: c("age", "sex", gPC1", ...)
    hide_bin_reference=T
    ) {
  
  ############################################################
  ##   Build required printing functions: mean_sd & n_pct   ##
  ############################################################
  
  # ======================================
  ## For continuous vars: mean_sd
  # ======================================
  mean_sd<-function(x, d=2) {
    sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
            round(sd(x, na.rm=T), digits = d))
  }
  
  
  # ==================================
  ## Print categorical vars as n (%)
  # ==================================
  n_pct <- function(x, level=F, d=1) {
    if(level==F) {
      sapply(as.list(names(table(x))), function(lvl) {
        paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100, d), "%)") }) } 
    else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100, d), "%)")}
  }
  
  # ===========================
  # Format P-value function
  # ===========================
  format_p.fun <- function(p, digits) {
    P<-ifelse(as.numeric(p) <0.01, format(as.numeric(p), scientific = T, digits=2), 
              round(as.numeric(p), digits))
    return(P)
  }
  
  
  ##############################################
  ##  Define inputs & Load required packages  ##
  ##############################################
  
  ## Packages ==============
  library(tidyverse) ; library(purrr)
  
  
  ## Inputs =================
  
  #Digits for printing (default = c(1,1))
  d_cont = digits[1] ; d_pct = digits[2] ; d_pval=digits[3]
  
  #Pvalue adjustments
  if(p_adjust[1] == "agesex") { p_adjust_vars = c("age", "sex") } else {
    p_adjust_vars = p_adjust }
  
  #Strata var order, if different from variable coding (default = F)
  if(var_strata==F) { data <- data %>% 
    mutate(random_strata = sample(c(0,1), size=nrow(.), replace=T)) ; 
  var_strata = "random_strata" ; var_strata_order=c(0,1) } else if(var_strata != F) {
    var_strata_order = if(var_strata_order[1] == F) {
    c(na.omit(unique(data[, ..var_strata][[1]]))) } else{ var_strata_order }
  }
  
  #Make dataframe with strata_var (no missing; ordered by var_strata_order) & vars_to_summarise
  dat_total <- data %>% 
    dplyr::select(strata=all_of(var_strata), all_of(names(vars_to_summarize)), all_of(p_adjust_vars)) %>%
    mutate(strata_ordered = factor(strata, levels=var_strata_order)) %>%
    filter(complete.cases(strata_ordered))
  
  #If P-trend should also be printed ... make a continuous strata exposure
  if("trend" %in% p_types) {
    dat_total <- dat_total %>%
      mutate(strata_ordered.cont = as.numeric(strata_ordered))
  }
  
  ## Check that all factor_vars are coded as factors
  if(factor_vars[1] != F) {
    dat_total <- dat_total %>% 
      mutate(across(all_of(factor_vars), ~as.factor(.)))
  } 
  
  #Make dataframe grouped by strata
  dat_grouped <- dat_total %>%
    group_by(strata_ordered)
  
  #Save strata level labels for table
  strata_var_lvls <- levels(dat_grouped$strata_ordered)
  
  footnotes <- list()
  
  
  
  #####################################################
  ##   Summarize each variable, grouping by STRATA   ##
  #####################################################

  summary_by_strata <-
    lapply(1:length(vars_to_summarize), function(v) {
    
    ## Make temporary vars for name, data, and Label
    var = names(vars_to_summarize)[v]
    var.total.dat = dat_total %>% dplyr::select(strata_ordered, Var=all_of(var), all_of(p_adjust_vars)) %>%
      filter(complete.cases(Var))
      #mutate(var_na=recode_na_as.fun(Var) ) %>% filter(var_na == 1)
    var.grouped.dat = dat_grouped %>% dplyr::select(strata_ordered, Var=all_of(var)) %>%
      filter(complete.cases(Var))
    var.Label=paste0(vars_to_summarize[v][[1]])
  
    # ==============================
    ## If var is continuous, mean/SD
    # ==============================
    
    if(is.numeric(var.total.dat$Var) == T) {
      var.summary <- cbind.data.frame(
        var.total.dat %>% reframe(Var=mean_sd(Var, d=d_cont)) %>% 
        t() %>% as.data.frame(), 
        (var.grouped.dat %>% reframe(Var=mean_sd(Var, d=d_cont)) %>% 
        t() %>% as.data.frame())[-1,]) 
      
      # Reformatting
      colnames(var.summary) <- c("Total", var_strata_order)
      rownames(var.summary) <- var.Label
      
      ## Add P-values using ANOVA (for ≥3-level strata)
      if(p_print == T) {
        p_to_print <- list()
        if(var %in% p_adjust_vars) {
          if("descriptive" %in% p_types) {
            P_Ftest=format_p.fun(anova(lm(formula(paste0(var, "~strata_ordered", "+", paste0(p_adjust_vars[-which(p_adjust_vars==var)], collapse = "+"))), data=dat_total))$'Pr(>F)'[1], d_pval) }
          if("trend" %in% p_types) {
            P_trend=format_p.fun(summary(lm(formula(paste0(var, "~strata_ordered.cont", "+", paste0(p_adjust_vars[-which(p_adjust_vars==var)], collapse = "+"))), data=dat_total))$coef[2,4], d_pval) }
        } else {
          if("descriptive" %in% p_types) {
            P_Ftest=format_p.fun(anova(lm(formula(paste0(var, "~strata_ordered", "+", paste0(p_adjust_vars, collapse = "+"))), data=dat_total))$'Pr(>F)'[1], d_pval) }
          if("trend" %in% p_types) {
            P_trend=format_p.fun(summary(lm(formula(paste0(var, "~strata_ordered.cont", "+", paste0(p_adjust_vars, collapse = "+"))), data=dat_total))$coef[2,4], d_pval) }
        } 
        
        # Select which P-values to print
        if("descriptive" %in% p_types & "trend" %in% p_types) {
          p_to_print <- c("P_value"=P_Ftest, "P_trend"=P_trend)
          } else if(p_types == "descriptive") {p_to_print <- c("P_value"=as.character(P_Ftest))} else if(p_types == "trend") {p_to_print <- c("P_trend"=as.character(P_trend)) }
        
        var.summary <- cbind.data.frame(var.summary, t(as.data.frame(p_to_print))) %>% mutate(across(starts_with("P_"), ~as.character(.)))
      } else {var.summary <- var.summary}
      
    }
    
    # ==============================
    ## If var is categorical, n (%)
    # ==============================
    
    if(is.numeric(var.total.dat$Var) == F) {
      if(is.factor(var.total.dat$Var) == T) {
        var_lvls <- levels(var.total.dat$Var) } else {
          var_lvls <- unique(var.total.dat$Var) 
        } ;  var.summary <- cbind.data.frame(rbind(
        NA, do.call(rbind.data.frame, lapply(var_lvls, function(lvl) {
          var.total.dat %>% 
            reframe(Total.Var=n_pct(Var, level=lvl, d=d_pct)) 
          }) )),
        lapply(var_lvls, function(lvl) {
          var.grouped.dat %>% reframe(Var=n_pct(Var, level=lvl, d=d_pct)) 
          }) %>% reduce(left_join, by = "strata_ordered") %>% 
          t() %>% as.data.frame() 
        ) 
      
      # Reformatting
      colnames(var.summary) <- c("Total", var_strata_order)
      rownames(var.summary) <- c(var.Label, paste0(" ", var_lvls))
      var.summary[1,]<-" "
      
      ## Add P-values
      if(p_print == T) {
        p_to_print <- list()
        P_ChiSq=format_p.fun(chisq.test(dat_total[, ..var][[1]], dat_total$strata_ordered)$p.value, d_pval)
        
        #Select which P-values to print
        if("descriptive" %in% p_types & "trend" %in% p_types) {
          p_to_print <- cbind.data.frame(P_value=c(P_ChiSq, rep("", length(var_lvls))), P_trend=c("-", rep("", length(var_lvls))))
        } else if(p_types == "descriptive" | p_types == "trend") { p_to_print <- data.frame(P_value=c(P_ChiSq, rep("", length(var_lvls)))) }
          
        #var.summary <- cbind.data.frame(var.summary, t(as.data.frame(P_to_print)))
        var.summary <- cbind.data.frame(var.summary, p_to_print)
      } else {var.summary <- var.summary}
      
      # If binary variable and hide_bin_reference == T
      if(hide_bin_reference==T) {
        if(length(var_lvls)==2 & var_lvls[1] %in% c("0", "1", F, T, "No", "Yes")){
          # Pick row corresponding to positive binary variable: 1 or T or Yes
          show <- var_lvls[length(var_lvls)] ; merge_cols <- c(which(startsWith(colnames(var.summary), "P_")==F))
          var.summary[1,c(merge_cols)] <- var.summary[which(rownames(var.summary) == paste0(" ",show)),c(merge_cols)]
          var.summary <- var.summary[1,] 
        }
      }
    } 
    
    return(var.summary)

  }) %>% do.call(rbind.data.frame, .)
  
  N=matrix(NA, 1, ncol(summary_by_strata), dimnames=list("N (%)", colnames(summary_by_strata))) 
  N[1,1:(1+length(strata_var_lvls))] <- rbind.data.frame(dat_total %>% reframe(N=round(n(), d=d_pct)), dat_total %>% reframe(N = gsub(".*[,]", "", n_pct(strata_ordered, d=d_pct)))) %>% t()
  N <- as.data.frame(N) %>% mutate(across(starts_with("P_"), ~as.character(ifelse(is.na(.), "-", .))))
  
  summary_by_strata <- bind_rows(N, summary_by_strata)
  
  if(var_strata == "random_strata") {
    return(summary_by_strata %>% dplyr::select(Total))
  } else {
    return(summary_by_strata)
  }
  
}

## END OF SCRIPT
