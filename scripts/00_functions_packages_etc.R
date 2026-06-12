
rm(list = ls())
gc()

library(tidyverse)
library(broom)
library(rstan)
library(bayesplot)
library(tidybayes)
library(posterior)
library(markdown)
library(flexsurv)

options(scipen = 999)
Sys.setenv(STANC3 = TRUE)

# Function to check dataset missingness quickly
check_nas <- function(dat) {
  dat %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>% 
    pivot_longer(everything())
}

# Function to find the nearest neighbour IPD trial to an aggregate trial
# Use the fu_ratio of the nearest neighbour as an approximate fu ratio for the
# respective aggregate trials
nn_ratio <- function(row, ipd) {
  
  # Calculated squared Euclidean distances
  distances <- (ipd$scale_arm_n - row$scale_arm_n) ^ 2 +
               (ipd$scale_arm_attr - row$scale_arm_attr) ^ 2 +
               (ipd$scale_dur - row$scale_dur) ^ 2
  
  # Get match (smallest distance)
  best_match <- which.min(distances)
  
  # Return mft ratio from best match
  return(ipd$mft_ratio[best_match])
  
}

# Function to produce a tidy count of each variable, informing coding
# Contingency for dummy variables, summary for numeric
tidy_sum <- function(df, var_name) {
  var <- df[[var_name]]
  cat("Variable:", var_name, "\n")
  if (is.integer(var)) {
    tbl <- table(var) %>% as.data.frame()
    print(tbl)
  } else if (is.numeric(var)) {
    tbl <- summary(var)
    print(tbl)
  } else {NA}
}

# Function for covariate posterior predictive checks
plot_subgroup_ppcs <- function(y, yrep, covariate_df, covariate_names, stat = "mean") {
  stopifnot(length(y) == ncol(yrep))
  
  plots <- list()
  
  for (var in covariate_names) {
    group <- covariate_df[[var]]
    
    # Skip if group has too many levels
    if (length(unique(group)) > 10) {
      warning(glue::glue("Skipping {var} (too many levels)"))
      next
    }
    
    p <- ppc_stat_grouped(
      y = y,
      yrep = yrep,
      group = group,
      stat = stat
    ) +
      ggtitle(glue::glue("PPC by {var}")) +
      theme_minimal()
    
    plots[[var]] <- p
  }
  
  return(plots)
}
