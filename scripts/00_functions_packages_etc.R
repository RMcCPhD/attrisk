
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
