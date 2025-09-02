
# Pool distributions for each variable based on best-fitting
# Helper functions written to pool for each trial

source("scripts/00_config.R")
library(tidyverse)

# Import stacked imputations
a_imp_df <- readRDS("processed_data/agg_bl_imputed_stack.rds") %>% 
  mutate(across(age_mu:former_never, ~ round(., 2)))

# Find best distribution per variable and pool accordingly ---------------------

# Parametric distributions
z_dists <- c("norm", "lnorm", "gamma", "beta", "weibull")

# Function to fit distribution to each variable and get best
find_best_dist <- function(x, families = z_dists) {
  x <- x[is.finite(x)]
  fits <- list()
  
  # Conditionally apply
  
  if ("norm" %in% families) {
    fits$norm  <- fitdistrplus::fitdist(x, "norm")
  }
  
  if ("lnorm" %in% families && all(x > 0)) {
    fits$lnorm <- fitdistrplus::fitdist(x, "lnorm")
  }
  
  if ("gamma" %in% families && all(x > 0)) {
    fits$gamma <- fitdistrplus::fitdist(x, "gamma")
  }
  
  if ("beta" %in% families && all(x > 0 & x < 1)) {
    # Avoid boundary issues (i.e. exactly 0 or 1)
    eps <- 1e-6
    x2  <- pmin(pmax(x, eps), 1 - eps)
    fits$beta <- fitdistrplus::fitdist(x2, "beta")
  }
  
  if ("weibull" %in% families && all(x > 0)) {
    fits$weibull <- fitdistrplus::fitdist(x, "weibull")
  }
  
  # Get AICs
  aics <- map_dbl(fits, ~ .x$aic)
  fam <- names(aics)[which.min(aics)]
  list(family = fam, fit = fits[[fam]])
}

# Apply function to continuous variables
# Stacked dataset (i.e. bind_rows each imputation iteration)
# Output is pooled estimates across iterations using best distribution
z_vars <- names(a_imp_df[c(6:28)])
# z_vars <- names(a_imp_df[c(6:22, 27:28)])

b_best_dists <- map(z_vars, function(v) {
  best <- find_best_dist(a_imp_df[[v]])
  tibble(
    var = v,
    family = best$family,
    aic = best$aic,
    estimates = list(if (is.null(best$fit)) NA else best$fit$estimate)
  )
})

# Tidy outputs
b_best_df <- bind_rows(b_best_dists) %>% 
  separate(estimates, into = c("est1", "est2"), sep = ", ") %>% 
  mutate(across(est1:est2, ~ gsub("c\\(|\\(|\\)", "", .))) %>% 
  pivot_longer(est1:est2) %>% 
  separate(value, into = c("dist", "value"), sep = " = ") %>% 
  select(-name)

# Convert non-normal parameters to original scale (i.e. mean/sd)
c_wide <- b_best_df %>% 
  pivot_wider(names_from = dist, values_from = value) %>% 
  mutate(across(meanlog:shape2, as.numeric))

# Get expected values for variables based on best-fitting distribution
get_values_family <- function(fam, p) {
  # Get distribution parameters
  switch(
    fam,
    "norm"    = p$mean,                           
    "lnorm"   = exp(p$meanlog + 0.5 * p$sdlog^2),     
    "gamma"   = if (!is.na(p$rate)) p$shape / p$rate else p$shape * p$scale,
    "weibull" = p$scale * gamma(1 + 1/p$shape), 
    "beta"    = p$shape1 / (p$shape1 + p$shape2),
    stop(sprintf("Unsupported family: %s", fam))
  )
}

# Pull variable row and compute expected value under distribution
get_exp_values <- function(df, varname) {
  row <- df %>% filter(var == varname)
  if (nrow(row) == 0) return(NA_real_)
  fam <- row$family[1]
  p <- as.list(row)
  get_values_family(fam, p)
}

# Build dataset of pooled variable parameters for continuous variables and male
d_df_exp <- tibble(vars = unique(c_wide$var)) %>% 
  mutate(
    values = map_dbl(vars, ~ get_exp_values(c_wide, .x)),
    values = round(values, 2)
  )

# Save
saveRDS(d_df_exp, "processed_data/agg_bl_dists_ready.rds")
