
# Pool distributions for each variable based on best-fitting
# Helper functions written to pool for each trial

source("scripts/00_config.R")
library(tidyverse)

# Import stacked imputations
a_imp_df <- readRDS("processed_data/agg_bl_imputed_stack.rds") %>% 
  mutate(across(age_mu:former_never, ~ round(., 2)))

# Nest
a_nest <- a_imp_df %>% 
  group_by(trial_id) %>% 
  nest() %>% 
  mutate(
    data_rest = map(data, function(df) {
      df %>% 
        select(age_mu:former_never) %>% 
        select(where(~n_distinct(.) == 1))
    })
  )

# Keep only variables that were imputed
a_keep <- a_nest %>% 
  mutate(
    data_filter = map(data, function(df) {
      df %>% 
        select(age_mu:former_never) %>% 
        select(where(~n_distinct(.) > 1))
    }),
  ) %>% 
  filter(
    map_lgl(data_filter, ~ ncol(.x) > 0)
  )

# Find best distribution per variable and pool accordingly ---------------------

# Parametric distributions
z_dists <- c("norm", "lnorm", "weibull", "beta")

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

# Variables
z_vars <- names(a_imp_df[c(6:28)])

# Apply function to continuous variables per trial
# Output is pooled estimates across iterations using best distribution
b_best_dists <- a_keep %>% 
  transmute(
    trial_id,
    best_dist = map(
      data_filter, function(df) {
        vars_pres <- intersect(z_vars, names(df))
        if (!length(vars_pres)) return(tibble())
        map_df(
          vars_pres, function(v) {
            best <- find_best_dist(df[[v]])
            tibble(
              var = v,
              family = best$family,
              aic = best$aic,
              estimates = list(if (is.null(best$fit)) NA else best$fit$estimate)
            )
          }
        )
      }
    )
  ) %>% 
  unnest(best_dist)

b_best_df <- b_best_dists %>% 
  unnest_wider(estimates, names_repair = "unique") %>% 
  mutate(across(meanlog:shape2, as.numeric))

# Get expected values for variables based on best-fitting distribution
get_values_family <- function(fam, p) {
  # Get distribution parameters
  switch(
    fam,
    "norm"    = p$mean,                           
    "lnorm"   = exp(p$meanlog + 0.5 * p$sdlog^2),     
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
# Scale to match participant-level covariate scaling
d_df_exp <- b_best_df %>% 
  group_by(trial_id) %>% 
  group_modify(~{
    dat <- .x
    out <- tibble(vars = unique(dat$var)) %>% 
      mutate(
        values = map_dbl(vars, ~ get_exp_values(dat, .x)),
        values = round(values, 2)
      )
  }) %>% 
  ungroup()

# Join with variable marginals that were not imputed
# All length 23
e_not_imp <- a_nest %>% 
  select(trial_id, data_rest) %>% 
  unnest(data_rest) %>% 
  distinct() %>% 
  pivot_longer(age_mu:egfr_d) %>% 
  filter(!is.na(value)) %>% 
  rename(vars = name, values = value) %>% 
  full_join(d_df_exp) %>% 
  arrange(trial_id, vars)

e_not_imp %>% count(trial_id) %>% filter(n != 23)

# Scale
f_df_scale <- e_not_imp %>% 
  mutate(
    values_scaled = case_when(
      grepl("age", vars) ~ values / 10,
      vars == "sbp_mu"   ~ values / values[vars == "sbp_d"],
      vars == "bmi_mu"   ~ values / values[vars == "bmi_d"],
      vars == "t2dur_mu" ~ values / values[vars == "t2dur_d"],
      vars == "egfr_mu"  ~ values / values[vars == "egfr_d"],
      vars %in% c("sbp_d", "bmi_d", "t2dur_d", "egfr_d") ~ 1,
      TRUE ~ values
    )
  )

# Save
saveRDS(f_df_scale, "processed_data/agg_bl_dists_ready.rds")
