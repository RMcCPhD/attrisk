
# Pool distributions for each variable based on best-fitting
# Helper functions written to pool for each trial

source("scripts/00_config.R")
library(tidyverse)

# Import stacked imputations
a_imp_df <- readRDS("processed_data/agg_bl_imputed_stack.rds") %>% 
  mutate(across(age_mu:former_never, ~ round(., 2)))

# Check t2dur marginals where negative values produced
a_chk_t2d <- a_imp_df %>%
  inner_join(
    read_csv("scratch_data/sobol_neg_vals.csv") %>% 
      select(trial_id)
  ) %>% 
  select(trial_id, contains("t2dur"))

# Nest
a_nest <- a_imp_df %>% 
  group_by(trial_id) %>% 
  nest() %>% 
  mutate(
    data_rest = map(data, function(df) {
      df %>% 
        select(age_mu:conmeds_size) %>% 
        select(where(~n_distinct(.) == 1))
    })
  )

# Keep only variables that were imputed
a_keep <- a_nest %>% 
  mutate(
    data_filter = map(data, function(df) {
      df %>% 
        select(age_mu:conmeds_size) %>% 
        select(where(~n_distinct(.) > 1))
    }),
  ) %>% 
  filter(
    map_lgl(data_filter, ~ ncol(.x) > 0)
  )

# Find best distribution per variable and pool accordingly ---------------------

# Variables that must be positive (negative values introduced using norm)
z_pos <- c("egfr", "hdl", "t2dur", "total_chol")
eps <- 1e-8

# Parametric distributions
z_dists <- c("norm", "lnorm", "gamma", "beta")

# Select appropriate families
fams_vars <- function(varname, families = z_dists) {
  base <- sub("_(mu|d)$", "", varname)
  if (base %in% z_pos) {
    setdiff(families, "norm")
  } else {
    families
  }
}

# Safely fit distributions
fit_safely <- function(dist, x) {
  tryCatch(
    suppressWarnings(fitdistrplus::fitdist(x, dist)),
    error = function(e) NULL
  )
}

# Function to fit distribution to each variable and get best
find_best_dist <- function(x, varname, families = z_dists) {
  x <- x[is.finite(x)]
  fams <- fams_vars(varname, families)
  fits <- list()

  # Conditionally fit
  if ("norm" %in% fams) fits$norm <- fit_safely("norm", x)
  if ("lnorm" %in% fams && all (x > 0)) fits$lnorm <- fit_safely("lnorm", x)
  if ("gamma" %in% fams && all (x > 0)) fits$gamma <- fit_safely("gamma", x)
  if ("beta" %in% fams && all(x > 0 & x < 1)) {
    x2 <- pmin(pmax(x, 1e-6), 1 - 1e-6)
    fits$beta <- fit_safely("beta", x2)
  }
  
  # Drop any failed fits
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop(sprintf("No dist fit for %s", varname))
  aics <- map_dbl(fits, ~ .x$aic)
  fam <- names(aics)[which.min(aics)]
  list(family = fam, fit = fits[[fam]], aic = unname(aics[fam]))
}

# Variables
z_vars <- names(a_imp_df[c(6:32)])

# Apply function per trial and variable
b_best_dists <- a_keep %>% 
  transmute(
    trial_id,
    best_dist = map(
      data_filter, function(df) {
        vars_pres <- intersect(z_vars, names(df))
        if (!length(vars_pres)) return(tibble())
        map_df(
          vars_pres, function(v) {
            best <- find_best_dist(df[[v]], varname = v)
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
  mutate(across(mean:shape2, as.numeric))

# Get expected values per family
get_values_family <- function(fam, p) {
  # Get distribution parameters
  switch(
    fam,
    "norm"    = p$mean,                           
    "lnorm"   = exp(p$meanlog + 0.5 * p$sdlog^2),   
    "gamma"   = p$shape / p$rate,
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
# All length 27
e_not_imp <- a_nest %>% 
  select(trial_id, data_rest) %>% 
  unnest(data_rest) %>% 
  distinct() %>% 
  pivot_longer(age_mu:egfr_d) %>% 
  filter(!is.na(value)) %>% 
  rename(vars = name, values = value) %>% 
  full_join(d_df_exp) %>% 
  arrange(trial_id, vars)

e_not_imp %>% count(trial_id) %>% filter(n != 27)

# Get marginals where any were categorical counts
f_marg_cat <- e_not_imp %>% 
  pivot_wider(names_from = "vars", values_from = "values") %>% 
  mutate(
    egfr_mu = if_else(trial_id == "NCT03689374", 91.8, egfr_mu),
    egfr_d = if_else(trial_id == "NCT03689374", 20.4, egfr_d)
  ) %>% 
  pivot_longer(age_d:white, names_to = "vars", values_to = "values")

# Scale
f_df_scale <- f_marg_cat %>% 
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

# Plot
plot_scaled <- f_df_scale %>% 
  mutate(
    vars = case_when(
      grepl("asian|black|other|white", vars) ~ paste0("race_", vars),
      grepl("curr|form", vars) ~ paste0("smoke_", vars),
      TRUE ~ vars
    )
  ) %>% 
  ggplot(aes(x = values)) +
  geom_density(alpha = 0.01) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~vars, scales = "free") +
  theme_classic()

plot_scaled

# Save
saveRDS(f_df_scale, "processed_data/agg_bl_dists_ready.rds")
