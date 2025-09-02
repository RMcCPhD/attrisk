
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import characteristics
a_imp_bl <- readRDS("processed_data/agg_bl_pooled_arms.rds")

# Helper functions for setting up proportion imputations
softmax <- function(x) {
  ex <- exp(x - max(x, na.rm = TRUE))
  ex / sum(ex)
}

renorm_simplex <- function(p, eps = 1e-3) {
  p <- pmax(eps, pmin(1 - eps, p))
  s <- sum(p)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(p)))
  p / s
}

logit     <- function(p) log(p/(1-p))
inv_logit <- function(z) 1/(1+exp(-z))
alr  <- function(p) log(p[-length(p)] / p[length(p)])
ialr <- function(z) softmax(c(z, 0))

# Format trttype, remove refs for categorical variables
# Remove old variables
# Split out blinding
b_format <- a_imp_bl %>% 
  mutate(
    hba1c_d = case_when(
      trial_id == "NCT00294723" & hba1c_d == 12 ~ 1.2,
      TRUE ~ hba1c_d
    ),
    trttype = factor(trttype),
    blind_ol = if_else(blinding == "open_label", 1L, 0L),
    blind_dbl = if_else(blinding == "double", 1L, 0L),
    blind_trip = if_else(blinding == "triple", 1L, 0L),
    blind_quad = if_else(blinding == "quadruple", 1L, 0L),
  ) %>% 
  select(trial_id:phase, blind_ol:blind_quad, everything()) %>% 
  select(-blinding, -source)

race_cols <- c("asian","black","other","white")

# Function to build dataset with race and smoking status prepped for imputation
prep_df <- function(df) {
  
  # Race log-transformed
  p_race <- as.matrix(df[, race_cols])
  mask_race <- is.finite(p_race)
  p_race[mask_race & p_race <= 0] <- 1e-3
  p_race[mask_race & p_race >= 1] <- 1 - 1e-3
  
  race_mat <- t(apply(p_race, 1, function(p) {
    if (all(is.finite(p))) alr(p) else rep(NA_real_, length(p) - 1)
  }))
  
  race_df <- as_tibble(
    race_mat, 
    .name_repair = ~ paste0("race_alr", seq_len(ncol(race_mat)))
  )
  
  # Smoke - logit
  smoke_cur <- df$current
  smoke_cur[is.finite(smoke_cur) & smoke_cur <= 0] <- 1e-3
  smoke_cur[is.finite(smoke_cur) & smoke_cur >= 1] <- 1 - 1e-3
  smoke_logit <- if_else(is.finite(smoke_cur), logit(smoke_cur), NA_real_)
  
  # Add transformed variables
  df %>% 
    select(-all_of(c(race_cols, "current", "former_never"))) %>% 
    bind_cols(race_df, tibble(smoke_logit = smoke_logit))
}

# Apply function to df
b_amelia_df <- prep_df(b_format)

# Bounds for male proportion
prop_male <- intersect(c("male"), names(b_amelia_df))

prop_bounds <- if (length(prop_male)) {
  do.call(
    rbind,
    lapply(prop_male, function(v) c(match(v, names(b_amelia_df)), 0, 1))
  )
} else {
  matrix(numeric(0), nrow = 0, ncol = 3)
}

# Bounds for numeric variables
b_format %>% 
  reframe(across(ends_with(c("_mu", "_d")), ~ unlist(
    list(min = min(., na.rm = TRUE), max = max(., na.rm = TRUE))
  ))) %>% 
  pivot_longer(everything()) %>% 
  arrange(name) %>% 
  print(n = 32)

marg_limits <- rbind(
  c(match("age_mu", names(b_format)), 39.5, 74.9),
  c(match("age_d", names(b_format)), 3.67, 13.4),
  c(match("bmi_mu", names(b_format)), 24.2, 37.4),
  c(match("bmi_d", names(b_format)), 1.83, 10.6),
  c(match("egfr_mu", names(b_format)), 24.1, 154),
  c(match("egfr_d", names(b_format)), 4.4, 70),
  c(match("hba1c_mu", names(b_format)), 6.7, 10.2),
  c(match("hba1c_d", names(b_format)), 0.4, 6.3),
  c(match("hdl_mu", names(b_format)), 0.94, 2.7),
  c(match("hdl_d", names(b_format)), 0.02, 0.7),
  c(match("sbp_mu", names(b_format)), 124, 151),
  c(match("sbp_d", names(b_format)), 6.7, 17.4),
  c(match("t2dur_mu", names(b_format)), 0.72, 20.7),
  c(match("t2dur_d", names(b_format)), 0.87, 18),
  c(match("total_chol_mu", names(b_format)), 4, 9.8),
  c(match("total_chol_d", names(b_format)), 0.07, 1.9)
)

# Combine bounds for male and numeric covariates
bounds <- rbind(
  prop_bounds, 
  marg_limits[!is.na(marg_limits[,1]), , drop = FALSE]
)

# Run amelia
set.seed(123)

b_impute <- amelia(
  b_amelia_df,
  m = 20,
  idvars  = "trial_id",
  noms    = "trttype",
  empri   = 0.01 * nrow(b_amelia_df),
  bounds  = bounds,
  incheck = TRUE
)

# Function to inverse transform race and smoking to original format
conv_transf <- function(df) {
  alr_cols <- grep("^race_alr", names(df), value = TRUE)
  stopifnot(length(alr_cols) == length(race_cols) - 1)
  
  Z <- as.matrix(df[, alr_cols, drop = FALSE])
  race_rec <- t(apply(Z, 1, function(z) {
    if (all(is.finite(z))) ialr(z) else rep(NA_real_, length(race_cols))
  }))
  race_rec <- t(apply(race_rec, 1, renorm_simplex, eps = 1e-3))
  colnames(race_rec) <- race_cols
  
  # Smoking back from sm_logit (current determines both)
  cur <- inv_logit(df$smoke_logit)
  cur[!is.finite(cur)] <- NA_real_
  cur <- pmin(1 - 1e-3, pmax(1e-3, cur))
  fn  <- 1 - cur
  
  df %>%
    select(-all_of(c(alr_cols, "smoke_logit"))) %>%
    bind_cols(
      as.data.frame(race_rec),
      tibble::tibble(current = cur, former_never = fn)
    )
}

# Apply inverse transformation
b_imp_transf <- lapply(b_impute$imputations, conv_transf)

# Add iteration to each imputed dataset
c_iter <- lapply(seq_along(b_imp_transf), \(j) mutate(b_imp_transf[[j]], imp = j))

# Stack, check race and smoking sum to 1
c_bind <- bind_rows(c_iter) %>% 
  select(imp, trial_id:male, asian:former_never, everything())

# Inspect
c_plot_stack <- c_bind %>% 
  pivot_longer(age_mu:former_never) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.01) +
  facet_wrap(~name, scales = "free")

c_plot_stack

ggsave(
  "outputs/conduct_model/agg_bl/dist_before_clp.png",
  c_plot_stack,
  height = 10,
  width = 10,
  units = "in"
)

# Save
saveRDS(c_bind, "processed_data/agg_bl_imputed_stack.rds")


