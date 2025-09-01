
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import bl chars
a_imp_bl <- readRDS("processed_data/agg_bl_fourth_ext.rds")
a_imp_set <- readRDS("processed_data/latest_dataset.rds")

# Add class
a_imp_class <- readRDS("processed_data/agg_bl_first_ext.rds")

a_add_class <- a_imp_class %>% 
  select(trial_id:trtcls5) %>% 
  left_join(a_imp_bl)

# Manually add baseline chars for two trials not in set
# write_csv(a_add_class, "scratch_data/add_n2_manually.csv")
a_imp_fix <- a_add_class %>% 
  filter(!trial_id %in% c("NCT02730377", "NCT02911948")) %>% 
  full_join(read_csv("scratch_data/added_n2_manually.csv")) %>% 
  arrange(trial_id, arm_id)

# Add treatment type
# Remove experimental arms
# Collapse trials with placebo and ac into one trttype
b_add_trttype <- a_imp_fix %>% 
  mutate(
    trttype = case_when(
      class %in% c("dpp4", "glp1", "sglt2") ~ NA,
      class == "placebo" ~ 0L,
      !class == "placebo" & !class %in% c("dpp4", "glp1", "sglt") ~ 1L
    )
  ) %>% 
  filter(!is.na(trttype)) %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      length(unique(trttype)) > 1 ~ 2L,
      TRUE ~ trttype
    )
  ) %>% 
  select(trial_id:arm_id, trttype, everything())

# Helper function - compute pooled SD when summarising group chars
pooledSD <- function(mu, delta, n) {
  # Get weighted mean
  keep_mu <- is.finite(mu) & is.finite(n)
  mu <- mu[keep_mu]
  n_mu <- n[keep_mu]
  if (length(n_mu) == 0) return(NA_real_)
  N <- sum(n_mu)
  if (N <= 1) return(NA_real_)
  Mu <- sum(n_mu * mu) / N
  
  # Get pooled SD
  keep_delta <- is.finite(delta) & is.finite(n)
  within_sd <- sum( (n[keep_delta] - 1) * (delta[keep_delta]^2) )
  bw_sd <- sum( n_mu * (mu - Mu)^2 )
  sqrt( (within_sd + bw_sd) / (N - 1) )
}

# Add sample size to get weighted means and sds collapsing trttype
b_imp_n <- b_add_trttype %>% 
  left_join(read_csv("scratch_data/aact_n.csv")) %>% 
  select(trial_id:trtcls5, n, everything())

# Manually add two missing
# write_csv(
#   b_imp_n %>% filter(is.na(n)),
#   "scratch_data/aact_n_manual.csv"
# )

# Add manual n
b_add_n2 <- b_imp_n %>% 
  filter(!trial_id %in% c("NCT02730377", "NCT02911948")) %>% 
  full_join(read_csv("scratch_data/aact_n_manual_added.csv")) %>% 
  arrange(trial_id, arm_id)

# Summarise across trttype per trial
# Weighted mean and sd accounting for within-group variation
# Sum race counts
b_sum_trttype <- b_add_n2 %>% 
  group_by(trial_id, trttype) %>% 
  reframe(
    arm_n = sum(n),
    across(asian:white, ~ sum(.)),
    male = sum(male),
    across(current:former_never, ~ sum(.)),
    across(contains("mu"), ~ weighted.mean(., n)),
    age_d = pooledSD(age_mu, age_d, n),
    sbp_d = pooledSD(sbp_mu, sbp_d, n),
    bmi_d = pooledSD(bmi_mu, bmi_d, n),
    t2dur_d = pooledSD(t2dur_mu, t2dur_d, n),
    egfr_d = pooledSD(egfr_mu, egfr_d, n),
    hba1c_d = pooledSD(hba1c_mu, hba1c_d, n),
    hdl_d = pooledSD(hdl_mu, hdl_d, n),
    total_chol_d = pooledSD(total_chol_mu, total_chol_d, n)
  ) %>% 
  mutate(across(age_mu:total_chol_d, ~ round(., 2))) %>% 
  select(
    trial_id:trttype, arm_n, contains("age"), contains("sbp"), contains("bmi"),
    contains("t2dur"), contains("egfr"), contains("hba1c"), contains("hdl"),
    contains("total_chol"), everything()
  )

# Join remaining groups from design dataset
b_add_vars <- b_sum_trttype %>% 
  arrange(trial_id) %>% 
  cbind(
    a_imp_set %>% 
      arrange(trial_id) %>% 
      select(-trial_id, -n, -arm_n, -trttype, -attr_prop)
  ) %>% 
  as_tibble() %>% 
  select(source, trial_id, trttype, attr, everything())

# Check missingness
b_add_vars %>% 
  reframe(across(trial_id:former_never, ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  arrange(desc(prop)) %>% 
  print(n = 25)

# Setup for amelia -------------------------------------------------------------

# Format trttype, remove refs for categorical variables
# Remove old variables
# Split out blinding
b_format <- b_add_vars %>% 
  mutate(
    trttype = factor(trttype),
    blind_ol = if_else(blinding == "open_label", 1L, 0L),
    blind_dbl = if_else(blinding == "double", 1L, 0L),
    blind_trip = if_else(blinding == "triple", 1L, 0L),
    blind_quad = if_else(blinding == "quadruple", 1L, 0L),
  ) %>% 
  select(-c(
    white, former_never, 
    age_lower_dbl:age_upper_dbl, 
    bmi_lower:bmi_upper,
    male_prop:egfr_mean
  )) %>% 
  select(trial_id:phase, blind_ol:blind_quad, everything()) %>% 
  select(-blinding, -source)

saveRDS(b_format, "processed_data/agg_bl_before_imp.rds")

# Set bounds for proportion imputations
prop_vars <- intersect(
  c("male", "asian", "black", "other", "current"), 
  names(b_format)
)

bounds <- do.call(
  rbind, 
  lapply(
    prop_vars, 
    function(v) c(match(v, names(b_format)), 0, 1)
  )
)

# Set bounds for mean and sd imputations
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
  c(match("hba1c_d", names(b_format)), 0.4, 12),
  c(match("hdl_mu", names(b_format)), 0.94, 2.7),
  c(match("hdl_d", names(b_format)), 0.02, 0.7),
  c(match("sbp_mu", names(b_format)), 124, 151),
  c(match("sbp_d", names(b_format)), 6.7, 17.4),
  c(match("t2dur_mu", names(b_format)), 0.72, 20.7),
  c(match("t2dur_d", names(b_format)), 0.87, 18),
  c(match("total_chol_mu", names(b_format)), 4, 9.8),
  c(match("total_chol_d", names(b_format)), 0.07, 1.9)
)

bounds <- rbind(
  bounds, 
  marg_limits[!is.na(marg_limits[,1]), , drop = FALSE]
)

# Run amelia -------------------------------------------------------------------

set.seed(123)

b_impute <- amelia(
  b_format,
  m = 20,
  idvars  = "trial_id",
  noms    = "trttype",
  empri   = 0.01 * nrow(b_format),
  bounds  = bounds,
  incheck = TRUE
)

# Fix race and smoking groups imputed, sum to 1
fix_grps <- function(df) {
  
  # Smoking groups
  df$current <- pmin(pmax(df$current, 0), 1)
  df$former_never <- 1 - df$current
  
  # Race groups
  race_cols <- c("asian", "black", "other")
  race_X <- as.matrix(df[race_cols])
  race_X[!is.finite(race_X)] <- NA
  race_X[is.na(race_X)] <- 0
  race_X[race_X < 0] <- 0
  race_X[race_X > 1] <- 1
  
  race_s <- rowSums(race_X)
  zero <- !is.finite(race_s) | race_s <= 0
  if (any(zero)) race_X[zero, ] <- 1/3
  
  race_s <- rowSums(race_X)
  race_big <- race_s > 1
  if (any(race_big)) race_X[race_big, ] <- race_X[race_big, ] / race_s[race_big]
  
  df[race_cols] <- race_X
  df$white <- 1 - rowSums(race_X)
  df$white <- pmin(pmax(df$white, 0), 1)
  
  df %>%
    select(trial_id:other, white, everything()) %>% 
    mutate(
      across(age_mu:former_never, ~ round(., 2)),
      race_total = asian + black + other + white,
      white = case_when(
        race_total == 0.99 ~ white + 0.01,
        race_total == 1.01 ~ white - 0.01,
        TRUE ~ white
      ),
      race_total = asian + black + other + white,
      smoke_total = current + former_never
    )
}

b_fixed <- lapply(b_impute$imputations, fix_grps)
saveRDS(b_fixed, "processed_data/agg_bl_imp_lst.rds")

# Add imputation iteration to datasets, bind rows
c_bind <- bind_rows(
  lapply(seq_along(b_fixed), \(j) mutate(b_fixed[[j]], imp = j))
)

# Final fix for race caused by rounding
c_fix <- c_bind %>% 
  mutate(
    race_total = round(race_total, 2),
    white = case_when(
      race_total == 1.01 ~ white - 0.01,
      race_total == 0.99 ~ white + 0.01,
      TRUE ~ white
    ),
    race_total = asian + black + other + white
  ) %>% 
  select(-race_total, -smoke_total) %>% 
  select(trial_id:current, former_never, everything())

# Save
saveRDS(c_fix, "processed_data/agg_bl_design_imputed_iters.rds")
