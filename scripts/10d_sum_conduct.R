
# Summary statistics and plots
# Conduct-stage model

rm(list = ls())
gc()

source("scripts/00_config.R")
library(tidyverse)
library(posterior)
library(bayesplot)
library(flexsurv)
library(metafor)
library(broom)

# Sobol: Arm-level linear predictors and hazard multipliers
sobol_imp_sum <- read_csv("outputs/conduct/sobol/sobol_linear_predictors.csv")

sobol_imp_sum %>% 
  mutate(across(-trial_id, ~ round(., 5))) %>% 
  write_csv("../../Attrisk/GitHub/conduct_sobol_linear_predictors_hazard_multipliers.csv")

# CITL: Posterior summaries
citl_imp_arm <- read_csv("outputs/conduct/citl/citl_arm_summary.csv")

citl_imp_arm %>% 
  left_join(sobol_imp_sum %>% distinct(trial_id, arm)) %>% 
  select(trial_id, arm, everything()) %>%
  mutate(across(-c(trial_id, arm, N, r_obs), ~ round(., 5))) %>% 
  write_csv("../../Attrisk/GitHub/conduct_citl_posterior_arm_summary.csv")

# Case scenario attrition risk prediction in remainder of script
# Predict attrition risk trajectories using the calibrated conduct-stage model
# Applied to hypothetical participant profiles in a fixed trial design
# Uses typical predictions as the model is not externally validated and
# cannot reliably make marginal predictions for a "new" trial

# Helper functions -------------------------------------------------------------

summ_vec <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE)
  )
}

gompertz_risk <- function(rate, shape, time) {
  # Vectorised over posterior draws.
  # rate = Gompertz rate parameter, shape = Gompertz shape parameter.
  # If shape is close to zero, use exponential limit.
  ifelse(
    abs(shape) < 1e-8,
    1 - exp(-rate * time),
    1 - exp(-(rate / shape) * (exp(shape * time) - 1))
  )
}

set_x_value <- function(X, nm, value) {
  if (nm %in% names(X)) {
    X[nm] <- value
  } else {
    warning("Design matrix column not found: ", nm)
  }
  X
}

# Load design model draws ------------------------------------------------------

f_design_mdl <- readRDS("outputs/design/design_model.rds")
f_design_df  <- readRDS("processed_data/stan_data.rds")

f_design_draws <- posterior::as_draws_matrix(as.array(f_design_mdl))

beta_cols <- grep("^beta\\[", colnames(f_design_draws), value = TRUE)

beta_draws   <- f_design_draws[, beta_cols, drop = FALSE]
mu_mu1_draws <- f_design_draws[, "mu_mu1"]

# Shape draws: use typical/average shape across AgD arms for the illustrative
# scenario, consistent with a typical-trial prediction rather than a new-trial
# prediction exercise.
shape_cols <- grep("^a\\[", colnames(f_design_draws), value = TRUE)
shape <- rowMeans(f_design_draws[, shape_cols, drop = FALSE], na.rm = TRUE)

# Load Sobol/conduct-stage participant-level effects ---------------------------

f_sobol <- readRDS("processed_data/sobol_lst.rds")

citl_beta <- f_sobol$loghr

# Harmonise older names if present
name_map <- c(
  hdl     = "chol_hdl",
  total   = "chol_total",
  current = "cursmoke",
  t2dur   = "durn",
  n_atc   = "conmeds"
)

names(citl_beta) <- recode(names(citl_beta), !!!name_map)

# Current expected conduct-stage predictor order
Pnm <- names(citl_beta)

# Load typical CITL intercept shift --------------------------------------------

# Prefer current CITL draw object if available. The typical shift is the
# posterior draw-level average across AgD trial arms.
if (file.exists("outputs/conduct/citl/citl_draws.rds")) {
  
  f_citl_draws <- readRDS("outputs/conduct/citl/citl_draws.rds")
  
  if ("delta_mat" %in% names(f_citl_draws)) {
    f_delta <- rowMeans(f_citl_draws$delta_mat, na.rm = TRUE)
  } else if ("delta" %in% names(f_citl_draws)) {
    f_delta <- rowMeans(f_citl_draws$delta, na.rm = TRUE)
  } else {
    stop("Could not find delta_mat or delta in outputs/conduct/citl/citl_draws.rds")
  }
  
} else if (file.exists("outputs/conduct_model/citl/delta_draws.rds")) {
  
  f_delta_draws <- readRDS("outputs/conduct_model/citl/delta_draws.rds")
  f_delta <- rowMeans(f_delta_draws, na.rm = TRUE)
  
} else {
  stop("Could not find CITL delta draws.")
}

# Optionally force a single fixed average shift rather than draw-level typical
# shift. Leave FALSE to retain posterior uncertainty in the typical CITL shift.
use_fixed_mean_delta <- FALSE

if (use_fixed_mean_delta) {
  f_delta <- rep(mean(f_delta, na.rm = TRUE), length(f_delta))
}

# Align posterior draw lengths -------------------------------------------------

n_draw <- min(
  nrow(beta_draws),
  length(mu_mu1_draws),
  length(shape),
  length(f_delta)
)

beta_draws   <- beta_draws[seq_len(n_draw), , drop = FALSE]
mu_mu1_draws <- mu_mu1_draws[seq_len(n_draw)]
shape        <- shape[seq_len(n_draw)]
f_delta      <- f_delta[seq_len(n_draw)]

# Define fixed trial design features ------------------------------------------

# Phase 3, double-blinded, 3 arms (experimental, placebo, active comparator),
# 52 weeks of follow-up, large sample size, one POM, >=100 facilities,
# conducted in the Americas and Europe.
X <- rep(0, ncol(f_design_df$X))
names(X) <- colnames(f_design_df$X)

# These names follow the previous design matrix coding.
# Warnings will be issued if any name has changed.
X <- set_x_value(X, "sample", 2)
X <- set_x_value(X, "blindingdouble", 1)
X <- set_x_value(X, "reg_am", 1)
X <- set_x_value(X, "reg_eur", 1)
X <- set_x_value(X, "n_arms", 2)
X <- set_x_value(X, "n_facilities", 2)

# Build design-stage typical log-rate -----------------------------------------

rate_lp <- as.numeric(beta_draws %*% X)
log_rate_base <- mu_mu1_draws + rate_lp

# Load IPD pooled covariate means/SDs for scaling ------------------------------

cov_scale_raw <- read_csv(
  "vivli/v2/Summaries/covariate_mu_d.csv",
  show_col_types = FALSE
) %>%
  rename_with(tolower)

# Robustly identify variable, mean, and SD columns
var_col <- intersect(names(cov_scale_raw), c("variable", "var", "param", "name", "covariate"))[1]
mu_col  <- intersect(names(cov_scale_raw), c("mean", "mu"))[1]
sd_col  <- intersect(names(cov_scale_raw), c("sd", "d"))[1]

if (any(is.na(c(var_col, mu_col, sd_col)))) {
  stop("Could not identify variable, mean/mu, and sd/d columns in covariate_mu_d.csv")
}

cov_scale <- cov_scale_raw %>%
  transmute(
    param = .data[[var_col]],
    mu = as.numeric(.data[[mu_col]]),
    sd = as.numeric(.data[[sd_col]])
  )

get_scale_row <- function(aliases) {
  out <- cov_scale %>%
    filter(param %in% aliases) %>%
    slice(1)
  
  if (nrow(out) == 0) {
    stop("No scaling row found for aliases: ", paste(aliases, collapse = ", "))
  }
  
  out
}

scale_value <- function(value, aliases) {
  sc <- get_scale_row(aliases)
  (value - sc$mu) / sc$sd
}

scale_age_to_age10 <- function(age_years) {
  # If scale table contains age10, scale age/10.
  # If it contains age, scale age in years and assign to age10 predictor.
  if ("age10" %in% cov_scale$param) {
    scale_value(age_years / 10, "age10")
  } else {
    scale_value(age_years, c("age", "age_years"))
  }
}

# Define hypothetical participant profiles -------------------------------------

# Raw values are stored for supplementary reporting.
# sex is coded as 1 = male, 0 = female.
profile_raw <- tribble(
  ~profile,       ~profile_label,   ~sex_label, ~age, ~sex, ~asian, ~black_aa, ~other, ~cursmoke, ~bmi, ~sbp, ~hba1c, ~chol_hdl, ~chol_total, ~egfr, ~durn, ~comorbs, ~conmeds,
  "low_risk",    "Participant 1",  "Female",   56,   0,    1,      0,         0,      0,         26,   125,  7.2,    1.3,       4.2,         90,    6,     0,       1,
  "low_risk",    "Participant 1",  "Male",     56,   1,    1,      0,         0,      0,         26,   125,  7.2,    1.3,       4.2,         90,    6,     0,       1,
  "med_risk",    "Participant 2",  "Female",   56,   0,    0,      0,         0,      0,         27,   135,  8.5,    1.1,       4.5,         75,    10,    2,       2,
  "med_risk",    "Participant 2",  "Male",     56,   1,    0,      0,         0,      0,         27,   135,  8.5,    1.1,       4.5,         75,    10,    2,       2,
  "higher_risk", "Participant 3",  "Female",   56,   0,    0,      1,         0,      1,         33,   150,  10.0,   0.9,       5.2,         50,    15,    4,       4,
  "higher_risk", "Participant 3",  "Male",     56,   1,    0,      1,         0,      1,         33,   150,  10.0,   0.9,       5.2,         50,    15,    4,       4
)

# Scale continuous covariates to match the current participant-level effects.
# Binary/categorical variables remain coded as 0/1.
profile_scaled <- profile_raw %>%
  rowwise() %>%
  mutate(
    age10 = scale_age_to_age10(age),
    bmi = scale_value(bmi, "bmi"),
    chol_hdl = scale_value(chol_hdl, c("chol_hdl", "hdl")),
    chol_total = scale_value(chol_total, c("chol_total", "total")),
    durn = scale_value(durn, c("durn", "t2dur", "duration", "duration_t2dm")),
    egfr = scale_value(egfr, "egfr"),
    hba1c = scale_value(hba1c, "hba1c"),
    sbp = scale_value(sbp, "sbp"),
    comorbs = comorbs,
    conmeds = conmeds
  ) %>%
  ungroup()

# Check all model predictors are available
missing_profile_cols <- setdiff(Pnm, names(profile_scaled))

if (length(missing_profile_cols) > 0) {
  stop(
    "The following participant-level predictors are missing from profile_scaled: ",
    paste(missing_profile_cols, collapse = ", ")
  )
}

# Save profile and design summaries -------------------------------------------

scenario_design <- tibble(
  design_feature = c(
    "Total follow-up",
    "Treatment type",
    "Eligibility: Lower age limit",
    "Eligibility: Upper age limit",
    "Eligibility: Lower BMI limit",
    "Eligibility: Upper BMI limit",
    "Eligibility: Lower eGFR limit",
    "Eligibility: Upper eGFR limit",
    "Trial size",
    "Trial phase",
    "Trial blinding",
    "Number of arms",
    "Number of facilities",
    "Number of primary outcomes",
    "WHO regions",
    "CITL shift used"
  ),
  specification = c(
    "52 weeks",
    "Placebo and active comparator",
    "18 years",
    "80 years",
    "None",
    "None",
    "None",
    "None",
    "Large (>500 participants)",
    "Phase 3",
    "Double blinding",
    "3 arms (experimental, placebo, active comparator)",
    ">=100 facilities",
    "One primary outcome",
    "Americas and Europe",
    "Average CITL intercept shift across AgD trial arms"
  )
)

write_csv(
  profile_raw,
  "outputs/conduct/scenario/conduct_scenario_profiles.csv"
)

write_csv(
  profile_scaled,
  "outputs/conduct/scenario/conduct_scenario_profiles_scaled.csv"
)

write_csv(
  scenario_design,
  "outputs/conduct/scenario/conduct_scenario_design.csv"
)

# Participant-level linear predictors -----------------------------------------

profile_lp <- profile_scaled %>%
  mutate(
    lp = as.numeric(as.matrix(select(., all_of(Pnm))) %*% citl_beta),
    hazard_multiplier = exp(lp)
  ) %>%
  select(profile, profile_label, sex_label, lp, hazard_multiplier)

profile_lp_summary <- profile_lp %>%
  mutate(
    lp = round(lp, 3),
    hazard_multiplier = round(hazard_multiplier, 3)
  )

write_csv(
  profile_lp_summary,
  "outputs/conduct/scenario/conduct_scenario_linear_predictor_summary.csv"
)

# Typical CITL shift summary ---------------------------------------------------

citl_shift_summary <- summ_vec(f_delta) %>%
  mutate(
    metric = "typical_citl_delta",
    rate_multiplier_mean = mean(exp(f_delta), na.rm = TRUE),
    rate_multiplier_median = median(exp(f_delta), na.rm = TRUE),
    rate_multiplier_q025 = quantile(exp(f_delta), 0.025, na.rm = TRUE),
    rate_multiplier_q975 = quantile(exp(f_delta), 0.975, na.rm = TRUE),
    .before = 1
  )

write_csv(
  citl_shift_summary,
  "outputs/conduct/scenario/conduct_scenario_citl_summary.csv"
)

# Estimate cumulative attrition risk trajectories ------------------------------

t_yrs <- seq(0, 1, length.out = 53)
t_wks <- t_yrs * 52

risk_list <- profile_scaled %>%
  mutate(
    risk = pmap(
      list(profile, profile_label, sex_label, row_number()),
      function(profile, profile_label, sex_label, row_id) {
        
        z <- profile_scaled %>%
          slice(row_id) %>%
          select(all_of(Pnm)) %>%
          as.matrix()
        
        lp <- as.numeric(z %*% citl_beta)
        
        log_rate <- log_rate_base[seq_len(n_draw)] + f_delta + lp
        rate <- exp(log_rate)
        
        risk_mat <- sapply(t_yrs, function(tt) {
          gompertz_risk(rate = rate, shape = shape, time = tt)
        })
        
        tibble(
          week = t_wks,
          risk_q50 = apply(risk_mat, 2, median, na.rm = TRUE),
          risk_q025 = apply(risk_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
          risk_q975 = apply(risk_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
          risk_mean = apply(risk_mat, 2, mean, na.rm = TRUE),
          risk_sd = apply(risk_mat, 2, sd, na.rm = TRUE)
        )
      }
    )
  )

risk_draws <- risk_list %>%
  select(profile, profile_label, sex_label, risk)

saveRDS(
  risk_draws,
  "outputs/conduct/scenario/conduct_scenario_risk_draws.rds"
)

# Risk summary outputs ---------------------------------------------------------

risk_summary <- risk_list %>%
  select(profile, profile_label, sex_label, risk) %>%
  unnest(risk)

risk_summary_out <- risk_summary %>%
  mutate(
    across(starts_with("risk_"), ~ .x * 100)
  )

write_csv(
  risk_summary_out,
  "outputs/conduct/scenario/conduct_scenario_risk_summary.csv"
)

selected_weeks <- c(0, 4, 12, 26, 52)

risk_summary_selected <- risk_summary_out %>%
  mutate(week = round(week)) %>%
  filter(week %in% selected_weeks) %>%
  mutate(
    across(starts_with("risk_"), ~ round(.x, 1)),
    out = paste0(risk_q50, " [", risk_q025, ", ", risk_q975, "]")
  )

write_csv(
  risk_summary_selected,
  "outputs/conduct/scenario/conduct_scenario_risk_summary_selected.csv"
)

# Legacy formatted table for manuscript supplement -----------------------------

risk_summary_tidy <- risk_summary_selected %>%
  select(profile_label, sex_label, week, risk_q50, risk_q025, risk_q975) %>%
  mutate(
    profile_label = factor(
      profile_label,
      levels = c("Participant 1", "Participant 2", "Participant 3")
    ),
    sex_label = factor(sex_label, levels = c("Female", "Male")),
    risk = paste0(risk_q50, " [", risk_q025, ", ", risk_q975, "]")
  ) %>%
  select(profile_label, week, sex_label, risk) %>%
  pivot_wider(
    names_from = sex_label,
    values_from = risk,
    names_prefix = "risk_"
  ) %>%
  arrange(profile_label, week) %>%
  rename(
    profile = profile_label,
    risk_female = risk_Female,
    risk_male = risk_Male
  )

write_csv(
  risk_summary_tidy,
  "outputs/conduct/scenario/conduct_scenario_old.csv"
)

# Plot -------------------------------------------------------------------------

risk_plot <- risk_summary_out %>%
  mutate(
    profile_label = factor(
      profile_label,
      levels = c("Participant 1", "Participant 2", "Participant 3")
    ),
    sex_label = factor(sex_label, levels = c("Female", "Male"))
  ) %>%
  ggplot(aes(x = week, colour = sex_label)) +
  geom_line(aes(y = risk_q50), linewidth = 0.8) +
  geom_line(aes(y = risk_q025), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = risk_q975), linetype = "dashed", linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 80), n.breaks = 5) +
  scale_x_continuous(limits = c(0, 52), n.breaks = 5) +
  facet_wrap(~ profile_label, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Follow-up time (weeks)",
    y = "Predicted attrition risk (%)",
    colour = "Sex"
  )

ggsave(
  "outputs/conduct/scenario/conduct_scenario.png",
  risk_plot,
  height = 4,
  width = 8,
  units = "in"
)

ggsave(
  "outputs/conduct/scenario/conduct_scenario.png",
  risk_plot,
  height = 4,
  width = 8,
  units = "in"
)