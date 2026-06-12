
# Summary statistics and plots
# Conduct-stage model

source("scripts/00_config.R")
library(tidyverse)
library(posterior)
library(bayesplot)
library(flexsurv)
library(metafor)
library(broom)

# Import trial-level Cox model outputs
a_imp_res <- read_csv("vivli/01_res.csv")
a_imp_diag <- read_csv("vivli/02_diag.csv")
a_imp_vcov <- read_csv("vivli/03_cor.csv")

# Tidy covar names
z_covars <- c(
  "Age (in decades)", "BMI (kg/m2)", "HDL cholesterol (mmol/L)",
  "Total cholesterol (mmol/L)", "N comorbidities",
  "Duration of T2DM (years)", "eGFR (ml/min/1.73m2)", "HbA1c (%)",
  "N conc. medication classes", "Asian", "Black/African American",
  "Other race/ethnicity", "Systolic BP (mm/Hg)", "Male", "Current smoker",
  "Treatment: Active comparator"
)

z_covars_lvl <- c(
  "Treatment: Active comparator",
  "Age (in decades)", "Male", "Asian", "Black/African American",
  "Other race/ethnicity", "Current smoker", "BMI (kg/m2)",
  "Systolic BP (mm/Hg)", "HDL cholesterol (mmol/L)",
  "Total cholesterol (mmol/L)",
  "Duration of T2DM (years)", "HbA1c (%)", "eGFR (ml/min/1.73m2)",
  "N conc. medication classes", "N comorbidities"
)

b_covars <- a_imp_res %>% 
  count(term) %>% 
  cbind(z_covars) %>% 
  as_tibble() %>% 
  select(-n) %>% 
  mutate(
    z_covars = factor(
      z_covars,
      levels = z_covars_lvl
    )
  ) %>%
  rename(covars = z_covars) %>%
  arrange(covars)

# IPD model: Trial-level log-hazard ratios -------------------------------------

# Convergences issues with covariates or missing category levels results in
# an incomplete set of covariates across trials
# Addressed by inserting via expansion grid, setting estimates to 0 and
# blown up SEs
# In MA, VCOVs blown up to promote shrinkage towards zero

# Tidy point estimates
b_res_tidy <- b_covars %>% 
  left_join(a_imp_res) %>%
  select(nct_id, covars, estimate, std.error, statistic, p.value) %>% 
  arrange(nct_id, covars)

# Get unique model terms
c_unq <- unique(b_res_tidy$covars)

# Complete terms set (i.e. were missing from IPD, convg issues)
c_covars_grid <- expand.grid(
  nct_id = unique(b_res_tidy$nct_id),
  covars = c_unq
  ) %>% 
  as_tibble()

# Join grid to estimates
# Add zero estimates and large SEs where missing
c_compl <- b_res_tidy %>%
  full_join(c_covars_grid) %>% 
  mutate(
    estimate = if_else(is.na(estimate), 0, estimate),
    std.error = if_else(is.na(std.error), 99999, std.error)
  ) %>% 
  arrange(nct_id, covars)

# Round to 3DP
# Calculate confidence intervals (set to NA where covariate absent from model)
c_rdy <- c_compl %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
    across(lci:uci, ~ if_else(std.error == 99999, NA, std.error)),
    across(estimate:uci, ~ round(., 3)),
  ) %>% 
  select(
    nct_id, covars, estimate, std.error, lci, uci,
    everything()
  )

# Save
write_csv(c_rdy, "outputs/summaries/conduct_model/ipd_cox_model_res.csv")

# Summarise model diagnostics
c_diag <- a_imp_diag
write_csv(c_diag, "outputs/summaries/conduct_model/ipd_cox_model_diag.csv")

# MA: Pooled log-hazard and hazard ratios --------------------------------------

# Correctly order estimates
d_lvl <- a_imp_res %>% count(term) %>% pull(term)

d_est <- c_rdy %>% 
  left_join(b_covars) %>% 
  mutate(term = factor(term, levels = c(d_lvl))) %>% 
  arrange(nct_id, term)

# Get VCOV
d_vcov <- readRDS("processed_data/conduct_cox_vcov.rds")

# Fit meta-regression with random effects for covariates and trials
# Random effects for term and trial
d_ma_id <- d_est$nct_id
d_ma_params <- d_est$term
d_ma_ests <- d_est$estimate
d_ma_vcov <- as.matrix(bdiag(d_vcov$cov_matrix))

d_ma <- rma.mv(
  yi = d_ma_ests,
  V = d_ma_vcov,
  mods = ~ factor(d_ma_params) - 1,
  random = ~ d_ma_params | d_ma_id
)

# Checks
summary(d_ma)
d_ma$tau2
d_ma$rho

# Tidy results
d_ma_tidy <- tidy(d_ma) %>% 
  mutate(term = gsub("factor\\(d_ma_params\\)", "", term)) %>% 
  mutate(
    lo = estimate - 1.96 * std.error,
    hi = estimate + 1.96 * std.error,
    hr = exp(estimate),
    lo_hr = exp(lo),
    hi_hr = exp(hi)
  ) %>% 
  select(-type) %>% 
  left_join(b_covars) %>% 
  arrange(covars) %>% 
  select(
    covars, estimate, lo, hi, hr, lo_hr, hi_hr,
    std.error, statistic, p.value
  ) %>% 
  rename(
    loghr = estimate,
    loghr_lo = lo,
    loghr_hi = hi
  ) %>% 
  mutate(
    across(loghr:statistic, ~ round(., 2)),
    p.value = round(p.value, 3)
  )

# Plot
d_ma_plot <- d_ma_tidy %>% 
  ggplot(aes(x = hr, xmin = lo_hr, xmax = hi_hr, y = fct_rev(covars))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0.4, 1.6), n.breaks = 7) +
  theme_classic() +
  labs(
    x = "Mean effect on attrition risk (95% confidence intervals)", 
    y = NULL
  )

# Save
write_csv(d_ma_tidy, "outputs/summaries/conduct_model/ipd_ma_estimates.csv")

ggsave(
  "outputs/plots/conduct_ipd_chars_risk.png",
  d_ma_plot,
  height = 6,
  width = 6,
  units = "in"
)

# Calibration-in-the-large -----------------------------------------------------

# Import materials from CITL
e_citl_draws <- readRDS("outputs/conduct_model/citl/delta_draws.rds")

# Import CITL summary
# Add trial ID
e_citl_sum <- read_csv("outputs/conduct_model/citl/delta_risk_sum.csv") %>% 
  cbind(
    tibble(
      trial_id = readRDS("processed_data/latest_dataset.rds") %>% 
        pull(trial_id)
    )
  ) %>% 
  as_tibble()

# Summarise distribution of delta (intercept shifts)
e_citl_overall <- tibble(
  mean   = mean(e_citl_draws),
  median = median(e_citl_draws),
  sd     = sd(e_citl_draws),
  lo     = quantile(e_citl_draws, 0.025),
  hi     = quantile(e_citl_draws, 0.975)
)

# Summarise trial-level delta draws
e_citl_trial <- tibble(
  trial_id = e_citl_sum$trial_id,
  trial_arm = seq_len(ncol(e_citl_draws)),
  delta_mean = colMeans(e_citl_draws),
  delta_median = apply(e_citl_draws, 2, median),
  delta_lo = apply(e_citl_draws, 2, quantile, 0.025),
  delta_hi = apply(e_citl_draws, 2, quantile, 0.975)
)

# Summarise distribution of calibration shifts
e_citl_dist <- tibble(delta = as.vector(e_citl_draws)) %>%
  ggplot(aes(x = delta)) +
  geom_density(fill = "steelblue", alpha = .4) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  labs(
    x = "Calibration intercept shift (δ)",
    y = "Density"
  )

# Summarise pre and post-CITL predicted risk
e_citl_risk <- e_citl_sum %>% 
  reframe(
    obs_mean  = mean(p_obs),
    pred_pre_mean  = mean(p_pred),
    pred_post_mean = mean(p_citl)
  )

# Summarise absolute calibration error
e_citl_error <- e_citl_sum %>%
  mutate(
    err_pre  = abs(p_pred - p_obs),
    err_post = abs(p_citl - p_obs)
  ) %>%
  reframe(
    mae_pre  = mean(err_pre),
    mae_post = mean(err_post)
  )

# Calibration plot
e_citl_plot_dat <- tibble(
  arm = rep(e_citl_sum$arm, 2),
  obs = rep(e_citl_sum$p_obs, 2),
  pred = c(e_citl_sum$p_pred, e_citl_sum$p_citl),
  stage = rep(c("Pre-CITL", "Post-CITL"), each = nrow(e_citl_sum))
)

e_citl_plot <- e_citl_plot_dat %>%
  ggplot(aes(x = obs, y = pred, colour = stage)) +
  geom_abline(linetype = 2) +
  geom_point(size = 0.5) +
  theme_classic() +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk",
    colour = NULL
  ) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5)

# Save summaries
write_csv(e_citl_overall, "outputs/summaries/conduct_model/citl_delta_sum_overall.csv")
write_csv(e_citl_sum, "outputs/summaries/conduct_model/citl_delta_sum_trial.csv")
write_csv(e_citl_risk, "outputs/summaries/conduct_model/citl_risk_sum.csv")
write_csv(e_citl_error, "outputs/summaries/conduct_model/citl_error.csv")

write_csv(e_citl_trial, "outputs/summaries/conduct_model/citl_tmp.csv")

# Save plots
ggsave(
  "outputs/plots/conduct_citl_delta.png",
  e_citl_dist,
  width = 4,
  height = 4,
  units = "in"
)

ggsave(
  "outputs/plots/conduct_citl_calibration.png",
  e_citl_plot,
  width = 4,
  height = 4,
  units = "in"
)

# Case scenario attrition risk prediction --------------------------------------

# Load design model draws
f_design_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")
f_design_df <- readRDS("mars/design_model/latest_stan_design.rds")

# f_design_df$X %>% 
#   as_tibble() %>% 
#   pivot_longer(everything()) %>% 
#   group_by(name, value) %>% 
#   distinct() %>% 
#   arrange(name, value) %>% 
#   print(n = 51)

# Get posterior draws
f_design_draws <- as_draws_matrix(as.array(f_design_mdl))

# Extract parameters
beta_cols   <- grep("^beta\\[", colnames(f_design_draws), value = TRUE)
mu_mu1_col  <- "mu_mu1"
tau_mu1_col <- "tau_mu1"

beta_draws    <- f_design_draws[, beta_cols]
mu_mu1_draws  <- f_design_draws[, mu_mu1_col]
tau_mu1_draws <- f_design_draws[, tau_mu1_col]

# Load CITL outputs (intercept shifts and Sobol-sequenced characteristics)
f_delta_draws <- readRDS("outputs/conduct_model/citl/delta_draws.rds")
f_sobol <- readRDS("processed_data/sobol_lst.rds")

# Extract log-rate and shape
f_log_b <- posterior::as_draws_matrix(f_design_draws[, grep("log_b", variables(f_design_draws))])
f_a <- posterior::as_draws_matrix(f_design_draws[, grep("^a\\[", variables(f_design_draws))])

# Collapse log-rate and shape across trials
f_log_b <- rowMeans(f_log_b)
f_a <- rowMeans(f_a)

# Intercept shift (delta)
f_delta <- rowMeans(f_delta_draws)

# Define trial design features
# Phase 3, double-blinded, 3 arms (one placebo, one active comparator)
# 52-weeks of follow-up, 500 participants
# Age limit <=80 years, no BMI or eGFR limits
# Situated in regions of the Americas and Europe
# One POM, 100 sites/facilities
X <- rep(0, ncol(f_design_df$X))
names(X) <- colnames(f_design_df$X)

# All set to 0 by default, change values of targets where necessary
X["sample"] <- 2
X["blindingdouble"] <- 1
X["reg_am"] <- 1
X["reg_eur"] <- 1
X["n_arms"] <- 2
X["n_facilities"] <- 2

# Build baseline risk
rate_lp <- beta_draws %*% X
mu1 <- mu_mu1_draws
log_rate_base <- mu1 + rate_lp
rate_base <- exp(log_rate_base)
shape <- rowMeans(f_design_draws[, grepl("^a\\[", colnames(f_design_draws))])

# Create participant profiles (low to high risk)
# Female and male participants, age and race fixed for each group
f_profile1 <- tibble(
  profile = "low_risk",
  param = names(f_sobol$loghr),
  sex = "female",
  value = c(
    age10 = 56 / 10,
    bmi = 26 / 5.57,
    hdl = 1.3,
    total = 4.2,
    current = 0,
    t2dur = 6 / 6.14,
    egfr = 90 / 20.59,
    hba1c = 7.2,
    asian = 1,
    black_aa = 0,
    other = 0,
    sbp = 125 / 15.94,
    sex = 0,
    comorbs = 0,
    n_atc = 1
  )
)

f_profile2 <- tibble(
  profile = "low_risk",
  param = names(f_sobol$loghr),
  sex = "male",
  value = c(
    age10 = 56 / 10,
    bmi = 26 / 5.57,
    hdl = 1.3,
    total = 4.2,
    current = 0,
    t2dur = 6 / 6.14,
    egfr = 90 / 20.59,
    hba1c = 7.2,
    asian = 1,
    black_aa = 0,
    other = 0,
    sbp = 125 / 15.94,
    sex = 1,
    comorbs = 0,
    n_atc = 1
  )
)

f_profile3 <- tibble(
  profile = "med_risk",
  param = names(f_sobol$loghr),
  sex = "female",
  value = c(
    age10 = 56 / 10,
    bmi = 27 / 5.57,
    hdl = 1.1,
    total = 4.5,
    current = 0,
    t2dur = 10,
    egfr = 75 / 20.59,
    hba1c = 8.5,
    asian = 0,
    black_aa = 0,
    other = 0,
    sbp = 135 / 15.94,
    sex = 0,
    comorbs = 2,
    n_atc = 2
  )
)

f_profile4 <- tibble(
  profile = "med_risk",
  param = names(f_sobol$loghr),
  sex = "male",
  value = c(
    age10 = 56 / 10,
    bmi = 27 / 5.57,
    hdl = 1.1,
    total = 4.5,
    current = 0,
    t2dur = 10,
    egfr = 75 / 20.59,
    hba1c = 8.5,
    asian = 0,
    black_aa = 0,
    other = 0,
    sbp = 135 / 15.94,
    sex = 1,
    comorbs = 2,
    n_atc = 2
  )
)

f_profile5 <- tibble(
  profile = "higher_risk",
  param = names(f_sobol$loghr),
  sex = "female",
  value = c(
    age10 = 56 / 10,
    bmi = 33 / 5.57,
    hdl = 0.9,
    total = 5.2,
    current = 1,
    t2dur = 15,
    egfr = 50 / 20.59,
    hba1c = 10,
    asian = 0,
    black_aa = 1,
    other = 0,
    sbp = 150 / 15.94,
    sex = 0,
    comorbs = 4,
    n_atc = 4
  )
)

f_profile6 <- tibble(
  profile = "higher_risk",
  param = names(f_sobol$loghr),
  sex = "male",
  value = c(
    age10 = 56 / 10,
    bmi = 33 / 5.57,
    hdl = 0.9,
    total = 5.2,
    current = 1,
    t2dur = 15,
    egfr = 50 / 20.59,
    hba1c = 10,
    asian = 0,
    black_aa = 1,
    other = 0,
    sbp = 150 / 15.94,
    sex = 1,
    comorbs = 4,
    n_atc = 4
  )
)

# Nest to automate estimation
f_join <- reduce(mget(ls(pattern = "f_profile")), bind_rows) %>% 
  group_by(profile, sex) %>% 
  nest()

# Compute linear predictors
# Build risk parameters and estimate risk
t_yrs <- seq(0, 1, length.out = 53)
t_wks <- t_yrs * 52

f_risk <- f_join %>% 
  mutate(
    risk = map(data, ~ {
      lp <- .x %>% 
        mutate(loghr = f_sobol$loghr[param]) %>%
        summarise(lp = sum(value * loghr)) %>%
        pull(lp)
      
      log_rate <- log_rate_base + f_delta + lp
      rate <- exp(log_rate)
      
      risk <- sapply(t_yrs, function(tt) {
        pgompertz(tt, shape = shape, rate = rate)
      })
    })
  )

# Unnest and summarise
f_risk_sum <- f_risk %>% 
  mutate(
    risk_sum = map(risk, ~ {
      tibble(
        week = t_wks,
        risk_q50 = apply(.x, 2, median),
        risk_q025 = apply(.x, 2, quantile, 0.025),
        risk_q975 = apply(.x, 2, quantile, 0.975) 
      )
    })
  )

# Plot
f_risk_plot <- f_risk_sum %>% 
  select(profile, sex, risk_sum) %>% 
  unnest(risk_sum) %>% 
  mutate(
    profile = case_match(
      profile,
      "low_risk" ~ "Participant 1",
      "med_risk" ~ "Participant 2",
      "higher_risk" ~ "Participant 3"
    ),
    profile = factor(
      profile,
      levels = c("Participant 1", "Participant 2", "Participant 3")
    ),
    sex = case_match(
      sex,
      "female" ~ "Female",
      "male" ~ "Male"
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    across(starts_with("risk_"), ~ . * 100)
  ) %>% 
  ggplot(aes(x = week, y = risk_q50, colour = sex)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = risk_q025, ymax = risk_q975, colour = sex),
    alpha = 0,
    fill = "steelblue",
    linetype = "dashed"
  ) +
  scale_y_continuous(limits = c(0, 100), n.breaks = 5) +
  scale_x_continuous(limits = c(0, 52), n.breaks = 5) + 
  facet_wrap(~profile, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Follow-up time (weeks)",
    y = "Predicted attrition risk (%)",
    colour = "Sex"
  )

ggsave(
  "outputs/plots/conduct_scenario.png",
  f_risk_plot,
  height = 4,
  width = 8,
  units = "in"
)

# Save risk summary
f_risk_sum_tidy <- f_risk_sum %>% 
  select(profile, sex, risk_sum) %>% 
  unnest(risk_sum) %>% 
  mutate(
    week = round(week),
    profile = case_match(
      profile,
      "low_risk" ~ "Participant 1",
      "med_risk" ~ "Participant 2",
      "higher_risk" ~ "Participant 3"
    ),
    profile = factor(
      profile,
      levels = c("Participant 1", "Participant 2", "Participant 3")
    ),
    sex = case_match(
      sex,
      "female" ~ "Female",
      "male" ~ "Male"
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    across(contains("risk_"), ~ . * 100),
    across(contains("risk_"), ~ round(., 1))
  ) %>% 
  pivot_wider(
    names_from = "sex", values_from = c("risk_q50", "risk_q025", "risk_q975")
  ) %>% 
  mutate(
    risk_female = paste0(
      risk_q50_Female, " [", risk_q025_Female, ", ", risk_q975_Female, "]"
    ),
    risk_male = paste0(
      risk_q50_Male, " [", risk_q025_Male, ", ", risk_q975_Male, "]"
    )
  ) %>% 
  select(profile, week, risk_female, risk_male)

write_csv(f_risk_sum_tidy, "outputs/summaries/conduct_model/conduct_scenario.csv")
