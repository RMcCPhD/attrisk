
rm(list = ls())
gc()

source("scripts/00_config.R")
library(tidyverse)
library(posterior)

# Output directory
z_out <- "outputs/conduct_stage/profile_risks/"
dir.create(z_out, recursive = TRUE, showWarnings = FALSE)

# Parameters -------------------------------------------------------------------

# Select representative profiles around the 10th, 50th and 90th percentiles of
# the participant prognostic score excluding sex
target_probs <- c(0.10, 0.50, 0.90)
selection_band <- 0.01

# Retain posterior uncertainty in the average CITL intercept shift
use_fixed_mean_delta <- FALSE

# Helper functions -------------------------------------------------------------

summ_vec <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

gompertz_risk <- function(rate, shape, time) {
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

f_design_mdl <- readRDS("outputs/design_stage/model/design_model.rds")
f_design_df <- readRDS("processed_data/stan_data.rds")
f_design_draws <- posterior::as_draws_matrix(as.array(f_design_mdl))
beta_cols <- grep("^beta\\[", colnames(f_design_draws), value = TRUE)
beta_draws <- f_design_draws[, beta_cols, drop = FALSE]
mu_mu1_draws <- f_design_draws[, "mu_mu1"]

# Use the average shape across AgD trial arms for the fixed illustrative design
shape_cols <- grep("^a\\[", colnames(f_design_draws), value = TRUE)
shape <- rowMeans(f_design_draws[, shape_cols, drop = FALSE], na.rm = TRUE)

# Load participant-level effects and generated profiles ------------------------

# Import Sobol inputs and nested profiles generated across all completed AgD
# datasets by 08c_sobol_generate_mi.R
sobol_mi <- readRDS("outputs/conduct_stage/profiles_imputation.rds")
sobol_nested_mi <- readRDS("outputs/conduct_stage/profiles_imputation_unnest.rds")

stopifnot(
  all(c(".imp", "sobol_lst") %in% names(sobol_mi)),
  all(c(".imp", "sobol_unnest") %in% names(sobol_nested_mi)),
  !anyDuplicated(sobol_mi$.imp),
  !anyDuplicated(sobol_nested_mi$.imp),
  identical(sobol_mi$.imp, sobol_nested_mi$.imp)
)

# Participant-level effects are fixed across AgD imputations. Verify that the
# coefficient vectors agree before using the common vector downstream.
citl_beta_mi <- map(sobol_mi$sobol_lst, "loghr")
citl_beta <- citl_beta_mi[[1]]

stopifnot(
  length(citl_beta_mi) == nrow(sobol_mi),
  all(map_lgl(citl_beta_mi, ~ identical(names(.x), names(citl_beta)))),
  all(map_lgl(
    citl_beta_mi,
    ~ isTRUE(all.equal(
      unname(.x),
      unname(citl_beta),
      tolerance = 1e-12
    ))
  ))
)

# Harmonise older names if present
name_map <- c(
  hdl = "chol_hdl",
  total = "chol_total",
  current = "cursmoke",
  t2dur = "durn",
  n_atc = "conmeds"
)

names(citl_beta) <- recode(names(citl_beta), !!!name_map)
Pnm <- names(citl_beta)

if (!"sex" %in% Pnm) {
  stop("The participant-level coefficient vector does not contain sex.")
}

nonsex_terms <- setdiff(Pnm, "sex")

# Combine one nested profile column across imputations and control arms. Each
# imputation-arm combination contributes the same number of Sobol profiles, as
# in the pooled-profile CITL exercise performed in 09ax.
bind_profile_sets <- function(profile_column) {
  map2_dfr(
    sobol_nested_mi$.imp,
    sobol_nested_mi$sobol_unnest,
    function(.imp, profile_data) {
      imp_id <- .imp
      profile_data %>%
        transmute(
          trial_id,
          profiles = .data[[profile_column]]
        ) %>%
        unnest(cols = profiles) %>%
        group_by(trial_id) %>%
        mutate(source_profile = row_number(), .after = trial_id) %>%
        ungroup() %>%
        mutate(source_imputation = imp_id, .before = 1)
    }
  )
}

# Reconstruct raw and scaled Sobol profiles across all completed datasets. Use
# explicit profile keys to guarantee that raw and scaled records remain aligned.
profile_keys <- c("source_imputation", "trial_id", "source_profile")

sobol_raw <- bind_profile_sets("sobol_unnest") %>%
  arrange(source_imputation, trial_id, source_profile) %>%
  mutate(candidate_id = row_number(), .before = 1)

sobol_scaled <- bind_profile_sets("sobol_rdy") %>%
  left_join(
    sobol_raw %>% select(candidate_id, all_of(profile_keys)),
    by = profile_keys
  ) %>%
  relocate(candidate_id) %>%
  arrange(candidate_id)

missing_raw_terms <- setdiff(
  c(
    "age", "male_p", "bmi", "durn", "sbp", "chol_hdl", "chol_total",
    "hba1c", "egfr", "cursmoke", "comorbs", "conmeds", "white",
    "black_aa", "asian", "other"
  ),
  names(sobol_raw)
)
missing_scaled_terms <- setdiff(Pnm, names(sobol_scaled))

if (length(missing_raw_terms) > 0 || length(missing_scaled_terms) > 0) {
  stop(
    "Generated Sobol profiles are missing required variables: ",
    paste(c(missing_raw_terms, missing_scaled_terms), collapse = ", ")
  )
}

stopifnot(
  nrow(sobol_raw) == nrow(sobol_scaled),
  n_distinct(sobol_raw$source_imputation) == nrow(sobol_nested_mi),
  n_distinct(sobol_scaled$source_imputation) == nrow(sobol_nested_mi),
  !anyDuplicated(sobol_raw[profile_keys]),
  !anyDuplicated(sobol_scaled[profile_keys]),
  identical(sobol_raw$candidate_id, sobol_scaled$candidate_id),
  identical(
    sobol_raw$source_imputation,
    sobol_scaled$source_imputation
  ),
  identical(sobol_raw$trial_id, sobol_scaled$trial_id),
  identical(sobol_raw$source_profile, sobol_scaled$source_profile),
  all(sobol_raw$male_p == sobol_scaled$sex)
)

# Calculate prognostic scores excluding sex -----------------------------------

nonsex_matrix <- sobol_scaled %>% select(all_of(nonsex_terms)) %>% as.matrix()

sobol_scaled <- sobol_scaled %>%
  mutate(
    lp_nonsex = as.numeric(nonsex_matrix %*% citl_beta[nonsex_terms])
  )

# Standardise candidate covariates before calculating medoid distances
medoid_centre <- colMeans(nonsex_matrix)
medoid_scale <- apply(nonsex_matrix, 2, sd)

if (any(!is.finite(medoid_scale)) || any(medoid_scale <= 0)) {
  stop("Cannot standardise one or more covariates for medoid selection.")
}

medoid_matrix <- sweep(
  sweep(nonsex_matrix, 2, medoid_centre, FUN = "-"),
  2,
  medoid_scale,
  FUN = "/"
)

score_ecdf <- ecdf(sobol_scaled$lp_nonsex)

# Select the medoid among profiles within one percentile point of each target
select_representative_profile <- function(target_probability) {
  
  lower_probability <- max(0, target_probability - selection_band)
  upper_probability <- min(1, target_probability + selection_band)
  
  target_score <- quantile(
    sobol_scaled$lp_nonsex,
    target_probability,
    names = FALSE
  )
  lower_score <- quantile(
    sobol_scaled$lp_nonsex,
    lower_probability,
    names = FALSE
  )
  upper_score <- quantile(
    sobol_scaled$lp_nonsex,
    upper_probability,
    names = FALSE
  )
  
  candidate_idx <- which(
    sobol_scaled$lp_nonsex >= lower_score &
      sobol_scaled$lp_nonsex <= upper_score
  )
  
  if (length(candidate_idx) == 0) {
    stop("No Sobol profiles were found around percentile ", target_probability)
  }
  
  candidate_distance <- as.matrix(
    dist(medoid_matrix[candidate_idx, , drop = FALSE])
  )
  medoid_position <- which.min(rowSums(candidate_distance))
  selected_idx <- candidate_idx[medoid_position]
  
  tibble(
    target_probability = target_probability,
    lower_probability = lower_probability,
    upper_probability = upper_probability,
    target_score = target_score,
    selected_score = sobol_scaled$lp_nonsex[selected_idx],
    achieved_percentile = score_ecdf(sobol_scaled$lp_nonsex[selected_idx]),
    absolute_score_difference = abs(
      sobol_scaled$lp_nonsex[selected_idx] - target_score
    ),
    n_candidate_profiles = length(candidate_idx),
    mean_medoid_distance = mean(candidate_distance[medoid_position, ]),
    candidate_id = sobol_scaled$candidate_id[selected_idx],
    source_imputation = sobol_scaled$source_imputation[selected_idx],
    source_trial_id = sobol_scaled$trial_id[selected_idx],
    source_profile = sobol_scaled$source_profile[selected_idx],
    source_sex = sobol_scaled$sex[selected_idx]
  )
}

profile_definitions <- tibble(
  target_probability = target_probs,
  profile = c("p10", "p50", "p90"),
  profile_label = c(
    "Lower-risk profile",
    "Moderate-risk profile",
    "Higher-risk profile"
  )
)

profile_selection <- map_dfr(
  target_probs,
  select_representative_profile
) %>%
  left_join(profile_definitions, by = "target_probability") %>%
  relocate(profile, profile_label)

stopifnot(
  profile_selection$achieved_percentile >=
    profile_selection$lower_probability,
  profile_selection$achieved_percentile <=
    profile_selection$upper_probability
)

# Join selected rows back to raw and scaled characteristics --------------------

selected_raw <- profile_selection %>%
  left_join(
    sobol_raw %>%
      select(-source_imputation, -trial_id, -source_profile),
    by = "candidate_id"
  )

selected_scaled <- profile_selection %>%
  left_join(
    sobol_scaled %>%
      select(candidate_id, all_of(nonsex_terms)),
    by = "candidate_id"
  )

# Duplicate each selected profile as female and male while holding all other
# participant characteristics constant
profile_raw <- bind_rows(
  selected_raw %>% mutate(sex_label = "Female", sex = 0),
  selected_raw %>% mutate(sex_label = "Male", sex = 1)
) %>%
  mutate(
    target_percentile = 100 * target_probability,
    achieved_percentile = 100 * achieved_percentile,
    race = case_when(
      white == 1 ~ "White",
      black_aa == 1 ~ "Black or African American",
      asian == 1 ~ "Asian",
      other == 1 ~ "Other"
    )
  ) %>%
  arrange(target_probability, sex) %>%
  select(
    profile, profile_label, 
    target_percentile, selected_score, achieved_percentile, 
    sex_label, sex, age, 
    race, white, black_aa, asian, other,
    cursmoke, bmi, sbp, hba1c, chol_hdl, chol_total, egfr, durn,
    comorbs, conmeds,
    source_trial_id, source_imputation, source_profile, source_sex
  )

profile_scaled <- bind_rows(
  selected_scaled %>% mutate(sex_label = "Female", sex = 0),
  selected_scaled %>% mutate(sex_label = "Male", sex = 1)
) %>%
  mutate(
    target_percentile = 100 * target_probability,
    achieved_percentile = 100 * achieved_percentile
  ) %>%
  arrange(target_probability, sex) %>%
  select(
    profile,
    profile_label,
    target_percentile,
    selected_score,
    achieved_percentile,
    sex_label,
    sex,
    all_of(nonsex_terms),
    source_trial_id,
    source_imputation,
    source_profile
  )

stopifnot(
  nrow(profile_raw) == 2L * length(target_probs),
  nrow(profile_scaled) == 2L * length(target_probs),
  !anyNA(profile_raw),
  !anyNA(profile_scaled),
  all(rowSums(as.matrix(
    profile_raw %>% select(white, black_aa, asian, other)
  )) == 1)
)

# Load CITL intercept shift ----------------------------------------------------

# Import the primary pooled-imputation CITL output promoted by 09ax. The delta
# matrix has posterior draws in rows and AgD control arms in columns; average
# across arms within each draw for the fixed illustrative trial design.
f_citl_draws <- readRDS("outputs/conduct_stage/citl/citl_draws.rds")

required_citl_objects <- c(
  "I", "J", "J_by_imputation", 
  "M", "imputation", "S", "trial_id",
  "delta"
)

missing_citl_objects <- setdiff(required_citl_objects, names(f_citl_draws))

if (length(missing_citl_objects) > 0) {
  stop(
    "Pooled-imputation CITL output is missing: ",
    paste(missing_citl_objects, collapse = ", ")
  )
}

delta_mat <- f_citl_draws$delta

profile_count_check <- sobol_scaled %>%
  count(source_imputation, trial_id, name = "n_profiles") %>%
  mutate(
    expected_profiles = f_citl_draws$J_by_imputation[
      match(source_imputation, f_citl_draws$imputation)
    ]
  )

stopifnot(
  is.matrix(delta_mat),
  all(dim(delta_mat) == c(f_citl_draws$S, f_citl_draws$I)),
  f_citl_draws$M == nrow(sobol_mi),
  identical(f_citl_draws$imputation, sobol_mi$.imp),
  setequal(f_citl_draws$trial_id, sobol_scaled$trial_id),
  length(f_citl_draws$J_by_imputation) == f_citl_draws$M,
  f_citl_draws$J == sum(f_citl_draws$J_by_imputation),
  nrow(sobol_scaled) == f_citl_draws$J * f_citl_draws$I,
  !anyNA(profile_count_check$expected_profiles),
  all(
    profile_count_check$n_profiles ==
      profile_count_check$expected_profiles
  ),
  all(is.finite(delta_mat))
)

f_delta <- rowMeans(delta_mat, na.rm = TRUE)

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

beta_draws <- beta_draws[seq_len(n_draw), , drop = FALSE]
mu_mu1_draws <- mu_mu1_draws[seq_len(n_draw)]
shape <- shape[seq_len(n_draw)]
f_delta <- f_delta[seq_len(n_draw)]

# Define fixed trial design features ------------------------------------------

# Phase 3, double-blinded, 3 arms (experimental, placebo, active comparator),
# 52 weeks of follow-up, large sample size, one POM, >=100 facilities,
# conducted in the Americas and Europe
X <- rep(0, ncol(f_design_df$X))
names(X) <- colnames(f_design_df$X)

X <- set_x_value(X, "n_cat500-749", 1)
X <- set_x_value(X, "maskingdouble", 1)
X <- set_x_value(X, "region_of_the_americas", 1)
X <- set_x_value(X, "european_region", 1)
X <- set_x_value(X, "n_arms_cat", 0)
X <- set_x_value(X, "n_sites_cat>=200", 1)

rate_lp <- as.numeric(beta_draws %*% X)
log_rate_base <- mu_mu1_draws + rate_lp

# Save profile-selection and design summaries ---------------------------------

profile_selection_out <- profile_selection %>%
  mutate(
    across(ends_with("probability"), ~ 100 * .x),
    achieved_percentile = 100 * achieved_percentile
  )

scenario_design <- tibble(
  design_feature = c(
    "Total follow-up",
    "Treatment type",
    "Eligibility: Lower age limit", "Eligibility: Upper age limit",
    "Eligibility: Lower BMI limit", "Eligibility: Upper BMI limit",
    "Eligibility: Lower eGFR limit", "Eligibility: Upper eGFR limit",
    "Trial size", "Trial phase", "Trial blinding", "Number of arms",
    "Number of facilities", "Number of primary outcomes",
    "WHO regions", "Participant profiles", "CITL shift used"
  ),
  specification = c(
    "52 weeks", "Placebo and active comparator", "18 years", "80 years",
    "None", "None", "None", "None", "Large (>500 participants)",
    "Phase 3", "Double blinding",
    "3 arms (experimental, placebo, active comparator)",
    ">=100 facilities", "One primary outcome", "Americas and Europe",
    paste0(
      "Sobol profiles pooled across ",
      nrow(sobol_mi),
      " completed AgD datasets; medoids selected within +/-",
      100 * selection_band,
      " percentile point of the 10th, 50th and 90th percentiles of the ",
      "non-sex prognostic score"
    ),
    paste0(
      "Average pooled-imputation CITL intercept shift across ",
      f_citl_draws$I,
      " AgD trial arms"
    )
  )
)

write_csv(
  profile_selection_out,
  file.path(z_out, "conduct_scenario_profile_selection.csv")
)

write_csv(
  profile_raw,
  file.path(z_out, "conduct_scenario_profiles.csv")
)

write_csv(
  profile_scaled,
  file.path(z_out, "conduct_scenario_profiles_scaled.csv")
)

write_csv(
  scenario_design,
  file.path(z_out, "conduct_scenario_design.csv")
)

# Plot the prognostic-score distribution and selected profiles
profile_selection_plot <- sobol_scaled %>%
  ggplot(aes(x = lp_nonsex)) +
  geom_density(fill = "grey80", colour = "grey40") +
  geom_vline(
    data = profile_selection,
    aes(xintercept = selected_score, colour = profile_label),
    linewidth = 0.8
  ) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(
    x = "Participant prognostic score excluding sex",
    y = "Density",
    colour = NULL
  )

ggsave(
  file.path(z_out, "conduct_scenario_profile_selection.png"),
  profile_selection_plot,
  width = 8,
  height = 4,
  units = "in"
)

# Participant-level linear predictors -----------------------------------------

profile_lp <- profile_scaled %>%
  mutate(
    lp_nonsex = as.numeric(
      as.matrix(select(., all_of(nonsex_terms))) %*%
        citl_beta[nonsex_terms]
    ),
    lp = as.numeric(as.matrix(select(., all_of(Pnm))) %*% citl_beta),
    hazard_multiplier = exp(lp)
  ) %>%
  select(
    profile,
    profile_label,
    target_percentile,
    sex_label,
    selected_score,
    lp_nonsex,
    lp,
    hazard_multiplier
  )

stopifnot(
  all(abs(profile_lp$selected_score - profile_lp$lp_nonsex) < 1e-10)
)

profile_lp_summary <- profile_lp %>%
  mutate(across(c(lp_nonsex, lp, hazard_multiplier), ~ round(.x, 3)))

write_csv(
  profile_lp_summary,
  file.path(z_out, "conduct_scenario_linear_predictor_summary.csv")
)

# CITL shift summary -----------------------------------------------------------

citl_shift_summary <- summ_vec(f_delta) %>%
  mutate(
    metric = "average_citl_delta",
    rate_multiplier_mean = mean(exp(f_delta), na.rm = TRUE),
    rate_multiplier_median = median(exp(f_delta), na.rm = TRUE),
    rate_multiplier_q025 = quantile(
      exp(f_delta),
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    rate_multiplier_q975 = quantile(
      exp(f_delta),
      0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .before = 1
  )

write_csv(
  citl_shift_summary,
  file.path(z_out, "conduct_scenario_citl_summary.csv")
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
        log_rate <- log_rate_base + f_delta + lp
        rate <- exp(log_rate)
        
        risk_mat <- sapply(t_yrs, function(tt) {
          gompertz_risk(rate = rate, shape = shape, time = tt)
        })
        
        tibble(
          week = t_wks,
          risk_q50 = apply(risk_mat, 2, median, na.rm = TRUE),
          risk_q025 = apply(
            risk_mat,
            2,
            quantile,
            probs = 0.025,
            na.rm = TRUE,
            names = FALSE
          ),
          risk_q975 = apply(
            risk_mat,
            2,
            quantile,
            probs = 0.975,
            na.rm = TRUE,
            names = FALSE
          ),
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
  file.path(z_out, "conduct_scenario_risk_draws.rds")
)

# Risk summary outputs ---------------------------------------------------------

risk_summary <- risk_list %>%
  select(profile, profile_label, target_percentile, sex_label, risk) %>%
  unnest(risk)

risk_summary_out <- risk_summary %>%
  mutate(across(starts_with("risk_"), ~ 100 * .x))

write_csv(
  risk_summary_out,
  file.path(z_out, "conduct_scenario_risk_summary.csv")
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
  file.path(z_out, "conduct_scenario_risk_summary_selected.csv")
)

# Formatted table for manuscript supplement -----------------------------------

risk_summary_tidy <- risk_summary_selected %>%
  select(profile_label, target_percentile, sex_label, week, out) %>%
  mutate(
    profile_label = factor(
      profile_label,
      levels = profile_definitions$profile_label
    ),
    sex_label = factor(sex_label, levels = c("Female", "Male"))
  ) %>%
  pivot_wider(
    names_from = sex_label,
    values_from = out,
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
  file.path(z_out, "conduct_scenario_table.csv")
)

# Plot trajectories ------------------------------------------------------------

risk_plot <- risk_summary_out %>%
  mutate(
    profile_label = factor(
      profile_label,
      levels = profile_definitions$profile_label
    ),
    sex_label = factor(sex_label, levels = c("Female", "Male"))
  ) %>%
  ggplot(aes(x = week, colour = sex_label)) +
  geom_line(aes(y = risk_q50), linewidth = 0.8) +
  geom_line(aes(y = risk_q025), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = risk_q975), linetype = "dashed", linewidth = 0.6) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 5,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(limits = c(0, 52), n.breaks = 5) +
  facet_wrap(~profile_label) +
  theme_classic() +
  labs(
    x = "Follow-up time (weeks)",
    y = "Predicted attrition risk (%)",
    colour = "Sex"
  )

ggsave(
  file.path(z_out, "conduct_scenario.png"),
  risk_plot,
  height = 4,
  width = 10,
  units = "in"
)