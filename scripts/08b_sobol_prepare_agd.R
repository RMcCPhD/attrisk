
# Prepare baseline characteristics for aggregate-level data
# Pulled in from previous project covariate extraction

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import tidied AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import IPD baseline characteristics
# Used as target schema and to fill missing AgD summaries where possible
imp_bl_ipd <- read_csv("outputs/manuscript_summaries/ipd_bl_sum.csv") %>% 
  rename(trial_id = nct_id)

ipd_schema <- names(imp_bl_ipd)

# Import AgD baseline characteristics
z_imp <- "../additional_covariates/Outputs/baseline_characteristics_trial_wide.csv"
imp_bl_agd <- read_csv(z_imp, show_col_types = FALSE) %>% 
  rename(trial_id = nct_id)

# Covariate sets ---------------------------------------------------------------

# Variables with mean/SD
num_vars <- c(
  "age",
  "bmi",
  "chol_hdl",
  "chol_total",
  "hba1c",
  "durn",
  "egfr",
  "sbp"
)

# Variables represented as trial-level proportions
prop_vars <- c(
  "cursmoke",
  "female_p",
  "white",
  "black_aa",
  "asian",
  "other"
)

# Race variables in extracted data
race_vars <- c(
  "race_ethnicity_amer_ind_alask_nat_proportion",
  "race_ethnicity_asi_proportion",
  "race_ethnicity_black_or_aa_proportion",
  "race_ethnicity_multiracial_proportion",
  "race_ethnicity_nat_hawa_pac_isla_proportion",
  "race_ethnicity_other_proportion",
  "race_ethnicity_whi_proportion"
)

# Race variables to collapse into other
other_race_vars <- c(
  "race_ethnicity_amer_ind_alask_nat_proportion",
  "race_ethnicity_multiracial_proportion",
  "race_ethnicity_nat_hawa_pac_isla_proportion",
  "race_ethnicity_other_proportion"
)

# AgD baseline characteristics -------------------------------------------------

# Get 342-trial AttRisk AgD set and retain target characteristics
agd_set <- imp_bl_agd %>% 
  inner_join(
    imp_agd %>% 
      transmute(
        trial_id,
        trttype = recode(
          as.character(trttype),
          placebo = "plc",
          active = "ac",
          placebo_active = "plc_ac"
        )
      ),
    by = "trial_id"
  ) %>% 
  select(
    trial_id, trttype, n_total,
    age_years_mean, age_years_sd,
    bmi_kg_m_2_mean, bmi_kg_m_2_sd,
    total_cholesterol_mmol_l_mean, total_cholesterol_mmol_l_sd,
    hdl_mg_dl_mean, hdl_mg_dl_sd,
    hdl_mmol_l_mean, hdl_mmol_l_sd,
    duration_t2dm_years_mean, duration_t2dm_years_sd,
    egfr_ml_min_1_73m_2_mean, egfr_ml_min_1_73m_2_sd,
    hba1c_percent_mean, hba1c_percent_sd,
    blood_pressure_systolic_mmhg_mean,
    blood_pressure_systolic_mmhg_sd,
    smoking_status_current_proportion,
    smoking_status_former_proportion,
    smoking_status_former_never_proportion,
    smoking_status_never_proportion,
    all_of(race_vars),
    sex_male_proportion
  )

# Check all AttRisk AgD trials retained
stopifnot(
  nrow(agd_set) == nrow(imp_agd),
  n_distinct(agd_set$trial_id) == nrow(imp_agd),
  setequal(agd_set$trial_id, imp_agd$trial_id)
)

# Harmonise names, units and categories with IPD summary
# Convert HDL from mg/dL to mmol/L where needed
# Set unreported race categories to zero only where race was reported
agd_tidy <- agd_set %>% 
  mutate(
    race_reported = if_any(all_of(race_vars), ~ !is.na(.x)),
    race_other = if_else(
      race_reported,
      rowSums(pick(all_of(other_race_vars)), na.rm = TRUE),
      NA_real_
    )
  ) %>% 
  transmute(
    trial_id,
    trttype,
    n = n_total,
    age_mu = age_years_mean,
    bmi_mu = bmi_kg_m_2_mean,
    chol_hdl_mu = coalesce(hdl_mmol_l_mean, hdl_mg_dl_mean / 38.67),
    chol_total_mu = total_cholesterol_mmol_l_mean,
    durn_mu = duration_t2dm_years_mean,
    egfr_mu = egfr_ml_min_1_73m_2_mean,
    hba1c_mu = hba1c_percent_mean,
    sbp_mu = blood_pressure_systolic_mmhg_mean,
    age_d = age_years_sd,
    bmi_d = bmi_kg_m_2_sd,
    chol_hdl_d = coalesce(hdl_mmol_l_sd, hdl_mg_dl_sd / 38.67),
    chol_total_d = total_cholesterol_mmol_l_sd,
    durn_d = duration_t2dm_years_sd,
    egfr_d = egfr_ml_min_1_73m_2_sd,
    hba1c_d = hba1c_percent_sd,
    sbp_d = blood_pressure_systolic_mmhg_sd,
    cursmoke = case_when(
      !is.na(smoking_status_current_proportion) ~ smoking_status_current_proportion,
      !is.na(smoking_status_former_never_proportion) ~ 1 - smoking_status_former_never_proportion,
      !is.na(smoking_status_former_proportion)
      & !is.na(smoking_status_never_proportion) ~ 1 - smoking_status_former_proportion - smoking_status_never_proportion,
      TRUE ~ NA_real_
    ),
    white = if_else(
      race_reported,
      coalesce(race_ethnicity_whi_proportion, 0),
      NA_real_
    ),
    black_aa = if_else(
      race_reported,
      coalesce(race_ethnicity_black_or_aa_proportion, 0),
      NA_real_
    ),
    asian = if_else(
      race_reported,
      coalesce(race_ethnicity_asi_proportion, 0),
      NA_real_
    ),
    other = race_other,
    female_p = 1 - sex_male_proportion
  )

# Most trials missing at least one baseline summary
check_nas(agd_tidy) %>% arrange(desc(value)) %>% print(n = Inf)

# Combine together -------------------------------------------------------------

# Merge IPD and AgD trial baseline characteristics
shared_vars <- intersect(names(agd_tidy), names(imp_bl_ipd)) %>% 
  setdiff("trial_id")

cmbn_join <- agd_tidy %>% 
  left_join(imp_bl_ipd, by = "trial_id")

# Coalesce duplicate variables
for (v in shared_vars) {
  cmbn_join[[v]] <- coalesce(
    cmbn_join[[paste0(v, ".x")]],
    cmbn_join[[paste0(v, ".y")]]
  )
}

cmbn_join <- cmbn_join %>% 
  select(
    -all_of(c(paste0(shared_vars, ".x"), paste0(shared_vars, ".y")))
  ) %>% 
  select(trial_id, all_of(shared_vars), everything())

# Check missingness
check_nas(cmbn_join) %>% arrange(desc(value)) %>% print(n = Inf)

# MICE setup -------------------------------------------------------------------

# Variables imputed directly using predictive mean matching
num_summary_vars <- c(
  paste0(num_vars, "_mu"),
  paste0(num_vars, "_d")
)

standalone_prop_vars <- c("cursmoke", "female_p")
race_comp_vars <- c("white", "black_aa", "asian", "other")

# Check join retained order
stopifnot(identical(cmbn_join$trial_id, agd_tidy$trial_id))

# Constrain observed proportions between 0 and 1
cmbn_join <- cmbn_join %>% 
  mutate(
    across(
      all_of(standalone_prop_vars),
      ~ pmin(pmax(.x, 0), 1)
    )
  )

# Prepare race composition
# Race is either jointly observed or jointly missing
race_mat <- cmbn_join %>% 
  select(all_of(race_comp_vars)) %>% 
  as.matrix()

race_n_missing <- rowSums(!is.finite(race_mat))
stopifnot(all(race_n_missing %in% c(0L, length(race_comp_vars))))

# Treat all-zero race compositions as missing
race_total_raw <- rowSums(race_mat, na.rm = TRUE)
race_missing <- race_n_missing == length(race_comp_vars) |
  !is.finite(race_total_raw) |
  race_total_raw <= 0

race_observed <- !race_missing
race_total <- race_total_raw[race_observed]

# Close observed race compositions
race_closed <- race_mat
race_closed[race_observed, ] <- sweep(
  race_mat[race_observed, , drop = FALSE],
  MARGIN = 1,
  STATS = race_total,
  FUN = "/"
)

# Transform race using other as additive log-ratio reference
# Small replacement allows zero proportions to be transformed
race_zero <- 1e-6
race_alr_vars <- c("white_alr", "black_aa_alr", "asian_alr")

race_alr <- matrix(
  NA_real_,
  nrow = nrow(cmbn_join),
  ncol = length(race_alr_vars),
  dimnames = list(NULL, race_alr_vars)
)

race_alr[race_observed, ] <- log(
  sweep(
    race_closed[race_observed, c("white", "black_aa", "asian"), drop = FALSE] +
      race_zero,
    MARGIN = 1,
    STATS = race_closed[race_observed, "other"] + race_zero,
    FUN = "/"
  )
)

# Prepare MICE data
# Exclude trial ID but retain all baseline characteristics as predictors
mice_data <- cmbn_join %>% 
  select(
    trttype,
    n,
    all_of(num_summary_vars),
    all_of(standalone_prop_vars)
  ) %>% 
  mutate(
    trttype = factor(trttype),
    across(
      where(is.numeric),
      ~ replace(.x, !is.finite(.x), NA_real_)
    )
  ) %>% 
  bind_cols(as_tibble(race_alr))

# Get previous imputation flag
get_prior_flag <- function(flag_nm) {
  
  if (!flag_nm %in% names(cmbn_join)) {
    return(rep(FALSE, nrow(cmbn_join)))
  }
  
  out <- as.logical(cmbn_join[[flag_nm]])
  out[is.na(out)] <- FALSE
  out
}

# Flag imputed numeric summaries
make_num_flag <- function(var) {
  
  mu_nm <- paste0(var, "_mu")
  d_nm <- paste0(var, "_d")
  flag_nm <- paste0(var, "_imputed")
  
  agd_missing <- !is.finite(agd_tidy[[mu_nm]]) |
    !is.finite(agd_tidy[[d_nm]])
  joined_missing <- !is.finite(cmbn_join[[mu_nm]]) |
    !is.finite(cmbn_join[[d_nm]])
  used_ipd <- agd_missing & !joined_missing
  
  joined_missing | (used_ipd & get_prior_flag(flag_nm))
}

# Flag imputed proportions
make_prop_flag <- function(var) {
  
  flag_nm <- paste0(var, "_imputed")
  agd_missing <- !is.finite(agd_tidy[[var]])
  joined_missing <- !is.finite(cmbn_join[[var]])
  used_ipd <- agd_missing & !joined_missing
  
  joined_missing | (used_ipd & get_prior_flag(flag_nm))
}

# Apply flags
imputation_flags <- tibble(.rows = nrow(cmbn_join))

for (v in num_vars) {
  imputation_flags[[paste0(v, "_imputed")]] <- make_num_flag(v)
}

for (v in prop_vars) {
  imputation_flags[[paste0(v, "_imputed")]] <- make_prop_flag(v)
}

# Flag all-zero race compositions as imputed
for (v in race_comp_vars) {
  flag_nm <- paste0(v, "_imputed")
  imputation_flags[[flag_nm]] <- imputation_flags[[flag_nm]] | race_missing
}

# Set MICE method and predictor matrix
mice_method <- rep("", ncol(mice_data))
names(mice_method) <- names(mice_data)

impute_vars <- names(mice_data)[colSums(is.na(mice_data)) > 0]
mice_method[impute_vars] <- "pmm"

predictor_matrix <- mice::make.predictorMatrix(mice_data)
predictor_matrix[mice_method == "", ] <- 0L
diag(predictor_matrix) <- 0L

# MICE parameters
n_imputations <- 20L
n_iterations <- 20L
imputation_seed <- 20260721L

# Impute
agd_mice <- mice::mice(
  data = mice_data,
  m = n_imputations,
  maxit = n_iterations,
  method = mice_method,
  predictorMatrix = predictor_matrix,
  visitSequence = impute_vars,
  donors = 5L,
  ridge = 1e-5,
  seed = imputation_seed,
  printFlag = TRUE
)

# Print logged events
if (!is.null(agd_mice$loggedEvents) && nrow(agd_mice$loggedEvents) > 0) {
  print(agd_mice$loggedEvents)
}

# Complete imputed datasets ----------------------------------------------------

# Get pooled IPD summaries for comorbidities and concomitant medications
# Used as Poisson rates to simulate participant-level counts in each trial
count_vars <- c("comorbs_mu", "conmeds_mu", "comorbs_d", "conmeds_d")
ipd_schema <- union(ipd_schema, count_vars)

imp_count_sum <- read_csv(
  "vivli/v2/Summaries/overall_summary.csv",
  show_col_types = FALSE
) %>% 
  filter(name %in% c("comorbs", "conmeds")) %>% 
  transmute(name, lambda = parse_number(value))

count_lambda <- set_names(imp_count_sum$lambda, imp_count_sum$name)

# Check count summaries
stopifnot(
  setequal(names(count_lambda), c("comorbs", "conmeds")),
  all(is.finite(count_lambda)),
  all(count_lambda > 0)
)

count_seed <- imputation_seed + 100000L

# Simulate count summary for one trial
simulate_count_summary <- function(n, lambda) {
  
  n <- max(as.integer(round(n)), 2L)
  draws <- rpois(n, lambda = lambda)
  
  c(mu = mean(draws), d = sd(draws))
}

# Simulate comorbidity and concomitant medication summaries for all trials
simulate_arm_count_summaries <- function(n, imp) {
  
  set.seed(count_seed + imp)
  
  tibble(n = n) %>% 
    mutate(
      comorbs = map(
        n,
        ~ simulate_count_summary(.x, count_lambda[["comorbs"]])
      ),
      conmeds = map(
        n,
        ~ simulate_count_summary(.x, count_lambda[["conmeds"]])
      )
    ) %>% 
    transmute(
      comorbs_mu = map_dbl(comorbs, ~ .x[["mu"]]),
      conmeds_mu = map_dbl(conmeds, ~ .x[["mu"]]),
      comorbs_d = map_dbl(comorbs, ~ .x[["d"]]),
      conmeds_d = map_dbl(conmeds, ~ .x[["d"]])
    )
}

# Complete one imputed dataset
complete_agd_imputation <- function(imp) {
  
  completed <- mice::complete(agd_mice, action = imp) %>% 
    as_tibble()
  
  count_complete <- simulate_arm_count_summaries(completed$n, imp)
  
  # Invert additive log-ratio transformation
  exp_alr <- exp(as.matrix(completed[, race_alr_vars]))
  alr_denominator <- 1 + rowSums(exp_alr)
  
  race_imputed <- cbind(
    white = exp_alr[, "white_alr"] / alr_denominator,
    black_aa = exp_alr[, "black_aa_alr"] / alr_denominator,
    asian = exp_alr[, "asian_alr"] / alr_denominator,
    other = 1 / alr_denominator
  )
  
  race_complete <- race_closed
  race_complete[race_missing, ] <- race_imputed[race_missing, ]
  
  # Combine
  tibble(trial_id = cmbn_join$trial_id) %>% 
    bind_cols(
      completed %>% 
        select(
          trttype,
          n,
          all_of(num_summary_vars),
          all_of(standalone_prop_vars)
        )
    ) %>% 
    bind_cols(as_tibble(race_complete)) %>% 
    bind_cols(imputation_flags) %>% 
    bind_cols(count_complete) %>% 
    mutate(
      across(
        all_of(standalone_prop_vars),
        ~ pmin(pmax(.x, 0), 1)
      )
    ) %>% 
    select(all_of(ipd_schema))
}

# Complete all imputed datasets
agd_bl_mi <- map_dfr(
  seq_len(n_imputations),
  ~ complete_agd_imputation(.x) %>% 
    mutate(.imp = .x, .before = 1)
)

# Validate completed datasets --------------------------------------------------

# Check race sums to 1
race_total_mi <- agd_bl_mi %>% 
  transmute(total = white + black_aa + asian + other) %>% 
  pull(total)

# Final checks
stopifnot(
  nrow(agd_bl_mi) == nrow(cmbn_join) * n_imputations,
  n_distinct(agd_bl_mi$.imp) == n_imputations,
  all(count(agd_bl_mi, .imp)$n == nrow(cmbn_join)),
  all(count(agd_bl_mi, .imp, trial_id)$n == 1L),
  !anyNA(agd_bl_mi),
  all(abs(race_total_mi - 1) < 1e-8),
  all(
    as.matrix(
      agd_bl_mi %>% select(all_of(c(standalone_prop_vars, race_comp_vars)))
    ) >= 0
  ),
  all(
    as.matrix(
      agd_bl_mi %>% select(all_of(c(standalone_prop_vars, race_comp_vars)))
    ) <= 1
  ),
  all(
    as.matrix(agd_bl_mi %>% select(ends_with("_d"))) >= 0
  )
)

# Use first imputed dataset as reference population
# Keep all imputations for sensitivity analyses
agd_bl_ref <- agd_bl_mi %>% 
  filter(.imp == 1L) %>% 
  select(-.imp)

# Check reference population
stopifnot(
  nrow(agd_bl_ref) == nrow(cmbn_join),
  !anyDuplicated(agd_bl_ref$trial_id),
  !anyNA(agd_bl_ref)
)

# Save outputs -----------------------------------------------------------------

# Save MICE objects
saveRDS(agd_mice, "processed_data/agd_bl_mice.rds")
saveRDS(agd_bl_mi, "processed_data/agd_bl_imputations.rds")
saveRDS(agd_bl_ref, "processed_data/agd_bl_reference.rds")

# Save manuscript summaries
write_csv(
  agd_bl_mi %>% rename(nct_id = trial_id),
  "outputs/manuscript_summaries/agd_bl_sum_mi.csv"
)

write_csv(
  agd_bl_ref %>% rename(nct_id = trial_id),
  "outputs/manuscript_summaries/agd_bl_sum.csv"
)
