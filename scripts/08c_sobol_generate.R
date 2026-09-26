
# Generate participant profiles for each control arm
# Uses baseline summaries pooled across imputations to capture uncertainty
# Pool profiles across completed datasets for conduct-stage calibration and
# comparison with pooled AgD reference distributions

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Conduct-stage path
conduct_path <- "outputs/conduct_stage/"

# Filepath to nma repo
nma_path2 <- "../nma_more_covariates/"

# Import completed AgD baseline summaries
imp_dist_mi <- readRDS(paste0(conduct_path, "agd_bl_imputations_nested.rds"))

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import pooled IPD-derived rank correlation matrix from previous NMA repo
imp_cor <- readRDS(paste0(nma_path2, "Scratch_data/correlation_summary.rds"))

# Import participant-level log-hazard estimates
# Keep participant-level coefficients fixed across imputations
imp_loghr <- readRDS("outputs/conduct/ma/out/ma_tidy.rds") %>%
  rename(loghr = estimate)

# Import covariate scaling
imp_scale <- read_csv("vivli/Summaries/covariate_mu_d.csv")

# Number of Sobol points per trial and imputation
n_int <- 64L

# Use common randomisation across imputations where applicable
sobol_seed <- 20260721L

# Check baseline summaries ---------------------------------------------------

# Required variables
required_dist_cols <- c(
  "trial_id", "n",
  "age_mu", "age_d",
  "bmi_mu", "bmi_d",
  "chol_hdl_mu", "chol_hdl_d",
  "chol_total_mu", "chol_total_d",
  "durn_mu", "durn_d",
  "egfr_mu", "egfr_d",
  "hba1c_mu", "hba1c_d",
  "sbp_mu", "sbp_d",
  "cursmoke", "female_p",
  "white", "black_aa", "asian", "other",
  "comorbs_mu", "comorbs_d",
  "conmeds_mu", "conmeds_d"
)

# Get AgD trial order
trial_order <- imp_agd %>% transmute(trial_id, .trial_order = row_number())

stopifnot(
  all(c(".imp", "agd_bl") %in% names(imp_dist_mi)),
  !anyDuplicated(imp_dist_mi$.imp),
  !anyDuplicated(trial_order$trial_id),
  all(
    map_lgl(
      imp_dist_mi$agd_bl,
      ~ setequal(.x$trial_id, trial_order$trial_id)
    )
  )
)

# Covariate sets -------------------------------------------------------------

race_vars <- c("white", "black_aa", "asian", "other")
prop_vars <- c("cursmoke", "female_p", race_vars)

positive_mu_vars <- c(
  "bmi_mu", "chol_hdl_mu", "chol_total_mu", "durn_mu",
  "egfr_mu", "sbp_mu", "comorbs_mu", "conmeds_mu"
)

# Names in the IPD export (left) and names used by add_integration (right)
ipd_cor_lookup <- c(
  age10 = "age",
  sex = "male_p",
  bmi = "bmi",
  durn = "durn",
  sbp = "sbp",
  chol_hdl = "chol_hdl",
  chol_total = "chol_total",
  hba1c_base = "hba1c",
  egfr = "egfr",
  cursmoke = "cursmoke"
)

# Race/ethnicity, comorbidities, and conmeds treated as independent
# i.e. no correlation with other characteristics
independent_cor_vars <- c("race_cat", "comorbs", "conmeds")
cor_vars <- c(unname(ipd_cor_lookup), independent_cor_vars)

# Prepare participant-level coefficients -------------------------------------

# Get pooled coefficients
citl_ma <- imp_loghr %>%
  select(term, loghr) %>%
  pivot_wider(names_from = "term", values_from = "loghr") %>%
  rename_with(~ paste0(.x, "_loghr"))

# Expected coefficient names
required_citl_cols <- c(
  "age10_loghr", "bmi_loghr", "chol_hdl_loghr", "chol_total_loghr",
  "smoke_loghr", "durn_loghr", "egfr_loghr", "hba1c_loghr", 
  "raceasian_loghr", "raceblack_aa_loghr", "raceother_na_unk_loghr",
  "sbp_loghr", "sex_loghr", "comorbs_loghr", "conmeds_loghr"
)

missing_citl_cols <- setdiff(required_citl_cols, names(citl_ma))

if (length(missing_citl_cols) > 0) {
  stop(
    "Pooled MA output missing: ",
    paste(missing_citl_cols, collapse = ", ")
  )
}

# Rename pooled betas
citl_beta <- c(
  age10 = citl_ma$age10_loghr,
  bmi = citl_ma$bmi_loghr,
  chol_hdl = citl_ma$chol_hdl_loghr,
  chol_total = citl_ma$chol_total_loghr,
  cursmoke = citl_ma$smoke_loghr,
  durn = citl_ma$durn_loghr,
  egfr = citl_ma$egfr_loghr,
  hba1c = citl_ma$hba1c_loghr,
  asian = citl_ma$raceasian_loghr,
  black_aa = citl_ma$raceblack_aa_loghr,
  other = citl_ma$raceother_na_unk_loghr,
  sbp = citl_ma$sbp_loghr,
  sex = citl_ma$sex_loghr,
  comorbs = citl_ma$comorbs_loghr,
  conmeds = citl_ma$conmeds_loghr
)

stopifnot(!anyNA(citl_beta), all(is.finite(citl_beta)))

# Helper functions -----------------------------------------------------------

# Centre variables on mean and scale by standard deviation
scale_centre_cov <- function(dat, imp_scale) {

  scale_lkp <- imp_scale %>%
    filter(name %in% names(dat)) %>%
    distinct(name, .keep_all = TRUE)

  stopifnot(
    all(is.finite(scale_lkp$mu)),
    all(is.finite(scale_lkp$d)),
    all(scale_lkp$d > 0)
  )

  mu <- setNames(scale_lkp$mu, scale_lkp$name)
  d <- setNames(scale_lkp$d, scale_lkp$name)

  dat %>%
    mutate(
      across(
        all_of(scale_lkp$name),
        ~ (.x - mu[[cur_column()]]) / d[[cur_column()]]
      )
    )

}

# Prepare IPD correlation matrix -----------------------------------------------

# Keep cor for dependent variables 
imp_cor <- as.matrix(imp_cor)
ipd_source_vars <- names(ipd_cor_lookup)
ipd_profile_vars <- unname(ipd_cor_lookup)
missing_ipd_rows <- setdiff(ipd_source_vars, rownames(imp_cor))
missing_ipd_cols <- setdiff(ipd_source_vars, colnames(imp_cor))

if (length(missing_ipd_rows) > 0L || length(missing_ipd_cols) > 0L) {
  stop(
    "IPD correlation matrix missing required variables: ",
    paste(unique(c(missing_ipd_rows, missing_ipd_cols)), collapse = ", ")
  )
}

stopifnot(
  !is.null(rownames(imp_cor)),
  !is.null(colnames(imp_cor)),
  identical(rownames(imp_cor), colnames(imp_cor)),
  !anyDuplicated(rownames(imp_cor)),
  all(is.finite(imp_cor)),
  all(abs(imp_cor) <= 1 + 1e-8),
  all(abs(diag(imp_cor) - 1) < 1e-8),
  isSymmetric(imp_cor, tol = 1e-10),
  all(eigen(imp_cor, symmetric = TRUE, only.values = TRUE)$values > 0)
)

ipd_cor_used <- imp_cor[ipd_source_vars, ipd_source_vars, drop = FALSE]
dimnames(ipd_cor_used) <- list(ipd_profile_vars, ipd_profile_vars)

# Ensure vars are ordered correctly
sobol_cor <- matrix(
  0,
  nrow = length(cor_vars),
  ncol = length(cor_vars),
  dimnames = list(cor_vars, cor_vars)
)
diag(sobol_cor) <- 1
sobol_cor[ipd_profile_vars, ipd_profile_vars] <- ipd_cor_used

# Check scales are correct for Gaussian-copula scale conversion used by
# multinma add_integration() with cor_adjust = "spearman"
sobol_gaussian_cor <- 2 * sin(pi * sobol_cor / 6)
diag(sobol_gaussian_cor) <- 1

stopifnot(
  identical(rownames(sobol_cor), cor_vars),
  identical(colnames(sobol_cor), cor_vars),
  isSymmetric(sobol_cor, tol = 1e-10),
  all(abs(diag(sobol_cor) - 1) < 1e-10),
  all(eigen(sobol_cor, symmetric = TRUE, only.values = TRUE)$values > 0),
  isSymmetric(sobol_gaussian_cor, tol = 1e-10),
  all(eigen(
    sobol_gaussian_cor,
    symmetric = TRUE,
    only.values = TRUE
  )$values > 0)
)

cor_diagnostics <- tibble(
  n_source_variables = nrow(imp_cor),
  n_ipd_variables_used = length(ipd_profile_vars),
  n_independent_variables = length(independent_cor_vars),
  independent_variables = paste(independent_cor_vars, collapse = ", "),
  min_eigen_source_rank = min(
    eigen(imp_cor, symmetric = TRUE, only.values = TRUE)$values
  ),
  min_eigen_used_rank = min(
    eigen(sobol_cor, symmetric = TRUE, only.values = TRUE)$values
  ),
  min_eigen_used_gaussian = min(
    eigen(sobol_gaussian_cor, symmetric = TRUE, only.values = TRUE)$values
  )
)

# Quantile function for mutually exclusive race categories
q_race <- function(p, white, black_aa, asian, other) {
  white <- pmax(as.numeric(white), 0)
  black_aa <- pmax(as.numeric(black_aa), 0)
  asian <- pmax(as.numeric(asian), 0)
  other <- pmax(as.numeric(other), 0)
  total <- white + black_aa + asian + other
  white <- white / total
  black_aa <- black_aa / total
  asian <- asian / total
  cut_white <- white
  cut_black <- white + black_aa
  cut_asian <- white + black_aa + asian
  ifelse(
    p < cut_white,
    0L,
    ifelse(p < cut_black, 1L, ifelse(p < cut_asian, 2L, 3L))
  )
}

# Generate profiles for each completed dataset ---------------------------------

generate_sobol_imputation <- function(.imp, imp_dist) {
  message(sprintf("Generating Sobol profiles for imputation %d", .imp))
  missing_dist_cols <- setdiff(required_dist_cols, names(imp_dist))
  if (length(missing_dist_cols) > 0) {
    stop(
      "Baseline summaries missing: ",
      paste(missing_dist_cols, collapse = ", ")
    )
  }

  # Match AgD trial order
  imp_dist <- trial_order %>%
    left_join(imp_dist, by = "trial_id") %>%
    arrange(.trial_order) %>%
    select(-.trial_order)

  race_total <- rowSums(as.matrix(imp_dist[race_vars]))

  stopifnot(
    identical(imp_dist$trial_id, trial_order$trial_id),
    !anyNA(imp_dist[required_dist_cols]),
    all(
      is.finite(as.matrix(
        imp_dist[setdiff(required_dist_cols, "trial_id")]
      ))
    ),
    all(as.matrix(imp_dist[prop_vars]) >= 0),
    all(as.matrix(imp_dist[prop_vars]) <= 1),
    all(abs(race_total - 1) < 1e-8),
    all(imp_dist$n > 0),
    all(as.matrix(imp_dist[positive_mu_vars]) > 0)
  )

  # Prepare distribution parameters
  prep_dist <- imp_dist %>%
    mutate(
      across(ends_with("_d"), ~ pmax(.x, .Machine$double.eps)),
      bmi_meanlog = get_lnorm_meanlog(bmi_mu, bmi_d),
      bmi_sdlog = get_lnorm_sdlog(bmi_mu, bmi_d),
      chol_hdl_meanlog = get_lnorm_meanlog(chol_hdl_mu, chol_hdl_d),
      chol_hdl_sdlog = get_lnorm_sdlog(chol_hdl_mu, chol_hdl_d),
      sbp_meanlog = get_lnorm_meanlog(sbp_mu, sbp_d),
      sbp_sdlog = get_lnorm_sdlog(sbp_mu, sbp_d),
      chol_total_shape = get_gamma_shape(chol_total_mu, chol_total_d),
      chol_total_rate = get_gamma_rate(chol_total_mu, chol_total_d),
      durn_shape = get_weibull_shape(durn_mu, durn_d),
      durn_scale = get_weibull_scale(durn_mu, durn_shape),
      egfr_shape = get_weibull_shape(egfr_mu, egfr_d),
      egfr_scale = get_weibull_scale(egfr_mu, egfr_shape)
    )

  # Apply the same randomisation seed to each completed dataset
  set.seed(sobol_seed)

  sobol_gen <- add_integration(
    prep_dist,
    age = distr(qnorm, mean = age_mu, sd = age_d),
    male_p = distr(qbern, prob = 1 - female_p),
    bmi = distr(qlnorm, meanlog = bmi_meanlog, sdlog = bmi_sdlog),
    durn = distr(qweibull, shape = durn_shape, scale = durn_scale),
    sbp = distr(qlnorm, meanlog = sbp_meanlog, sdlog = sbp_sdlog),
    chol_hdl = distr(
      qlnorm,
      meanlog = chol_hdl_meanlog,
      sdlog = chol_hdl_sdlog
    ),
    chol_total = distr(
      qgamma,
      shape = chol_total_shape,
      rate = chol_total_rate
    ),
    hba1c = distr(qnorm, mean = hba1c_mu, sd = hba1c_d),
    egfr = distr(qweibull, shape = egfr_shape, scale = egfr_scale),
    cursmoke = distr(qbern, prob = cursmoke),
    race_cat = distr(
      q_race,
      white = white,
      black_aa = black_aa,
      asian = asian,
      other = other
    ),
    comorbs = distr(qpois, lambda = comorbs_mu),
    conmeds = distr(qpois, lambda = conmeds_mu),
    cor = sobol_cor,
    cor_adjust = "spearman",
    n_int = n_int
  )

  int_cols <- names(sobol_gen)[startsWith(names(sobol_gen), ".int_")]
  stopifnot(length(int_cols) == length(cor_vars))

  # Unnest raw profiles
  sobol_raw <- sobol_gen %>%
    select(trial_id, all_of(int_cols)) %>%
    rename_with(~ sub("^\\.int_", "", .x)) %>%
    unnest(cols = -trial_id) %>%
    group_by(trial_id) %>%
    mutate(profile = row_number()) %>%
    ungroup() %>%
    mutate(
      race_cat = as.integer(race_cat),
      white = as.integer(race_cat == 0L),
      black_aa = as.integer(race_cat == 1L),
      asian = as.integer(race_cat == 2L),
      other = as.integer(race_cat == 3L)
    ) %>%
    select(-race_cat) %>%
    relocate(trial_id, profile)

  profiles_per_trial <- sobol_raw %>%
    count(trial_id, name = "n_profiles")

  stopifnot(
    nrow(sobol_raw) == nrow(trial_order) * n_int,
    all(profiles_per_trial$n_profiles == n_int),
    !anyNA(sobol_raw),
    all(rowSums(as.matrix(sobol_raw[race_vars])) == 1)
  )

  # Scale profiles for participant-level linear predictor
  sobol_scaled <- sobol_raw %>%
    rename(
      age10 = age,
      sex = male_p
    ) %>%
    mutate(age10 = age10 / 10) %>%
    scale_centre_cov(imp_scale)

  missing_predictors <- setdiff(names(citl_beta), names(sobol_scaled))
  if (length(missing_predictors) > 0) {
    stop(
      "Generated profiles missing: ",
      paste(missing_predictors, collapse = ", ")
    )
  }

  # Match profiles to AgD trial order
  sobol_scaled_ordered <- sobol_scaled %>%
    left_join(trial_order, by = "trial_id") %>%
    arrange(.trial_order, profile)

  I <- nrow(trial_order)
  J <- n_int
  Pnm <- names(citl_beta)
  pK <- length(Pnm)

  Z <- sobol_scaled_ordered %>%
    select(all_of(Pnm)) %>%
    as.matrix()

  stopifnot(
    nrow(Z) == I * J,
    ncol(Z) == pK,
    all(is.finite(Z)),
    identical(
      unique(sobol_scaled_ordered$trial_id),
      trial_order$trial_id
    )
  )

  # Nest raw and scaled profiles by trial
  raw_nested <- tibble(
    trial_id = trial_order$trial_id
  ) %>%
    mutate(
      sobol_unnest = map(
        trial_id,
        ~ sobol_raw %>%
          filter(trial_id == .x) %>%
          select(-trial_id, -profile)
      )
    )

  scaled_nested <- tibble(
    trial_id = trial_order$trial_id
  ) %>%
    mutate(
      sobol_rdy = map(
        trial_id,
        ~ sobol_scaled %>%
          filter(trial_id == .x) %>%
          select(-trial_id, -profile)
      )
    )

  sobol_unnest <- trial_order %>%
    select(-.trial_order) %>%
    left_join(raw_nested, by = "trial_id") %>%
    left_join(scaled_nested, by = "trial_id")

  # Prepare CITL input
  sobol_lst <- list(
    .imp = .imp,
    I = I,
    J = J,
    Pnm = Pnm,
    pK = pK,
    trial_id = trial_order$trial_id,
    Z = Z,
    loghr = citl_beta
  )

  # Participant-level linear predictors and hazard multipliers
  sobol_lp <- as.vector(Z %*% citl_beta)
  sobol_hr <- exp(sobol_lp)

  sobol_lp_arm <- sobol_scaled_ordered %>%
    transmute(
      arm = .trial_order,
      trial_id,
      profile,
      lp = sobol_lp,
      hr_multiplier = sobol_hr
    ) %>%
    group_by(arm, trial_id) %>%
    summarise(
      mean_lp = mean(lp),
      median_lp = median(lp),
      mean_hr_multiplier = mean(hr_multiplier),
      median_hr_multiplier = median(hr_multiplier),
      .groups = "drop"
    )

  # Profile checks
  profile_checks <- sobol_raw %>%
    summarise(
      n_profiles = n(),
      n_trials = n_distinct(trial_id),
      profiles_per_trial = n_int,
      prop_white = mean(white),
      prop_asian = mean(asian),
      prop_black_aa = mean(black_aa),
      prop_other = mean(other),
      prop_current_smoker = mean(cursmoke),
      prop_female = mean(1 - male_p),
      n_missing_values = sum(is.na(as.matrix(
            select(sobol_raw, where(is.numeric))
      )))
    )

  list(
    dist = imp_dist,
    cor_raw = sobol_cor,
    cor = sobol_cor,
    cor_diagnostics = cor_diagnostics,
    sobol_unnest = sobol_unnest,
    sobol_lst = sobol_lst,
    sobol_lp_arm = sobol_lp_arm,
    profile_checks = profile_checks
  )

}

# Generate profiles across imputations ---------------------------------------

sobol_mi <- imp_dist_mi %>%
  mutate(
    out = map2(.imp, agd_bl, generate_sobol_imputation)
  ) %>%
  transmute(
    .imp,
    dist = map(out, "dist"),
    cor_raw = map(out, "cor_raw"),
    cor = map(out, "cor"),
    cor_diagnostics = map(out, "cor_diagnostics"),
    sobol_unnest = map(out, "sobol_unnest"),
    sobol_lst = map(out, "sobol_lst"),
    sobol_lp_arm = map(out, "sobol_lp_arm"),
    profile_checks = map(out, "profile_checks")
  )

stopifnot(
  nrow(sobol_mi) == nrow(imp_dist_mi),
  identical(sobol_mi$.imp, imp_dist_mi$.imp),
  all(map_int(sobol_mi$sobol_lst, "I") == nrow(trial_order)),
  all(map_int(sobol_mi$sobol_lst, "J") == n_int)
)

# Combine diagnostics --------------------------------------------------------

cor_diagnostics_mi <- map2_dfr(
  sobol_mi$.imp,
  sobol_mi$cor_diagnostics,
  ~ mutate(.y, .imp = .x, .before = 1)
)

profile_checks_mi <- map2_dfr(
  sobol_mi$.imp,
  sobol_mi$profile_checks,
  ~ mutate(.y, .imp = .x, .before = 1)
)

sobol_lp_arm_mi <- map2_dfr(
  sobol_mi$.imp,
  sobol_mi$sobol_lp_arm,
  ~ mutate(.y, .imp = .x, .before = 1)
)

# Save outputs
saveRDS(sobol_mi, paste0(conduct_path, "profiles_imputations_full_data.rds"))

saveRDS(
  sobol_mi %>% select(.imp, sobol_lst),
  paste0(conduct_path, "profiles_imputation.rds")
)

saveRDS(
  sobol_mi %>% select(.imp, sobol_unnest),
  paste0(conduct_path, "profiles_imputation_unnest.rds")
)

write_csv(
  sobol_mi %>%
    select(.imp, sobol_unnest) %>% 
    unnest(sobol_unnest) %>% 
    ungroup() %>% 
    select(-sobol_unnest) %>% 
    unnest(sobol_rdy) %>% 
    ungroup(),
  paste0(conduct_path, "profiles_unnested_per_imputation.csv")
)

# Pool participant profiles across completed datasets ------------------------

# Combine profiles across imputations and control arms
# Each imputation and arm contributes the same number of profiles
bind_profile_sets <- function(profile_column) {
  map2_dfr(
    sobol_mi$.imp,
    sobol_mi$sobol_unnest,
    function(.imp, profile_data) {
      imp_id <- .imp
      profile_data %>%
        transmute(
          trial_id,
          profiles = .data[[profile_column]]
        ) %>%
        unnest(cols = profiles) %>%
        group_by(trial_id) %>%
        mutate(profile = row_number(), .after = trial_id) %>%
        ungroup() %>%
        mutate(.imp = imp_id, .before = 1)
    }
  )
}

# Raw profiles retain the participant-characteristic measurement scales;
# scaled profiles are those supplied to the participant linear predictor
sobol_raw_pooled <- bind_profile_sets("sobol_unnest")
sobol_scaled_pooled <- bind_profile_sets("sobol_rdy")
expected_profile_rows <- nrow(sobol_mi) * nrow(trial_order) * n_int

stopifnot(
  nrow(sobol_raw_pooled) == expected_profile_rows,
  nrow(sobol_scaled_pooled) == expected_profile_rows,
  n_distinct(sobol_raw_pooled$.imp) == nrow(sobol_mi),
  n_distinct(sobol_scaled_pooled$.imp) == nrow(sobol_mi),
  all(count(sobol_raw_pooled, .imp, trial_id)$n == n_int),
  all(count(sobol_scaled_pooled, .imp, trial_id)$n == n_int),
  !anyNA(sobol_raw_pooled),
  !anyNA(sobol_scaled_pooled)
)

# Summarise a pooled participant-profile distribution
summarise_pooled_profiles <- function(data) {
  data %>%
    select(-any_of(c(".imp", "trial_id", "profile"))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    group_by(name) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      q025 = quantile(value, 0.025, na.rm = TRUE, names = FALSE),
      q975 = quantile(value, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

sobol_raw_summary_pooled <- summarise_pooled_profiles(sobol_raw_pooled)
sobol_scaled_summary_pooled <- summarise_pooled_profiles(sobol_scaled_pooled)

write_csv(
  sobol_raw_summary_pooled,
  paste0(conduct_path, "sobol_raw_summary_pooled.csv")
)

write_csv(
  sobol_scaled_summary_pooled,
  paste0(conduct_path, "sobol_scaled_summary_pooled.csv")
)

# Pool AgD reference distributions across completed datasets -----------------

# Give each imputation-arm combination equal weight in the reference mixture
# Matches the fixed number of Sobol profiles generated per arm and imputation
reference_dist_pooled <- map2_dfr(
  imp_dist_mi$.imp,
  imp_dist_mi$agd_bl,
  ~ mutate(.y, .imp = .x, .before = 1)
)

stopifnot(
  nrow(reference_dist_pooled) == nrow(sobol_mi) * nrow(trial_order),
  all(count(reference_dist_pooled, .imp, trial_id)$n == 1L),
  !anyNA(reference_dist_pooled)
)

# Mixture moments for continuous characteristics represented by arm-specific
# means and standard deviations
summarise_continuous_reference <- function(name, component_mean, component_sd) {
  pooled_mean <- mean(component_mean)
  pooled_variance <- mean(component_sd^2 + component_mean^2) - pooled_mean^2
  tibble(
    name = name,
    mean = pooled_mean,
    sd = sqrt(pmax(pooled_variance, 0))
  )
}

# Mixture moments for binary characteristics represented by arm-specific
# probabilities
summarise_binary_reference <- function(name, probability) {
  pooled_probability <- mean(probability)
  tibble(
    name = name,
    mean = pooled_probability,
    sd = sqrt(pooled_probability * (1 - pooled_probability))
  )
}

# Count profiles were generated from arm-specific Poisson distributions
summarise_count_reference <- function(name, lambda) {
  pooled_mean <- mean(lambda)
  pooled_variance <- mean(lambda + lambda^2) - pooled_mean^2
  tibble(
    name = name,
    mean = pooled_mean,
    sd = sqrt(pmax(pooled_variance, 0))
  )
}

continuous_names <- c(
  "age", "bmi", "chol_hdl", "chol_total", "durn", "egfr",
  "hba1c", "sbp"
)

agd_reference_continuous_pooled <- map_dfr(
  continuous_names,
  ~ summarise_continuous_reference(
    name = .x,
    component_mean = reference_dist_pooled[[paste0(.x, "_mu")]],
    component_sd = reference_dist_pooled[[paste0(.x, "_d")]]
  )
)

agd_reference_binary_pooled <- bind_rows(
  summarise_binary_reference(
    "female",
    reference_dist_pooled$female_p
  ),
  summarise_binary_reference(
    "asian",
    reference_dist_pooled$asian
  ),
  summarise_binary_reference(
    "black_aa",
    reference_dist_pooled$black_aa
  ),
  summarise_binary_reference(
    "other",
    reference_dist_pooled$other
  ),
  summarise_binary_reference(
    "white",
    reference_dist_pooled$white
  ),
  summarise_binary_reference(
    "cursmoke",
    reference_dist_pooled$cursmoke
  )
)

agd_reference_count_pooled <- bind_rows(
  summarise_count_reference(
    "comorbs",
    reference_dist_pooled$comorbs_mu
  ),
  summarise_count_reference(
    "conmeds",
    reference_dist_pooled$conmeds_mu
  )
)

table_variable_order <- c(
  "age", "female", "asian",  "black_aa", "other",  "white",  "cursmoke",
  "durn",  "hba1c", "bmi", "egfr", "sbp", "chol_hdl", "chol_total",
  "comorbs", "conmeds"
)

agd_reference_summary_pooled <- bind_rows(
  agd_reference_continuous_pooled,
  agd_reference_binary_pooled,
  agd_reference_count_pooled
) %>%
  arrange(match(name, table_variable_order))

stopifnot(
  identical(agd_reference_summary_pooled$name, table_variable_order),
  !anyNA(agd_reference_summary_pooled),
  all(is.finite(as.matrix(
        agd_reference_summary_pooled %>% select(-name)
  )))
)

write_csv(
  agd_reference_summary_pooled,
  paste0(conduct_path, "agd_reference_summary_pooled.csv")
)

# Compare pooled profiles with pooled AgD reference distributions ------------

# Use female rather than male for presentation in Supplementary Table 11
sobol_table_profile_summary <- sobol_raw_pooled %>%
  mutate(female = 1 - male_p) %>%
  select(all_of(table_variable_order)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "name",
    values_to = "value"
  ) %>%
  group_by(name) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(name, table_variable_order))

variable_lookup <- tribble(
  ~name,
  ~characteristic,
  ~variable_type,
  "age", "Age (years)", "continuous",
  "female", "Female", "binary",
  "asian", "Asian", "binary",
  "black_aa", "Black/African American", "binary",
  "other", "Other race/ethnicity", "binary",
  "white", "White/Caucasian", "binary",
  "cursmoke", "Current smoker", "binary",
  "durn", "Duration of type 2 diabetes (years)", "continuous",
  "hba1c", "HbA1c (%)", "continuous",
  "bmi", "BMI (kg/m^2)", "continuous",
  "egfr", "eGFR (mL/min/1.73 m^2)", "continuous",
  "sbp", "Systolic blood pressure (mmHg)", "continuous",
  "chol_hdl", "HDL cholesterol (mmol/L)", "continuous",
  "chol_total", "Total cholesterol (mmol/L)", "continuous",
  "comorbs", "Number of comorbidities", "count",
  "conmeds", "Number of concomitant medications", "count"
)

profile_reference_comparison_pooled <- agd_reference_summary_pooled %>%
  rename(
    reference_mean = mean,
    reference_sd = sd
  ) %>%
  left_join(
    sobol_table_profile_summary %>%
      rename(
        profile_mean = mean,
        profile_sd = sd
      ),
    by = "name"
  ) %>%
  left_join(variable_lookup, by = "name") %>%
  mutate(
    mean_difference_reference_minus_profiles =
    reference_mean - profile_mean,
    sd_difference_reference_minus_profiles =
    reference_sd - profile_sd,
    difference_percentage_points = if_else(
      variable_type == "binary",
      100 * mean_difference_reference_minus_profiles,
      NA_real_
    )
  ) %>%
  arrange(match(name, table_variable_order))

stopifnot(
  nrow(profile_reference_comparison_pooled) == length(table_variable_order),
  identical(
    profile_reference_comparison_pooled$name,
    table_variable_order
  ),
  !anyNA(
    profile_reference_comparison_pooled %>%
      select(-difference_percentage_points)
  )
)

# Format values for direct insertion into Supplementary Table 11
format_one_decimal <- function(x) {
  rounded <- round(x, 1)
  rounded[abs(rounded) < 0.05] <- 0
  sprintf("%.1f", rounded)
}

compare_diff <- profile_reference_comparison_pooled %>%
  mutate(
    reference_value = case_when(
      variable_type == "binary" ~ paste0(
        format_one_decimal(100 * reference_mean),
        "%"
      ),
      TRUE ~ paste0(
        format_one_decimal(reference_mean),
        " [",
        format_one_decimal(reference_sd),
        "]"
      )
    ),
    profile_value = case_when(
      variable_type == "binary" ~ paste0(
        format_one_decimal(100 * profile_mean),
        "%"
      ),
      TRUE ~ paste0(
        format_one_decimal(profile_mean),
        " [",
        format_one_decimal(profile_sd),
        "]"
      )
    ),
    difference_value = case_when(
      variable_type == "binary" ~ paste0(
        format_one_decimal(difference_percentage_points),
        "%"
      ),
      TRUE ~ paste0(
        format_one_decimal(mean_difference_reference_minus_profiles),
        " [",
        format_one_decimal(sd_difference_reference_minus_profiles),
        "]"
      )
    )
  ) %>%
  transmute(
    `Baseline characteristic` = characteristic,
    `Pooled control-arm reference` = reference_value,
    `Pooled Sobol-generated participant profiles` = profile_value,
    `Difference (reference - profiles)` = difference_value
  )

write_csv(
  compare_diff,
  paste0(conduct_path, "compare_profiles_reference.csv")
)

