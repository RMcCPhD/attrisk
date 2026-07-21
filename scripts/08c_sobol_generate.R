
# Generate Sobol sequences representing synthetic IPD participant sets
# Sequences generated for each AgD trial using trial-specific summaries

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import AgD baseline characteristics
# First imputed dataset used as definitive reference population
imp_dist <- readRDS("processed_data/agd_bl_reference.rds")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import pooled log-hazard estimates from IPD models
imp_loghr <- readRDS("outputs/conduct/ma/out/ma_tidy.rds") %>% 
  rename(loghr = estimate)

# Import covariate scaling
imp_scale <- read_csv(
  "vivli/v2/Summaries/covariate_mu_d.csv",
  show_col_types = FALSE
)

# Output directory
z_out <- "outputs/conduct/sobol"
dir.create(z_out, recursive = TRUE, showWarnings = FALSE)

# Number of Sobol points per trial
n_int <- 64L

# Check baseline summaries -----------------------------------------------------

# Required variables
required_dist_cols <- c(
  "trial_id", "n",
  "age_mu", "age_d", "bmi_mu", "bmi_d",
  "chol_hdl_mu", "chol_hdl_d", "chol_total_mu", "chol_total_d",
  "durn_mu", "durn_d", "egfr_mu", "egfr_d",
  "hba1c_mu", "hba1c_d", "sbp_mu", "sbp_d",
  "cursmoke", "female_p", "white", "black_aa", "asian", "other",
  "comorbs_mu", "comorbs_d", "conmeds_mu", "conmeds_d"
)

# Stop if variables missing
missing_dist_cols <- setdiff(required_dist_cols, names(imp_dist))

if (length(missing_dist_cols) > 0) {
  stop(
    "Baseline reference set missing: ",
    paste(missing_dist_cols, collapse = ", ")
  )
}

# Get AgD trial order
trial_order <- imp_agd %>% 
  transmute(trial_id, .trial_order = row_number())

# Check one summary per AgD trial
stopifnot(
  !anyDuplicated(trial_order$trial_id),
  !anyDuplicated(imp_dist$trial_id),
  setequal(trial_order$trial_id, imp_dist$trial_id)
)

# Match baseline summaries to AgD trial order
imp_dist <- trial_order %>% 
  left_join(imp_dist, by = "trial_id") %>% 
  arrange(.trial_order) %>% 
  select(-.trial_order)

# Define proportions
race_vars <- c("white", "black_aa", "asian", "other")
prop_vars <- c("cursmoke", "female_p", race_vars)
race_total <- rowSums(as.matrix(imp_dist[race_vars]))

# Final checks
stopifnot(
  identical(imp_dist$trial_id, trial_order$trial_id),
  !anyNA(imp_dist[required_dist_cols]),
  all(is.finite(as.matrix(imp_dist[setdiff(required_dist_cols, "trial_id")]))),
  all(as.matrix(imp_dist[prop_vars]) >= 0),
  all(as.matrix(imp_dist[prop_vars]) <= 1),
  all(abs(race_total - 1) < 1e-8),
  all(imp_dist$n > 0)
)

# Helper functions -------------------------------------------------------------

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

# Weighted correlation
weighted_cor <- function(x, y, w) {
  
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]
  
  if (length(x) < 3L) return(NA_real_)
  
  w <- w / sum(w)
  
  mx <- sum(w * x)
  my <- sum(w * y)
  
  vx <- sum(w * (x - mx)^2)
  vy <- sum(w * (y - my)^2)
  
  if (vx <= .Machine$double.eps || vy <= .Machine$double.eps) {
    return(0)
  }
  
  cxy <- sum(w * (x - mx) * (y - my))
  cxy / sqrt(vx * vy)
}

# Weighted Spearman correlation matrix
weighted_cor_mat <- function(data, vars, weight = "n") {
  
  x <- data %>% 
    select(all_of(vars)) %>% 
    mutate(
      across(
        everything(),
        ~ rank(.x, na.last = "keep", ties.method = "average")
      )
    )
  
  out <- matrix(
    NA_real_,
    nrow = length(vars),
    ncol = length(vars),
    dimnames = list(vars, vars)
  )
  
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      out[i, j] <- weighted_cor(
        x[[vars[i]]],
        x[[vars[j]]],
        data[[weight]]
      )
    }
  }
  
  diag(out) <- 1
  out
}

# Ensure correlation matrix is positive definite
make_pd_cor <- function(cor_mat, posd.tol = 1e-6) {
  
  cor_mat <- as.matrix(cor_mat)
  cor_mat[!is.finite(cor_mat)] <- 0
  cor_mat <- (cor_mat + t(cor_mat)) / 2
  diag(cor_mat) <- 1
  
  out <- Matrix::nearPD(
    cor_mat,
    corr = TRUE,
    do2eigen = TRUE,
    posd.tol = posd.tol
  )$mat %>% 
    as.matrix()
  
  out <- (out + t(out)) / 2
  diag(out) <- 1
  dimnames(out) <- dimnames(cor_mat)
  
  out
}

# Distribution parameter helper functions
get_lnorm_meanlog <- function(mu, d) {
  log(mu^2 / sqrt(d^2 + mu^2))
}

get_lnorm_sdlog <- function(mu, d) {
  sqrt(log(1 + d^2 / mu^2))
}

get_gamma_shape <- function(mu, d) {
  (mu / d)^2
}

get_gamma_rate <- function(mu, d) {
  mu / d^2
}

get_weibull_shape <- function(mu, d) {
  
  cv <- d / mu
  
  map_dbl(cv, ~ {
    uniroot(
      function(k) {
        sqrt(gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1) - .x
      },
      interval = c(0.1, 50)
    )$root
  })
}

get_weibull_scale <- function(mu, shape) {
  mu / gamma(1 + 1 / shape)
}

# Quantile function for mutually exclusive race categories
# 0 = White, 1 = Black/AA, 2 = Asian, 3 = Other
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

# Estimate ecological correlation matrix --------------------------------------

# Could not export participant-level correlation matrix due to disclosure risk
# Generate weighted Spearman ecological correlation across all AttRisk trials
# male_p matches downstream model coding
# race_cat is mean categorical race coding in each trial
cor_vars <- c(
  "age",
  "male_p",
  "bmi",
  "durn",
  "sbp",
  "chol_hdl",
  "chol_total",
  "hba1c",
  "egfr",
  "cursmoke",
  "race_cat",
  "comorbs",
  "conmeds"
)

# Transmute summaries
cor_dat <- imp_dist %>% 
  transmute(
    trial_id,
    n,
    age = age_mu,
    male_p = 1 - female_p,
    bmi = bmi_mu,
    durn = durn_mu,
    sbp = sbp_mu,
    chol_hdl = chol_hdl_mu,
    chol_total = chol_total_mu,
    hba1c = hba1c_mu,
    egfr = egfr_mu,
    cursmoke,
    race_cat = black_aa + 2 * asian + 3 * other,
    comorbs = comorbs_mu,
    conmeds = conmeds_mu
  )

# Estimate weighted Spearman correlation
eco_cor_raw <- weighted_cor_mat(
  data = cor_dat,
  vars = cor_vars,
  weight = "n"
)

# Ensure positive definite
eco_cor <- make_pd_cor(eco_cor_raw)

# Get eigenvalues
raw_eigen <- eigen(
  eco_cor_raw,
  symmetric = TRUE,
  only.values = TRUE
)$values

pd_eigen <- eigen(
  eco_cor,
  symmetric = TRUE,
  only.values = TRUE
)$values

# Final checks
stopifnot(
  identical(rownames(eco_cor), cor_vars),
  identical(colnames(eco_cor), cor_vars),
  isSymmetric(eco_cor, tol = 1e-10),
  all(abs(diag(eco_cor) - 1) < 1e-10),
  all(pd_eigen > 0)
)

# Correlation matrix diagnostics
cor_matrix_diagnostics <- tibble(
  n_trials = nrow(cor_dat),
  n_variables = ncol(eco_cor),
  min_eigen_raw_spearman = min(raw_eigen),
  min_eigen_pd_spearman = min(pd_eigen),
  max_abs_pd_adjustment = max(abs(eco_cor - eco_cor_raw)),
  mean_abs_pd_adjustment = mean(abs(eco_cor - eco_cor_raw)),
  n_nonpositive_raw_eigenvalues = sum(raw_eigen <= 0)
)

write_csv(
  cor_matrix_diagnostics,
  file.path(z_out, "sobol_correlation_matrix_diagnostics.csv")
)

# Nest matrices for compatibility with other Sobol scripts
cor_out <- tibble(
  out_raw = list(eco_cor_raw),
  out = list(eco_cor),
  min_eigen = min(pd_eigen)
)

# Save correlation inputs
saveRDS(cor_vars, "created_metadata/sobol_covariates_order.rds")
saveRDS(cor_out, "processed_data/sobol_cor_mat.rds")
saveRDS(imp_dist, "processed_data/sobol_dist.rds")

# Prepare distribution parameters ---------------------------------------------

# Check positive means for non-negative distributions
positive_mu_vars <- c(
  "bmi_mu",
  "chol_hdl_mu",
  "chol_total_mu",
  "durn_mu",
  "egfr_mu",
  "sbp_mu",
  "comorbs_mu",
  "conmeds_mu"
)

stopifnot(all(as.matrix(imp_dist[positive_mu_vars]) > 0))

# Add distribution parameters
prep_dist <- imp_dist %>% 
  mutate(
    across(ends_with("_d"), ~ pmax(.x, .Machine$double.eps)),
    
    # Log-normal: BMI, HDL cholesterol and SBP
    bmi_meanlog = get_lnorm_meanlog(bmi_mu, bmi_d),
    bmi_sdlog   = get_lnorm_sdlog(bmi_mu, bmi_d),
    
    chol_hdl_meanlog = get_lnorm_meanlog(chol_hdl_mu, chol_hdl_d),
    chol_hdl_sdlog   = get_lnorm_sdlog(chol_hdl_mu, chol_hdl_d),
    
    sbp_meanlog = get_lnorm_meanlog(sbp_mu, sbp_d),
    sbp_sdlog   = get_lnorm_sdlog(sbp_mu, sbp_d),
    
    # Gamma: total cholesterol
    chol_total_shape = get_gamma_shape(chol_total_mu, chol_total_d),
    chol_total_rate  = get_gamma_rate(chol_total_mu, chol_total_d),
    
    # Weibull: diabetes duration and eGFR
    durn_shape = get_weibull_shape(durn_mu, durn_d),
    durn_scale = get_weibull_scale(durn_mu, durn_shape),
    
    egfr_shape = get_weibull_shape(egfr_mu, egfr_d),
    egfr_scale = get_weibull_scale(egfr_mu, egfr_shape)
  )

# Generate Sobol sequences -----------------------------------------------------

# Generate 64 trial-specific Sobol profiles per AgD trial
# Correlation arguments follow order in cor_vars
# Spearman adjustment applied within add_integration()
sobol_gen <- add_integration(
  prep_dist,
  age        = distr(qnorm,  mean = age_mu, sd = age_d),
  male_p     = distr(qbern,  prob = 1 - female_p),
  bmi        = distr(qlnorm, meanlog = bmi_meanlog, sdlog = bmi_sdlog),
  durn       = distr(qweibull, shape = durn_shape, scale = durn_scale),
  sbp        = distr(qlnorm, meanlog = sbp_meanlog, sdlog = sbp_sdlog),
  chol_hdl   = distr(qlnorm, meanlog = chol_hdl_meanlog, sdlog = chol_hdl_sdlog),
  chol_total = distr(qgamma, shape = chol_total_shape, rate = chol_total_rate),
  hba1c      = distr(qnorm, mean = hba1c_mu, sd = hba1c_d),
  egfr       = distr(qweibull, shape = egfr_shape, scale = egfr_scale),
  cursmoke   = distr(qbern, prob = cursmoke),
  race_cat   = distr(
    q_race,
    white = white,
    black_aa = black_aa,
    asian = asian,
    other = other
  ),
  comorbs    = distr(qpois, lambda = comorbs_mu),
  conmeds    = distr(qpois, lambda = conmeds_mu),
  cor = eco_cor,
  cor_adjust = "spearman",
  n_int = n_int
)

# Get integration columns
int_cols <- names(sobol_gen)[startsWith(names(sobol_gen), ".int_")]
stopifnot(length(int_cols) == length(cor_vars))

# Unnest Sobol profiles
sobol_raw <- sobol_gen %>% 
  select(trial_id, all_of(int_cols)) %>% 
  rename_with(., ~ gsub(".int_", "", .)) %>% 
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

# Check number of profiles per trial
profiles_per_trial <- sobol_raw %>% 
  count(trial_id, name = "n_profiles")

# Final checks
stopifnot(
  nrow(sobol_raw) == nrow(trial_order) * n_int,
  nrow(profiles_per_trial) == nrow(trial_order),
  all(profiles_per_trial$n_profiles == n_int),
  !anyNA(sobol_raw),
  all(rowSums(as.matrix(sobol_raw[race_vars])) == 1)
)

# Transform Sobol points
# Centre on pooled IPD mean and scale by standard deviation
sobol_scaled <- sobol_raw %>% 
  rename(
    age10 = age,
    sex = male_p
  ) %>% 
  mutate(age10 = age10 / 10) %>% 
  scale_centre_cov(imp_scale)

# Nest raw Sobol profiles
raw_nested <- sobol_raw %>% 
  select(-profile) %>% 
  group_by(trial_id) %>% 
  nest(.key = "sobol_unnest") %>% 
  ungroup()

# Nest scaled Sobol profiles
scaled_nested <- sobol_scaled %>% 
  select(-profile) %>% 
  group_by(trial_id) %>% 
  nest(.key = "sobol_rdy") %>% 
  ungroup()

# Combine
sobol_unnest <- trial_order %>% 
  select(-.trial_order) %>% 
  left_join(raw_nested, by = "trial_id") %>% 
  left_join(scaled_nested, by = "trial_id")

# Save outputs
saveRDS(sobol_unnest, "processed_data/sobol_out.rds")
saveRDS(sobol_unnest, "processed_data/sobol_unnest.rds")

# Inspect distributions
sobol_plot <- sobol_scaled %>% 
  pivot_longer(-c(trial_id, profile)) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  theme_classic()

ggsave(
  file.path(z_out, "sobol_distributions.png"),
  plot = sobol_plot,
  width = 8,
  height = 6,
  units = "in"
)

# Prepare CITL inputs ----------------------------------------------------------

# Get pooled coefficients
citl_ma <- imp_loghr %>% 
  select(term, loghr) %>% 
  pivot_wider(names_from = "term", values_from = "loghr") %>% 
  rename_with(., ~ paste0(., "_loghr"))

# Expected coefficient names
required_citl_cols <- c(
  "age10_loghr",
  "bmi_loghr",
  "chol_hdl_loghr",
  "chol_total_loghr",
  "smoke_loghr",
  "durn_loghr",
  "egfr_loghr",
  "hba1c_loghr",
  "raceasian_loghr",
  "raceblack_aa_loghr",
  "raceother_na_unk_loghr",
  "sbp_loghr",
  "sex_loghr",
  "comorbs_loghr",
  "conmeds_loghr"
)

# Stop if coefficients missing
missing_citl_cols <- setdiff(required_citl_cols, names(citl_ma))

if (length(missing_citl_cols) > 0) {
  stop(
    "Pooled MA output missing: ",
    paste(missing_citl_cols, collapse = ", ")
  )
}

# Rename pooled betas
citl_beta <- c(
  age10       = citl_ma$age10_loghr,
  bmi         = citl_ma$bmi_loghr,
  chol_hdl    = citl_ma$chol_hdl_loghr,
  chol_total  = citl_ma$chol_total_loghr,
  cursmoke    = citl_ma$smoke_loghr,
  durn        = citl_ma$durn_loghr,
  egfr        = citl_ma$egfr_loghr,
  hba1c       = citl_ma$hba1c_loghr,
  asian       = citl_ma$raceasian_loghr,
  black_aa    = citl_ma$raceblack_aa_loghr,
  other       = citl_ma$raceother_na_unk_loghr,
  sbp         = citl_ma$sbp_loghr,
  sex         = citl_ma$sex_loghr,
  comorbs     = citl_ma$comorbs_loghr,
  conmeds     = citl_ma$conmeds_loghr
)

# Check predictors present
missing_predictors <- setdiff(names(citl_beta), names(sobol_scaled))

if (length(missing_predictors) > 0) {
  stop(
    "Generated profiles missing: ",
    paste(missing_predictors, collapse = ", ")
  )
}

stopifnot(!anyNA(citl_beta), all(is.finite(citl_beta)))

# Match scaled profiles to AgD trial order
sobol_scaled_ordered <- sobol_scaled %>% 
  left_join(trial_order, by = "trial_id") %>% 
  arrange(.trial_order, profile)

# Parameters
I <- nrow(trial_order) # Number of trials
J <- n_int             # Number of Sobol points
Pnm <- names(citl_beta) # Predictor names
pK <- length(Pnm)       # Number of predictors

# Construct matrix of Sobol sequences per trial: (I x J) x pK
Z <- sobol_scaled_ordered %>% 
  select(all_of(Pnm)) %>% 
  as.matrix()

# Final checks
stopifnot(
  nrow(Z) == I * J,
  ncol(Z) == pK,
  all(is.finite(Z)),
  identical(
    unique(sobol_scaled_ordered$trial_id),
    trial_order$trial_id
  )
)

# List of outputs
sobol_lst <- list(
  I = I,
  J = J,
  Pnm = Pnm,
  pK = pK,
  Z = Z,
  loghr = citl_beta
)

str(sobol_lst)

# Save
saveRDS(sobol_lst, "processed_data/sobol_lst.rds")

# Supplementary Sobol summaries -----------------------------------------------

# Summarise vector
summ_vec <- function(x, metric) {
  
  tibble(
    metric = metric,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

# Summarise overall distributions
summ_long <- function(dat) {
  
  dat %>% 
    select(-profile) %>% 
    pivot_longer(-trial_id) %>% 
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

# Summarise distributions by trial
summ_long_by_trial <- function(dat) {
  
  dat %>% 
    select(-profile) %>% 
    pivot_longer(-trial_id) %>% 
    group_by(trial_id, name) %>% 
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      q025 = quantile(value, 0.025, na.rm = TRUE, names = FALSE),
      q975 = quantile(value, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

# Summarise generated distributions
sobol_raw_summary <- summ_long(sobol_raw)
sobol_scaled_summary <- summ_long(sobol_scaled)
sobol_raw_by_trial <- summ_long_by_trial(sobol_raw)

# Numeric target summaries
continuous_names <- c(
  "age",
  "bmi",
  "chol_hdl",
  "chol_total",
  "hba1c",
  "durn",
  "egfr",
  "sbp"
)

continuous_targets <- map_dfr(continuous_names, ~ {
  tibble(
    trial_id = imp_dist$trial_id,
    name = .x,
    target_mean = imp_dist[[paste0(.x, "_mu")]],
    target_sd = imp_dist[[paste0(.x, "_d")]],
    target_proportion = NA_real_
  )
})

# Proportion target summaries
proportion_targets <- bind_rows(
  tibble(
    trial_id = imp_dist$trial_id,
    name = "male_p",
    target_proportion = 1 - imp_dist$female_p
  ),
  map_dfr(c("cursmoke", race_vars), ~ {
    tibble(
      trial_id = imp_dist$trial_id,
      name = .x,
      target_proportion = imp_dist[[.x]]
    )
  })
) %>% 
  mutate(
    target_mean = NA_real_,
    target_sd = NA_real_
  ) %>% 
  select(trial_id, name, target_mean, target_sd, target_proportion)

# Count target summaries
count_targets <- map_dfr(c("comorbs", "conmeds"), ~ {
  tibble(
    trial_id = imp_dist$trial_id,
    name = .x,
    target_mean = imp_dist[[paste0(.x, "_mu")]],
    target_sd = imp_dist[[paste0(.x, "_d")]],
    target_proportion = NA_real_
  )
})

# Combine target summaries
sobol_target_summary <- bind_rows(
  continuous_targets,
  proportion_targets,
  count_targets
)

# Compare generated and target distributions
sobol_generation_check <- sobol_raw_by_trial %>% 
  left_join(sobol_target_summary, by = c("trial_id", "name")) %>% 
  mutate(
    mean_minus_target = case_when(
      !is.na(target_mean) ~ mean - target_mean,
      !is.na(target_proportion) ~ mean - target_proportion,
      TRUE ~ NA_real_
    ),
    sd_minus_target = if_else(
      !is.na(target_sd),
      sd - target_sd,
      NA_real_
    )
  )

# Check generated profiles
sobol_profile_check_dat <- sobol_raw %>% 
  mutate(
    n_race_indicators = rowSums(pick(all_of(race_vars))),
    multiple_race_indicators = n_race_indicators > 1,
    male = male_p,
    female = 1 - male_p
  )

sobol_profile_checks <- sobol_profile_check_dat %>% 
  summarise(
    n_profiles = n(),
    n_trials = n_distinct(trial_id),
    profiles_per_trial = n_int,
    prop_multiple_race_indicators = mean(multiple_race_indicators),
    prop_white = mean(white),
    prop_asian = mean(asian),
    prop_black_aa = mean(black_aa),
    prop_other = mean(other),
    prop_current_smoker = mean(cursmoke),
    prop_female_generated = mean(female),
    prop_male_generated = mean(male),
    n_missing_values = sum(is.na(as.matrix(
      select(sobol_profile_check_dat, where(is.numeric))
    )))
  )

# Save summaries
write_csv(sobol_raw_summary, file.path(z_out, "sobol_raw_summary.csv"))
write_csv(sobol_scaled_summary, file.path(z_out, "sobol_scaled_summary.csv"))
write_csv(
  sobol_generation_check,
  file.path(z_out, "sobol_generation_check.csv")
)
write_csv(
  sobol_profile_checks,
  file.path(z_out, "sobol_profile_checks.csv")
)

# Summarise linear predictors --------------------------------------------------

# Participant-level linear predictors and hazard multipliers
sobol_lp <- as.vector(Z %*% citl_beta)
sobol_hr_multiplier <- exp(sobol_lp)

sobol_lp_dat <- sobol_scaled_ordered %>% 
  transmute(
    arm = .trial_order,
    trial_id,
    profile,
    lp = sobol_lp,
    hr_multiplier = sobol_hr_multiplier
  )

# Summarise by trial
sobol_lp_arm <- sobol_lp_dat %>% 
  group_by(trial_id) %>% 
  summarise(
    mean_lp = mean(lp),
    median_lp = median(lp),
    mean_hr_multiplier = mean(hr_multiplier),
    median_hr_multiplier = median(hr_multiplier),
    .groups = "drop"
  )

# Summarise overall
sobol_lp_summary <- bind_rows(
  summ_vec(sobol_lp_dat$lp, "participant_linear_predictor"),
  summ_vec(sobol_lp_dat$hr_multiplier, "participant_hazard_multiplier"),
  summ_vec(sobol_lp_arm$mean_lp, "arm_mean_linear_predictor"),
  summ_vec(sobol_lp_arm$mean_hr_multiplier, "arm_mean_hazard_multiplier")
)

# Save summaries
write_csv(
  sobol_lp_dat,
  file.path(z_out, "sobol_linear_predictors.csv")
)

write_csv(
  sobol_lp_arm,
  file.path(z_out, "sobol_linear_predictors_by_arm.csv")
)

write_csv(
  sobol_lp_summary,
  file.path(z_out, "sobol_linear_predictor_summary.csv")
)
