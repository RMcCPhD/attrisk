
# Check stability of Sobol integration grid size
# Generate alternative trial-specific grids and compare conduct-stage summaries

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import Sobol distribution summaries
imp_dist <- readRDS("processed_data/sobol_dist.rds")

# Import ecological correlation matrix
imp_cor <- readRDS("processed_data/sobol_cor_mat.rds")

# Import primary Sobol outputs
imp_sobol <- readRDS("processed_data/sobol_lst.rds")
imp_sobol_profiles <- readRDS("processed_data/sobol_unnest.rds")

# Import covariate scaling
imp_scale <- read_csv(
  "vivli/v2/Summaries/covariate_mu_d.csv",
  show_col_types = FALSE
)

# Import primary CITL outputs
imp_citl <- readRDS("outputs/conduct/citl/citl_draws.rds")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Output directory
z_out <- "outputs/sensitivity"
dir.create(z_out, recursive = TRUE, showWarnings = FALSE)

# Sobol grid sizes
n_sobol_grid <- c(64L, 128L, 256L, 512L)

# Check inputs -----------------------------------------------------------------

# Required CITL outputs
required_citl_elements <- c(
  "I",
  "J",
  "S",
  "trial_id",
  "sobol_lp",
  "n_citl_grid",
  "target_risk",
  "avg_risk_pre",
  "delta"
)

missing_citl_elements <- setdiff(
  required_citl_elements,
  names(imp_citl)
)

if (length(missing_citl_elements) > 0) {
  stop(
    "CITL outputs missing: ",
    paste(missing_citl_elements, collapse = ", "),
    ". Rerun scripts/09_conduct_citl.R"
  )
}

# Number of trials
I <- nrow(imp_dist)
n_citl_grid <- imp_citl$n_citl_grid

# Check trial order and dimensions
stopifnot(
  I == nrow(imp_agd),
  I == imp_sobol$I,
  I == imp_citl$I,
  identical(imp_dist$trial_id, imp_agd$trial_id),
  identical(imp_dist$trial_id, imp_sobol_profiles$trial_id),
  identical(imp_dist$trial_id, imp_citl$trial_id),
  ncol(imp_citl$target_risk) == I,
  ncol(imp_citl$avg_risk_pre) == I,
  ncol(imp_citl$delta) == I,
  length(imp_citl$sobol_lp) == I * imp_sobol$J,
  length(n_citl_grid) == 1L,
  is.finite(n_citl_grid),
  n_citl_grid > 0,
  imp_sobol$J == 64L,
  imp_citl$J == 64L
)

# Get correlation matrix
cor_mat <- imp_cor$out[[1]]

stopifnot(
  is.matrix(cor_mat),
  identical(rownames(cor_mat), colnames(cor_mat)),
  isSymmetric(cor_mat),
  all(eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values > 0)
)

# Helper functions -------------------------------------------------------------

# Quantiles
q025 <- function(x) {
  quantile(x, 0.025, na.rm = TRUE, names = FALSE)
}

q975 <- function(x) {
  quantile(x, 0.975, na.rm = TRUE, names = FALSE)
}

# Summarise vector
summ_vec <- function(x) {
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = q025(x),
    q975 = q975(x)
  )
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

# Add distribution parameters
make_distribution_parameters <- function(dist) {
  
  dist %>% 
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
}

# Centre variables on mean and scale by standard deviation
scale_centre_cov <- function(dat, imp_scale) {
  
  scale_lkp <- imp_scale %>% 
    filter(name %in% names(dat))
  
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

# Generate trial-specific Sobol profiles
generate_sobol_profiles <- function(n_sobol, dist, cor_mat, imp_scale, loghr) {
  
  # Generate integration points using same setup as primary analysis
  sobol_gen <- add_integration(
    dist,
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
    cor = cor_mat,
    cor_adjust = "spearman",
    n_int = n_sobol
  )
  
  # Unnest integration points
  int_cols <- names(sobol_gen)[startsWith(names(sobol_gen), ".int_")]
  
  raw <- sobol_gen %>% 
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
  
  # Scale profiles
  scaled <- raw %>% 
    rename(
      age10 = age,
      sex = male_p
    ) %>% 
    mutate(age10 = age10 / 10) %>% 
    scale_centre_cov(imp_scale)
  
  # Construct linear predictors
  Z <- scaled %>% 
    select(all_of(names(loghr))) %>% 
    as.matrix()
  
  lp <- as.vector(Z %*% loghr)
  lp_mat <- matrix(lp, nrow = n_sobol, ncol = nrow(dist))
  
  # Final checks
  stopifnot(
    nrow(raw) == nrow(dist) * n_sobol,
    all(count(raw, trial_id)$n == n_sobol),
    identical(unique(raw$trial_id), dist$trial_id),
    all(rowSums(as.matrix(raw[c("white", "black_aa", "asian", "other")])) == 1),
    all(is.finite(Z)),
    all(is.finite(lp_mat))
  )
  
  list(
    raw = raw,
    scaled = scaled,
    Z = Z,
    lp = lp,
    lp_mat = lp_mat,
    hr = exp(lp),
    hr_mat = exp(lp_mat)
  )
}

# Summarise generated covariates
summarise_covariates <- function(dat, n_sobol, scale) {
  
  dat %>% 
    select(-any_of(c("trial_id", "profile"))) %>% 
    pivot_longer(everything(), names_to = "covariate") %>% 
    group_by(covariate) %>% 
    summarise(
      n_sobol = n_sobol,
      scale = scale,
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      q025 = q025(value),
      q975 = q975(value),
      .groups = "drop"
    ) %>% 
    select(n_sobol, scale, covariate, everything())
}

# Participant-averaged survival
laplace_survival <- function(h, hr) {
  rowMeans(exp(-outer(h, hr)))
}

# Create interpolation functions for pre-CITL risk and CITL shift
make_citl_functions <- function(H, hr, n_grid = n_citl_grid, buffer = 6) {
  
  H <- H[is.finite(H) & H > 0]
  
  log_h_grid <- seq(
    log(min(H)) - buffer,
    log(max(H)) + buffer,
    length.out = n_grid
  )
  
  h_grid <- exp(log_h_grid)
  survival_grid <- laplace_survival(h_grid, hr)
  risk_grid <- 1 - survival_grid
  marginal_h_grid <- -log(pmax(survival_grid, .Machine$double.xmin))
  
  # Remove duplicated values following numerical underflow
  keep <- is.finite(marginal_h_grid) & !duplicated(marginal_h_grid)
  
  pre_interp <- approxfun(
    log_h_grid,
    risk_grid,
    rule = 2
  )
  
  delta_interp <- approxfun(
    marginal_h_grid[keep],
    log_h_grid[keep],
    rule = 2
  )
  
  list(
    pre = function(H) pre_interp(log(H)),
    delta = function(H) delta_interp(H) - log(H)
  )
}

# Recalculate CITL summaries for each trial-specific population
calculate_citl <- function(H_mat, hr_mat) {
  
  S <- nrow(H_mat)
  I <- ncol(H_mat)
  
  pre <- matrix(NA_real_, nrow = S, ncol = I)
  delta <- matrix(NA_real_, nrow = S, ncol = I)
  
  for (i in seq_len(I)) {
    
    if (i %% 25 == 1 || i == I) {
      message(sprintf("Trial %d/%d", i, I))
    }
    
    citl_fun <- make_citl_functions(
      H = H_mat[, i],
      hr = hr_mat[, i]
    )
    
    pre[, i] <- citl_fun$pre(H_mat[, i])
    delta[, i] <- citl_fun$delta(H_mat[, i])
  }
  
  list(
    pre = pre,
    delta = delta
  )
}

# Prepare Sobol inputs ---------------------------------------------------------

# Add distribution parameters
dist_par <- make_distribution_parameters(imp_dist)

# Pooled participant-level effects
citl_beta <- imp_sobol$loghr

# Check coefficient names
stopifnot(
  identical(colnames(imp_sobol$Z), names(citl_beta)),
  all(is.finite(citl_beta))
)

# Unnest primary raw profiles
primary_raw <- imp_sobol_profiles %>% 
  select(trial_id, sobol_unnest) %>% 
  unnest(sobol_unnest) %>% 
  group_by(trial_id) %>% 
  mutate(profile = row_number()) %>% 
  ungroup() %>% 
  relocate(trial_id, profile)

# Unnest primary scaled profiles
primary_scaled <- imp_sobol_profiles %>% 
  select(trial_id, sobol_rdy) %>% 
  unnest(sobol_rdy) %>% 
  group_by(trial_id) %>% 
  mutate(profile = row_number()) %>% 
  ungroup() %>% 
  relocate(trial_id, profile)

# Get primary linear predictors
primary_lp <- as.vector(imp_sobol$Z %*% citl_beta)
primary_lp_mat <- matrix(primary_lp, nrow = imp_sobol$J, ncol = I)
primary_hr <- exp(primary_lp)
primary_hr_mat <- exp(primary_lp_mat)

# Check CITL used current primary Sobol profiles
stopifnot(
  max(abs(primary_lp - imp_citl$sobol_lp)) < 1e-12
)

# Design-stage cumulative hazard
target_risk <- imp_citl$target_risk

stopifnot(
  all(target_risk > 0),
  all(target_risk < 1)
)

H_mat <- -log1p(-target_risk)

# Final checks
stopifnot(
  nrow(primary_raw) == I * imp_sobol$J,
  nrow(primary_scaled) == I * imp_sobol$J,
  identical(unique(primary_raw$trial_id), imp_dist$trial_id),
  identical(unique(primary_scaled$trial_id), imp_dist$trial_id),
  all(is.finite(primary_lp_mat)),
  all(is.finite(H_mat))
)

# Check saved CITL outputs match current primary Sobol profiles
check_draws <- unique(
  round(seq(1, nrow(H_mat), length.out = 8))
)

primary_match_error <- map_dbl(seq_len(I), ~ {
  
  pre_check <- 1 - laplace_survival(
    H_mat[check_draws, .x],
    primary_hr_mat[, .x]
  )
  
  max(
    abs(pre_check - imp_citl$avg_risk_pre[check_draws, .x]),
    na.rm = TRUE
  )
}) %>% 
  max()

if (primary_match_error > 1e-8) {
  stop(
    "Saved CITL outputs do not match current Sobol profiles. ",
    "Rerun scripts/09_conduct_citl.R before this script."
  )
}

# Generate alternative grids --------------------------------------------------

grid_objects <- list()
covariate_summaries <- list()
stability_work <- list()

for (n_sobol in n_sobol_grid) {
  
  message(sprintf("Processing n_sobol = %d", n_sobol))
  
  if (n_sobol == imp_sobol$J) {
    
    # Use saved primary Sobol and CITL outputs
    raw <- primary_raw
    scaled <- primary_scaled
    lp <- primary_lp
    hr <- primary_hr
    pre <- imp_citl$avg_risk_pre
    delta <- imp_citl$delta
    source <- "primary_saved_grid"
    
  } else {
    
    # Generate trial-specific alternative grid
    grid <- generate_sobol_profiles(
      n_sobol = n_sobol,
      dist = dist_par,
      cor_mat = cor_mat,
      imp_scale = imp_scale,
      loghr = citl_beta
    )
    
    # Recalculate CITL for each trial
    citl_grid <- calculate_citl(
      H_mat = H_mat,
      hr_mat = grid$hr_mat
    )
    
    raw <- grid$raw
    scaled <- grid$scaled
    lp <- grid$lp
    hr <- grid$hr
    pre <- citl_grid$pre
    delta <- citl_grid$delta
    source <- "generated_grid"
  }
  
  # Store grids
  grid_objects[[as.character(n_sobol)]] <- list(
    n_sobol = n_sobol,
    source = source,
    raw = raw,
    scaled = scaled,
    lp = lp,
    hr = hr
  )
  
  # Summarise covariates
  covariate_summaries[[as.character(n_sobol)]] <- bind_rows(
    summarise_covariates(raw, n_sobol, "raw"),
    summarise_covariates(scaled, n_sobol, "scaled")
  )
  
  # Store CITL outputs
  stability_work[[as.character(n_sobol)]] <- list(
    n_sobol = n_sobol,
    lp = lp,
    hr = hr,
    pre = pre,
    delta = delta
  )
}

# Summarise stability ----------------------------------------------------------

# Use largest grid as reference
ref <- stability_work[[as.character(max(n_sobol_grid))]]

# Overall summaries
stability_table <- map_dfr(stability_work, ~ {
  
  lp_sum <- summ_vec(.x$lp)
  hr_sum <- summ_vec(.x$hr)
  
  pre_diff <- abs(.x$pre - ref$pre)
  delta_diff <- abs(.x$delta - ref$delta)
  
  tibble(
    n_sobol = .x$n_sobol,
    mean_participant_linear_predictor = lp_sum$mean,
    sd_participant_linear_predictor = lp_sum$sd,
    q025_participant_linear_predictor = lp_sum$q025,
    q975_participant_linear_predictor = lp_sum$q975,
    mean_hazard_multiplier = hr_sum$mean,
    sd_hazard_multiplier = hr_sum$sd,
    q025_hazard_multiplier = hr_sum$q025,
    q975_hazard_multiplier = hr_sum$q975,
    mean_pre_citl_participant_marginalised_risk = mean(.x$pre, na.rm = TRUE),
    mean_post_citl_participant_marginalised_risk = mean(target_risk, na.rm = TRUE),
    mean_citl_intercept_shift = mean(.x$delta, na.rm = TRUE),
    mean_abs_difference_pre_citl_risk_vs_n512 = mean(pre_diff, na.rm = TRUE),
    max_abs_difference_pre_citl_risk_vs_n512 = max(pre_diff, na.rm = TRUE),
    mean_abs_difference_citl_intercept_shift_vs_n512 = mean(delta_diff, na.rm = TRUE),
    max_abs_difference_citl_intercept_shift_vs_n512 = max(delta_diff, na.rm = TRUE)
  )
}) %>% 
  arrange(n_sobol)

# Trial-level summaries
stability_by_trial <- map_dfr(stability_work, ~ {
  
  pre_diff <- abs(.x$pre - ref$pre)
  delta_diff <- abs(.x$delta - ref$delta)
  
  tibble(
    n_sobol = .x$n_sobol,
    trial_id = imp_dist$trial_id,
    mean_abs_difference_pre_citl_risk_vs_n512 = colMeans(pre_diff, na.rm = TRUE),
    max_abs_difference_pre_citl_risk_vs_n512 = apply(pre_diff, 2, max, na.rm = TRUE),
    mean_abs_difference_citl_intercept_shift_vs_n512 = colMeans(delta_diff, na.rm = TRUE),
    max_abs_difference_citl_intercept_shift_vs_n512 = apply(
      delta_diff,
      2,
      max,
      na.rm = TRUE
    )
  )
}) %>% 
  arrange(n_sobol, trial_id)

# Covariate summaries
covariate_summary <- bind_rows(covariate_summaries) %>% 
  arrange(n_sobol, scale, covariate)

# Validate interpolation -------------------------------------------------------

# Recalculate primary 64-grid CITL outputs using interpolation
message("Validating primary n_sobol = 64 grid")

primary_check <- calculate_citl(
  H_mat = H_mat,
  hr_mat = primary_hr_mat
)

# Compare to saved primary outputs
primary_validation <- tibble(
  check = c(
    "max_abs_pre_citl_risk_interpolation_error",
    "mean_abs_pre_citl_risk_interpolation_error",
    "max_abs_delta_interpolation_error",
    "mean_abs_delta_interpolation_error"
  ),
  value = c(
    max(abs(primary_check$pre - imp_citl$avg_risk_pre), na.rm = TRUE),
    mean(abs(primary_check$pre - imp_citl$avg_risk_pre), na.rm = TRUE),
    max(abs(primary_check$delta - imp_citl$delta), na.rm = TRUE),
    mean(abs(primary_check$delta - imp_citl$delta), na.rm = TRUE)
  )
)

# Save outputs -----------------------------------------------------------------

# Save summaries
write_csv(
  stability_table,
  file.path(z_out, "sobol_grid_stability.csv")
)

write_csv(
  stability_by_trial,
  file.path(z_out, "sobol_grid_stability_by_trial.csv")
)

write_csv(
  covariate_summary,
  file.path(z_out, "sobol_grid_covariate_summaries.csv")
)

write_csv(
  primary_validation,
  file.path(z_out, "sobol_grid_stability_primary_validation.csv")
)

# Save full object
saveRDS(
  list(
    stability = stability_table,
    stability_by_trial = stability_by_trial,
    covariate_summary = covariate_summary,
    primary_validation = primary_validation,
    n_citl_grid = n_citl_grid,
    grids = grid_objects,
    note = paste(
      "n_sobol = 64 uses saved primary Sobol and CITL outputs.",
      "Larger grids use trial-specific baseline distributions and the",
      "shared ecological Spearman correlation matrix."
    )
  ),
  file.path(z_out, "sobol_grid_stability.rds")
)

# Log --------------------------------------------------------------------------

# Get primary comparison to largest grid
mat_change <- stability_table %>% 
  filter(n_sobol == imp_sobol$J)

# Log lines
log_lines <- c(
  "Sobol grid stability sensitivity analysis",
  "",
  paste0("Grid sizes: ", paste(n_sobol_grid, collapse = ", ")),
  "n_sobol = 64 uses saved primary Sobol and CITL outputs.",
  "Larger grids use trial-specific baseline distributions.",
  "The ecological Spearman correlation matrix is shared across trials.",
  "No design-stage or conduct-stage model refitting was performed.",
  "",
  "Outputs produced:",
  file.path(z_out, "sobol_grid_stability.csv"),
  file.path(z_out, "sobol_grid_stability_by_trial.csv"),
  file.path(z_out, "sobol_grid_covariate_summaries.csv"),
  file.path(z_out, "sobol_grid_stability_primary_validation.csv"),
  file.path(z_out, "sobol_grid_stability.rds"),
  "",
  paste0(
    "Primary 64-grid mean absolute pre-CITL risk difference versus n = 512: ",
    signif(mat_change$mean_abs_difference_pre_citl_risk_vs_n512, 4)
  ),
  paste0(
    "Primary 64-grid maximum absolute pre-CITL risk difference versus n = 512: ",
    signif(mat_change$max_abs_difference_pre_citl_risk_vs_n512, 4)
  ),
  paste0(
    "Primary 64-grid mean absolute CITL shift difference versus n = 512: ",
    signif(mat_change$mean_abs_difference_citl_intercept_shift_vs_n512, 4)
  ),
  paste0(
    "Primary 64-grid maximum absolute CITL shift difference versus n = 512: ",
    signif(mat_change$max_abs_difference_citl_intercept_shift_vs_n512, 4)
  )
)

writeLines(
  log_lines,
  file.path(z_out, "sobol_grid_stability_README.txt")
)
