# Sobol integration stability check for conduct-stage construction
# Regenerates alternative Sobol grids and recomputes participant-marginalised
# conduct-stage summaries without changing the primary 64-point analysis outputs

rm(list = ls())
gc()

source("scripts/00_config.R")

# Optional local package library used when running from a clean command-line R
if (dir.exists("r_libs")) {
  .libPaths(c(normalizePath("r_libs"), .libPaths()))
}

if (!requireNamespace("randtoolbox", quietly = TRUE)) {
  stop(
    "Package 'randtoolbox' is required for this follow-on Sobol stability ",
    "script. Install it or run install.packages('randtoolbox')."
  )
}

# Inputs and outputs -----------------------------------------------------------

z_out_dir <- "outputs/sensitivity"

dir.create(z_out_dir, recursive = TRUE, showWarnings = FALSE)

z_sobol_dist <- "processed_data/sobol_dist.rds"
z_sobol_cor <- "processed_data/sobol_cor_mat.rds"
z_sobol_primary <- "processed_data/sobol_lst.rds"
z_sobol_primary_profiles <- "processed_data/sobol_unnest.rds"
z_scale <- "vivli/v2/Summaries/covariate_mu_d.csv"
z_citl <- "outputs/conduct/citl/citl_draws.rds"

for (z in c(
  z_sobol_dist,
  z_sobol_cor,
  z_sobol_primary,
  z_sobol_primary_profiles,
  z_scale,
  z_citl
)) {
  if (!file.exists(z)) stop("Missing required input: ", z)
}

imp_dist <- readRDS(z_sobol_dist)
imp_cor <- readRDS(z_sobol_cor)
imp_sobol_primary <- readRDS(z_sobol_primary)
imp_sobol_primary_profiles <- readRDS(z_sobol_primary_profiles)
imp_scale <- read.csv(z_scale, stringsAsFactors = FALSE)
imp_citl <- readRDS(z_citl)

n_sobol_grid <- c(64L, 128L, 256L, 512L)

# Helper functions -------------------------------------------------------------

q025 <- function(x) {
  quantile(x, 0.025, na.rm = TRUE, names = FALSE)
}

q975 <- function(x) {
  quantile(x, 0.975, na.rm = TRUE, names = FALSE)
}

summ_vec <- function(x) {
  data.frame(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = q025(x),
    q975 = q975(x)
  )
}

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
  
  vapply(
    cv,
    function(x) {
      uniroot(
        function(k) {
          sqrt(gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1) - x
        },
        interval = c(0.1, 50)
      )$root
    },
    numeric(1)
  )
}

get_weibull_scale <- function(mu, shape) {
  mu / gamma(1 + 1 / shape)
}

scale_centre_cov <- function(dat, imp_scale) {
  
  out <- dat
  scale_names <- intersect(imp_scale$name, names(out))
  
  for (nm in scale_names) {
    idx <- match(nm, imp_scale$name)
    out[[nm]] <- (out[[nm]] - imp_scale$mu[idx]) / imp_scale$d[idx]
  }
  
  out
}

assign_race_indicators <- function(dat, p_asian, p_black_aa, p_other, shift = 0) {
  
  p_asian <- max(as.numeric(p_asian), 0, na.rm = TRUE)
  p_black_aa <- max(as.numeric(p_black_aa), 0, na.rm = TRUE)
  p_other <- max(as.numeric(p_other), 0, na.rm = TRUE)
  p_white <- max(1 - p_asian - p_black_aa - p_other, 0)
  
  p_cat <- c(
    asian = p_asian,
    black_aa = p_black_aa,
    other = p_other,
    white = p_white
  )
  
  p_cat <- p_cat / sum(p_cat)
  cuts <- cumsum(p_cat)
  
  n <- nrow(dat)
  race_u <- (((seq_len(n) - 0.5) / n) + shift) %% 1
  
  race <- ifelse(
    race_u < cuts[["asian"]],
    "asian",
    ifelse(
      race_u < cuts[["black_aa"]],
      "black_aa",
      ifelse(race_u < cuts[["other"]], "other", "white")
    )
  )
  
  dat$asian <- as.integer(race == "asian")
  dat$black_aa <- as.integer(race == "black_aa")
  dat$other <- as.integer(race == "other")
  
  dat
}

make_distribution_parameters <- function(dist) {
  
  dist <- as.data.frame(dist)
  
  d_cols <- grep("_d$", names(dist), value = TRUE)
  dist[d_cols] <- lapply(dist[d_cols], function(x) pmax(x, .Machine$double.eps))
  
  dist$bmi_meanlog <- get_lnorm_meanlog(dist$bmi_mu, dist$bmi_d)
  dist$bmi_sdlog <- get_lnorm_sdlog(dist$bmi_mu, dist$bmi_d)
  
  dist$chol_hdl_meanlog <- get_lnorm_meanlog(dist$chol_hdl_mu, dist$chol_hdl_d)
  dist$chol_hdl_sdlog <- get_lnorm_sdlog(dist$chol_hdl_mu, dist$chol_hdl_d)
  
  dist$sbp_meanlog <- get_lnorm_meanlog(dist$sbp_mu, dist$sbp_d)
  dist$sbp_sdlog <- get_lnorm_sdlog(dist$sbp_mu, dist$sbp_d)
  
  dist$chol_total_shape <- get_gamma_shape(dist$chol_total_mu, dist$chol_total_d)
  dist$chol_total_rate <- get_gamma_rate(dist$chol_total_mu, dist$chol_total_d)
  
  dist$durn_shape <- get_weibull_shape(dist$durn_mu, dist$durn_d)
  dist$durn_scale <- get_weibull_scale(dist$durn_mu, dist$durn_shape)
  
  dist$egfr_shape <- get_weibull_shape(dist$egfr_mu, dist$egfr_d)
  dist$egfr_scale <- get_weibull_scale(dist$egfr_mu, dist$egfr_shape)
  
  dist
}

safe_chol <- function(cor_mat) {
  
  out <- tryCatch(
    chol(cor_mat),
    error = function(e) NULL
  )
  
  if (!is.null(out)) return(out)
  
  chol(as.matrix(Matrix::nearPD(cor_mat, corr = TRUE)$mat))
}

generate_sobol_profiles <- function(n_sobol, dist, cor_mat, imp_scale, loghr) {
  
  race_cols <- c("asian", "black_aa", "other")
  keep_cols <- setdiff(rownames(cor_mat), race_cols)
  cor_no_race <- cor_mat[keep_cols, keep_cols, drop = FALSE]
  
  u <- randtoolbox::sobol(
    n = n_sobol,
    dim = ncol(cor_no_race),
    scrambling = 0,
    seed = 1,
    init = TRUE,
    normal = FALSE
  )
  
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  z <- qnorm(u) %*% safe_chol(cor_no_race)
  u_cor <- pnorm(z)
  colnames(u_cor) <- colnames(cor_no_race)
  
  get_u <- function(name) {
    if (!name %in% colnames(u_cor)) {
      stop("Correlation matrix column not found: ", name)
    }
    u_cor[, name]
  }
  
  raw <- data.frame(
    age = qnorm(get_u("age"), mean = dist$age_mu, sd = dist$age_d),
    male_p = qbinom(get_u("female_p"), size = 1, prob = 1 - dist$female_p),
    bmi = qlnorm(get_u("bmi"), meanlog = dist$bmi_meanlog, sdlog = dist$bmi_sdlog),
    durn = qweibull(get_u("durn"), shape = dist$durn_shape, scale = dist$durn_scale),
    sbp = qlnorm(get_u("sbp"), meanlog = dist$sbp_meanlog, sdlog = dist$sbp_sdlog),
    chol_hdl = qlnorm(
      get_u("chol_hdl"),
      meanlog = dist$chol_hdl_meanlog,
      sdlog = dist$chol_hdl_sdlog
    ),
    chol_total = qgamma(
      get_u("chol_total"),
      shape = dist$chol_total_shape,
      rate = dist$chol_total_rate
    ),
    hba1c = qnorm(get_u("hba1c"), mean = dist$hba1c_mu, sd = dist$hba1c_d),
    egfr = qweibull(get_u("egfr"), shape = dist$egfr_shape, scale = dist$egfr_scale),
    cursmoke = qbinom(get_u("cursmoke"), size = 1, prob = dist$cursmoke),
    comorbs = qpois(get_u("comorbs"), lambda = dist$comorbs_mu),
    conmeds = qpois(get_u("conmeds"), lambda = dist$conmeds_mu)
  )
  
  raw <- assign_race_indicators(
    dat = raw,
    p_asian = dist$asian,
    p_black_aa = dist$black_aa,
    p_other = dist$other
  )
  
  scaled <- raw
  names(scaled)[names(scaled) == "age"] <- "age10"
  names(scaled)[names(scaled) == "male_p"] <- "sex"
  scaled$age10 <- scaled$age10 / 10
  scaled <- scale_centre_cov(scaled, imp_scale)
  
  Z <- as.matrix(scaled[, names(loghr), drop = FALSE])
  lp <- as.vector(Z %*% loghr)
  
  list(
    raw = raw,
    scaled = scaled,
    Z = Z,
    lp = lp,
    hr = exp(lp)
  )
}

summarise_covariates <- function(dat, n_sobol, scale) {
  
  nums <- names(dat)[vapply(dat, is.numeric, logical(1))]
  
  do.call(
    rbind,
    lapply(
      nums,
      function(nm) {
        x <- dat[[nm]]
        data.frame(
          n_sobol = n_sobol,
          scale = scale,
          covariate = nm,
          mean = mean(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          sd = sd(x, na.rm = TRUE),
          q025 = q025(x),
          q975 = q975(x)
        )
      }
    )
  )
}

laplace_survival <- function(h, hr) {
  rowMeans(exp(-outer(h, hr)))
}

make_h_grid <- function(H, n_grid = 4096L) {
  
  H <- H[is.finite(H) & H > 0]
  rng <- range(H)
  
  exp(seq(log(rng[1]), log(rng[2]), length.out = n_grid))
}

make_citl_functions <- function(h_grid, hr, n_iter = 48L) {
  
  pre_grid <- 1 - laplace_survival(h_grid, hr)
  
  target_survival <- exp(-h_grid)
  lo <- rep(-30, length(h_grid))
  hi <- rep(30, length(h_grid))
  
  for (iter in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    survival_mid <- laplace_survival(h_grid * exp(mid), hr)
    f_mid <- survival_mid - target_survival
    
    lo[f_mid > 0] <- mid[f_mid > 0]
    hi[f_mid <= 0] <- mid[f_mid <= 0]
  }
  
  delta_grid <- (lo + hi) / 2
  
  list(
    pre = approxfun(h_grid, pre_grid, rule = 2),
    delta = approxfun(h_grid, delta_grid, rule = 2)
  )
}

get_primary_hr <- function(imp_sobol) {
  
  J <- imp_sobol$J
  Z <- imp_sobol$Z[seq_len(J), names(imp_sobol$loghr), drop = FALSE]
  
  as.vector(exp(Z %*% imp_sobol$loghr))
}

get_delta_matrix <- function(delta_vec, template) {
  matrix(delta_vec, nrow = nrow(template), ncol = ncol(template))
}

# Prepare conduct-stage inputs -------------------------------------------------

citl_beta <- imp_sobol_primary$loghr

dist_par <- make_distribution_parameters(imp_dist)
cor_mat <- imp_cor$out[[1]]

target_risk <- imp_citl$target_risk
H_mat <- -log1p(-pmin(pmax(target_risk, 0), 1 - 1e-12))
H_vec <- as.vector(H_mat)
h_grid <- make_h_grid(H_vec)

# Generate grids and CITL stability summaries ---------------------------------

grid_objects <- list()
covariate_summaries <- list()
stability_work <- list()

for (n_sobol in n_sobol_grid) {
  
  message(sprintf("Processing n_sobol = %d", n_sobol))
  
  if (n_sobol == 64L) {
    hr <- get_primary_hr(imp_sobol_primary)
    lp <- log(hr)
    
    pre_vec <- as.vector(imp_citl$avg_risk_pre)
    delta_vec <- as.vector(imp_citl$delta)
    
    grid_objects[[as.character(n_sobol)]] <- list(
      n_sobol = n_sobol,
      source = "primary_saved_grid",
      raw = imp_sobol_primary_profiles$sobol_unnest[[1]],
      scaled = imp_sobol_primary_profiles$sobol_rdy[[1]],
      lp = lp,
      hr = hr
    )
    
    covariate_summaries[[as.character(n_sobol)]] <- rbind(
      summarise_covariates(
        imp_sobol_primary_profiles$sobol_unnest[[1]],
        n_sobol,
        "raw"
      ),
      summarise_covariates(
        imp_sobol_primary_profiles$sobol_rdy[[1]],
        n_sobol,
        "scaled"
      )
    )
  } else {
    grid <- generate_sobol_profiles(
      n_sobol = n_sobol,
      dist = dist_par,
      cor_mat = cor_mat,
      imp_scale = imp_scale,
      loghr = citl_beta
    )
    
    citl_fun <- make_citl_functions(h_grid, grid$hr)
    
    pre_vec <- citl_fun$pre(H_vec)
    delta_vec <- citl_fun$delta(H_vec)
    
    covariate_summaries[[as.character(n_sobol)]] <- rbind(
      summarise_covariates(grid$raw, n_sobol, "raw"),
      summarise_covariates(grid$scaled, n_sobol, "scaled")
    )
    
    grid_objects[[as.character(n_sobol)]] <- list(
      n_sobol = n_sobol,
      source = "generated_grid",
      raw = grid$raw,
      scaled = grid$scaled,
      lp = grid$lp,
      hr = grid$hr
    )
  }
  
  stability_work[[as.character(n_sobol)]] <- list(
    n_sobol = n_sobol,
    lp = grid_objects[[as.character(n_sobol)]]$lp,
    hr = grid_objects[[as.character(n_sobol)]]$hr,
    pre = pre_vec,
    delta = delta_vec
  )
}

ref <- stability_work[["512"]]

stability_table <- do.call(
  rbind,
  lapply(
    stability_work,
    function(x) {
      lp_sum <- summ_vec(x$lp)
      hr_sum <- summ_vec(x$hr)
      
      pre_diff <- abs(x$pre - ref$pre)
      delta_diff <- abs(x$delta - ref$delta)
      
      data.frame(
        n_sobol = x$n_sobol,
        mean_participant_linear_predictor = lp_sum$mean,
        sd_participant_linear_predictor = lp_sum$sd,
        q025_participant_linear_predictor = lp_sum$q025,
        q975_participant_linear_predictor = lp_sum$q975,
        mean_hazard_multiplier = hr_sum$mean,
        sd_hazard_multiplier = hr_sum$sd,
        q025_hazard_multiplier = hr_sum$q025,
        q975_hazard_multiplier = hr_sum$q975,
        mean_pre_citl_participant_marginalised_risk = mean(x$pre, na.rm = TRUE),
        mean_post_citl_participant_marginalised_risk = mean(target_risk, na.rm = TRUE),
        mean_citl_intercept_shift = mean(x$delta, na.rm = TRUE),
        mean_abs_difference_pre_citl_risk_vs_n512 = mean(pre_diff, na.rm = TRUE),
        max_abs_difference_pre_citl_risk_vs_n512 = max(pre_diff, na.rm = TRUE),
        mean_abs_difference_citl_intercept_shift_vs_n512 = mean(delta_diff, na.rm = TRUE),
        max_abs_difference_citl_intercept_shift_vs_n512 = max(delta_diff, na.rm = TRUE)
      )
    }
  )
)

stability_table <- stability_table[order(stability_table$n_sobol), ]
covariate_summary <- do.call(rbind, covariate_summaries)

# Validate interpolation against the saved primary 64-grid CITL results --------

primary_fun <- make_citl_functions(h_grid, get_primary_hr(imp_sobol_primary))
primary_pre_check <- primary_fun$pre(H_vec)
primary_delta_check <- primary_fun$delta(H_vec)

primary_validation <- data.frame(
  check = c(
    "max_abs_pre_citl_risk_interpolation_error",
    "mean_abs_pre_citl_risk_interpolation_error",
    "max_abs_delta_interpolation_error",
    "mean_abs_delta_interpolation_error"
  ),
  value = c(
    max(abs(primary_pre_check - as.vector(imp_citl$avg_risk_pre)), na.rm = TRUE),
    mean(abs(primary_pre_check - as.vector(imp_citl$avg_risk_pre)), na.rm = TRUE),
    max(abs(primary_delta_check - as.vector(imp_citl$delta)), na.rm = TRUE),
    mean(abs(primary_delta_check - as.vector(imp_citl$delta)), na.rm = TRUE)
  )
)

# Save outputs -----------------------------------------------------------------

write.csv(
  stability_table,
  file.path(z_out_dir, "sobol_grid_stability.csv"),
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

write.csv(
  covariate_summary,
  file.path(z_out_dir, "sobol_grid_covariate_summaries.csv"),
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

write.csv(
  primary_validation,
  file.path(z_out_dir, "sobol_grid_stability_primary_validation.csv"),
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

saveRDS(
  list(
    stability = stability_table,
    covariate_summary = covariate_summary,
    primary_validation = primary_validation,
    grids = grid_objects,
    note = paste(
      "n_sobol = 64 uses the saved primary grid and saved primary CITL",
      "draws. Larger grids are generated from the same IPD-informed",
      "distribution and ecological Gaussian-copula correlation matrix."
    )
  ),
  file.path(z_out_dir, "sobol_grid_stability.rds")
)

# Log --------------------------------------------------------------------------

mat_change <- stability_table[
  stability_table$n_sobol == 64L,
  c(
    "mean_abs_difference_pre_citl_risk_vs_n512",
    "max_abs_difference_pre_citl_risk_vs_n512",
    "mean_abs_difference_citl_intercept_shift_vs_n512",
    "max_abs_difference_citl_intercept_shift_vs_n512"
  )
]

log_lines <- c(
  "Sobol grid stability sensitivity analysis",
  "",
  paste0("Grid sizes: ", paste(n_sobol_grid, collapse = ", ")),
  "n_sobol = 64 uses the saved primary Sobol and CITL outputs.",
  "n_sobol = 128, 256, and 512 were regenerated from the saved reference distribution and correlation matrix.",
  "No design-stage or conduct-stage model refitting was performed.",
  "",
  "Outputs produced:",
  file.path(z_out_dir, "sobol_grid_stability.csv"),
  file.path(z_out_dir, "sobol_grid_covariate_summaries.csv"),
  file.path(z_out_dir, "sobol_grid_stability_primary_validation.csv"),
  file.path(z_out_dir, "sobol_grid_stability.rds"),
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
  file.path(z_out_dir, "sobol_grid_stability_README.txt")
)
