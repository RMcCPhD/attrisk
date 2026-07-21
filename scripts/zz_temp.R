
# Complete missing baseline summaries ------------------------------------------

# Helper functions to pool numeric and categorical
pool_mu_d <- function(mu, d, w) {
  
  keep <- is.finite(mu) & is.finite(d) & is.finite(w) & w > 0
  
  if (!any(keep)) {
    return(tibble(mu_ref = NA_real_, d_ref = NA_real_))
  }
  
  mu <- mu[keep]
  d <- d[keep]
  w <- w[keep]
  
  mu_pool <- sum(w * mu) / sum(w)
  
  d_pool <- sqrt(
    pmax(
      sum(w * (d^2 + mu^2)) / sum(w) - mu_pool^2,
      0
    )
  )
  
  tibble(mu_ref = mu_pool, d_ref = d_pool)
}

pool_prop <- function(p, w) {
  
  keep <- is.finite(p) & is.finite(w) & w > 0
  
  if (!any(keep)) return(NA_real_)
  
  weighted.mean(p[keep], w[keep])
}

# Helper functions to fill missingness
fill_num_var <- function(data, var, weight = "n") {
  
  mu_nm <- paste0(var, "_mu")
  d_nm <- paste0(var, "_d")
  flag_nm <- paste0(var, "_imputed")
  
  overall_ref <- pool_mu_d(
    mu = data[[mu_nm]],
    d = data[[d_nm]],
    w = data[[weight]]
  )
  
  trial_ref <- data %>%
    group_by(nct_id) %>%
    reframe(
      pool_mu_d(
        mu = .data[[mu_nm]],
        d = .data[[d_nm]],
        w = .data[[weight]]
      )
    ) %>%
    rename(
      "{mu_nm}_trial" := mu_ref,
      "{d_nm}_trial" := d_ref
    )
  
  trt_ref <- data %>%
    group_by(trttype) %>%
    reframe(
      pool_mu_d(
        mu = .data[[mu_nm]],
        d = .data[[d_nm]],
        w = .data[[weight]]
      )
    ) %>%
    rename(
      "{mu_nm}_trt" := mu_ref,
      "{d_nm}_trt" := d_ref
    )
  
  data %>%
    mutate(
      "{flag_nm}" := is.na(.data[[mu_nm]]) | is.na(.data[[d_nm]])
    ) %>%
    left_join(trial_ref, by = "nct_id") %>%
    left_join(trt_ref, by = c("trttype")) %>%
    mutate(
      "{mu_nm}" := coalesce(
        .data[[mu_nm]],
        .data[[paste0(mu_nm, "_trial")]],
        .data[[paste0(mu_nm, "_trt")]],
        overall_ref$mu_ref
      ),
      "{d_nm}" := coalesce(
        .data[[d_nm]],
        .data[[paste0(d_nm, "_trial")]],
        .data[[paste0(d_nm, "_trt")]],
        overall_ref$d_ref
      )
    ) %>%
    select(
      -any_of(c(
        paste0(mu_nm, "_trial"),
        paste0(d_nm, "_trial"),
        paste0(mu_nm, "_trt"),
        paste0(d_nm, "_trt")
      ))
    )
}

fill_prop_var <- function(data, var, weight = "n") {
  
  flag_nm <- paste0(var, "_imputed")
  
  overall_ref <- pool_prop(
    p = data[[var]],
    w = data[[weight]]
  )
  
  trial_ref <- data %>%
    group_by(nct_id) %>%
    summarise(
      "{var}_trial" := pool_prop(
        p = .data[[var]],
        w = .data[[weight]]
      ),
      .groups = "drop"
    )
  
  trt_ref <- data %>%
    group_by(trttype) %>%
    summarise(
      "{var}_trt" := pool_prop(
        p = .data[[var]],
        w = .data[[weight]]
      ),
      .groups = "drop"
    )
  
  data %>%
    mutate(
      "{flag_nm}" := is.na(.data[[var]])
    ) %>%
    left_join(trial_ref, by = "nct_id") %>%
    left_join(trt_ref, by = c("trttype")) %>%
    mutate(
      "{var}" := coalesce(
        .data[[var]],
        .data[[paste0(var, "_trial")]],
        .data[[paste0(var, "_trt")]],
        overall_ref
      ),
      "{var}" := pmin(pmax(.data[[var]], 0), 1)
    ) %>%
    select(
      -any_of(c(
        paste0(var, "_trial"),
        paste0(var, "_trt")
      ))
    )
}

# Helper function to address missingness
complete_bl_summaries <- function(data, num_vars, prop_vars, weight = "n") {
  
  out <- data
  
  for (v in num_vars) {
    out <- fill_num_var(out, v, weight = weight)
  }
  
  for (v in prop_vars) {
    out <- fill_prop_var(out, v, weight = weight)
  }
  
  # Renormalise race proportions
  out <- out %>%
    mutate(
      race_total = white + black_aa + asian + other,
      across(
        c(white, black_aa, asian, other),
        ~ if_else(race_total > 0, .x / race_total, .x)
      )
    ) %>%
    select(-race_total)
  
  out
}

# Apply function
compl_bl <- complete_bl_summaries(
  data = cmbn_all,
  num_vars = num_vars,
  prop_vars = prop_vars,
  weight = "n"
)

# No missingness remaining
check_nas(compl_bl) %>% print(n = Inf)

# Add comorbidities and concomitant medications as independent counts
counts_sum <- imp_count_sum %>% 
  mutate(med = parse_number(value)) %>% 
  select(name, med)

compl_all <- compl_bl %>% 
  mutate(
    comorbs_mu = counts_sum$med[counts_sum$name == "comorbs"],
    conmeds_mu = counts_sum$med[counts_sum$name == "conmeds"],
    comorbs_d = sqrt(comorbs_mu),
    conmeds_d = sqrt(conmeds_mu)
  )

# Double-check all AttRisk covariates are now present
glimpse(compl_all)
imp_covars

# Save for manuscript summaries
write_csv(compl_all, "outputs/manuscript_summaries/ipd_bl_sum.csv")

# Sobol generation requires distribution summaries and correlation matrices
# Since correlation matrices could not be exported (risk of deidentification),
# ecological approximation is used to generate (using covariate means and 
# arm size)

# Prepare single-row Sobol distribution summary --------------------------------

# Numeric variables are pooled using mixture moments across all completed
# IPD-derived arms. Categorical variables are pooled using arm-size weights.
sobol_num <- tibble(var = num_vars) %>% 
  mutate(
    pool = map(
      var,
      ~ pool_mu_d(
        mu = compl_all[[paste0(.x, "_mu")]],
        d = compl_all[[paste0(.x, "_d")]],
        w = compl_all$n
      )
    )
  ) %>% 
  unnest(pool) %>% 
  rename(mu = mu_ref, d = d_ref) %>% 
  pivot_wider(
    names_from = "var",
    values_from = c("mu", "d"),
    names_glue = "{var}_{.value}"
  )

sobol_prop <- tibble(var = prop_vars) %>% 
  mutate(
    value = map_dbl(
      var,
      ~ pool_prop(
        p = compl_all[[.x]],
        w = compl_all$n
      )
    )
  ) %>% 
  pivot_wider(
    names_from = "var",
    values_from = "value"
  )

sobol_dist <- tibble(n = sum(compl_all$n, na.rm = TRUE)) %>% 
  bind_cols(sobol_num, sobol_prop) %>% 
  mutate(
    race_total = white + black_aa + asian + other,
    across(
      c(white, black_aa, asian, other),
      ~ if_else(race_total > 0, .x / race_total, .x)
    ),
    comorbs_mu = counts_sum$med[counts_sum$name == "comorbs"],
    conmeds_mu = counts_sum$med[counts_sum$name == "conmeds"],
    comorbs_d = sqrt(comorbs_mu),
    conmeds_d = sqrt(conmeds_mu)
  ) %>% 
  select(-race_total)

# One row only, for generating 64 Sobol participants in add_integration()
nrow(sobol_dist)
glimpse(sobol_dist)
check_nas(sobol_dist) %>% print(n = Inf)

# Estimate correlation matrix at arm-level -------------------------------------

# Could not export covariate correlation because this risks deidentifying
# IPD trial participants, and a correlation matrix is required to generate
# Sobol sequences using add_integration() from multinma
# Solution is to generate a weighted Spearman ecological correlation matrix 

# Define covariates
cor_vars <- c(
  "age",
  "female_p",
  "bmi",
  "durn",
  "sbp",
  "chol_hdl",
  "chol_total",
  "hba1c",
  "egfr",
  "cursmoke",
  "asian",
  "other",
  "black_aa"
)

# Counts (comorbidities, conmeds) treated as independent because a summary
# of the distribution was not extracted
indep_count_vars <- c("comorbs", "conmeds")

# Transmute summaries
cor_dat <- compl_all %>%
  transmute(
    nct_id,
    n,
    age = age_mu,
    bmi = bmi_mu,
    chol_hdl = chol_hdl_mu,
    chol_total = chol_total_mu,
    durn = durn_mu,
    hba1c = hba1c_mu,
    egfr = egfr_mu,
    sbp = sbp_mu,
    female_p = female_p,
    cursmoke = cursmoke,
    black_aa = black_aa,
    asian = asian,
    other = other
  )

# Weighted Spearman correlation
weighted_cor <- function(x, y, w) {
  
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]
  
  if (length(x) < 3) return(NA_real_)
  
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

weighted_cor_mat <- function(data, vars, weight = "n", method = "spearman") {
  
  x <- data %>%
    select(all_of(vars))
  
  if (method == "spearman") {
    x <- x %>%
      mutate(
        across(
          everything(),
          ~ rank(.x, na.last = "keep", ties.method = "average")
        )
      )
  }
  
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

# Convert Spearman correlation to Gaussian-copula Pearson-scale correlation
spearman_to_copula_cor <- function(cor_mat) {
  2 * sin(pi * cor_mat / 6)
}

# Ensure correlation matrix is valid and positive definite
make_pd_copula_cor <- function(cor_mat, posd.tol = 1e-6) {
  
  cor_mat <- as.matrix(cor_mat)
  cor_mat[!is.finite(cor_mat)] <- 0
  cor_mat <- (cor_mat + t(cor_mat)) / 2
  diag(cor_mat) <- 1
  
  # Convert Spearman target correlations to Gaussian-copula scale
  copula_cor <- spearman_to_copula_cor(cor_mat)
  copula_cor <- (copula_cor + t(copula_cor)) / 2
  diag(copula_cor) <- 1
  
  # Force positive definiteness on the matrix add_integration will actually use
  out <- Matrix::nearPD(
    copula_cor,
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

eco_cor_spearman <- weighted_cor_mat(
  data = cor_dat,
  vars = cor_vars,
  weight = "n",
  method = "spearman"
)

# Append comorbidities and concomitant medications as independent dimensions
add_independent_cor_vars <- function(cor_mat, new_vars) {
  
  old_vars <- colnames(cor_mat)
  all_vars <- c(old_vars, new_vars)
  
  out <- matrix(
    0,
    nrow = length(all_vars),
    ncol = length(all_vars),
    dimnames = list(all_vars, all_vars)
  )
  
  out[old_vars, old_vars] <- cor_mat
  diag(out) <- 1
  
  out
}

cor_vars_final <- c(cor_vars, indep_count_vars)

cor_add <- tibble(
  out_raw = list(
    add_independent_cor_vars(
      cor_mat = eco_cor_spearman,
      new_vars = indep_count_vars
    )
  )
) %>% 
  mutate(
    out = map(out_raw, make_pd_copula_cor),
    min_eigen = map_dbl(
      out,
      ~ min(eigen(.x, symmetric = TRUE, only.values = TRUE)$values)
    )
  )

# Final checks
# All symmetric
# Diagonals equal 1
map(cor_add$out, isSymmetric)
map(cor_add$out, ~ all(abs(diag(.x) - 1) < 1e-8))
map_lgl(
  cor_add$out, 
  ~ all(eigen(.x, symmetric = TRUE, only.values = TRUE)$values > 0)
)

# Correlation matrix diagnostics for supplementary reporting ------------------

dir.create("outputs/conduct/sobol", recursive = TRUE, showWarnings = FALSE)

cor_copula_raw <- map(cor_add$out_raw, ~ {
  x <- as.matrix(.x)
  x[!is.finite(x)] <- 0
  x <- (x + t(x)) / 2
  diag(x) <- 1
  x <- spearman_to_copula_cor(x)
  x <- (x + t(x)) / 2
  diag(x) <- 1
  dimnames(x) <- dimnames(.x)
  x
})

cor_matrix_diagnostics <- cor_add %>%
  mutate(
    out_copula_raw = cor_copula_raw,
    min_eigen_raw_copula = map_dbl(
      out_copula_raw,
      ~ min(eigen(.x, symmetric = TRUE, only.values = TRUE)$values)
    ),
    min_eigen_pd = map_dbl(
      out,
      ~ min(eigen(.x, symmetric = TRUE, only.values = TRUE)$values)
    ),
    max_abs_pd_adjustment = map2_dbl(
      out_copula_raw,
      out,
      ~ max(abs(.y - .x), na.rm = TRUE)
    ),
    mean_abs_pd_adjustment = map2_dbl(
      out_copula_raw,
      out,
      ~ mean(abs(.y - .x), na.rm = TRUE)
    ),
    n_negative_raw_eigenvalues = map_int(
      out_copula_raw,
      ~ sum(eigen(.x, symmetric = TRUE, only.values = TRUE)$values <= 0)
    )
  ) %>%
  transmute(
    matrix = row_number(),
    n_variables = map_int(out, ncol),
    min_eigen_raw_copula,
    min_eigen_pd,
    max_abs_pd_adjustment,
    mean_abs_pd_adjustment,
    n_negative_raw_eigenvalues
  )

write_csv(
  cor_matrix_diagnostics,
  "outputs/conduct/sobol/sobol_correlation_matrix_diagnostics.csv"
)

# Save outputs -----------------------------------------------------------------

# Covariate names order
saveRDS(
  colnames(cor_add$out[[1]]),
  "created_metadata/sobol_covariates_order.rds"
)

# Baseline covariate summaries
saveRDS(sobol_dist, "processed_data/sobol_dist.rds")

# Correlation matrices
saveRDS(cor_add, "processed_data/sobol_cor_mat.rds")































# Complete missing baseline summaries ------------------------------------------

# Pool means and SDs using mixture moments
pool_mu_d <- function(mu, d, w) {
  keep <- is.finite(mu) & is.finite(d) & is.finite(w) & w > 0
  
  if (!any(keep)) {
    return(tibble(mu_ref = NA_real_, d_ref = NA_real_))
  }
  
  mu <- mu[keep]
  d <- d[keep]
  w <- w[keep]
  
  mu_pool <- sum(w * mu) / sum(w)
  d_pool <- sqrt(
    pmax(
      sum(w * (d^2 + mu^2)) / sum(w) - mu_pool^2,
      0
    )
  )
  
  tibble(mu_ref = mu_pool, d_ref = d_pool)
}

# Pool proportions using trial-size weights
pool_prop <- function(p, w) {
  keep <- is.finite(p) & is.finite(w) & w > 0
  
  if (!any(keep)) return(NA_real_)
  
  weighted.mean(p[keep], w[keep])
}

# Fill a mean/SD pair from the treatment-type reference, then the overall
# reference. The flag records whether either source summary was incomplete.
fill_num_var <- function(data, var, weight = "n") {
  mu_nm <- paste0(var, "_mu")
  d_nm <- paste0(var, "_d")
  flag_nm <- paste0(var, "_imputed")
  
  overall_ref <- pool_mu_d(
    mu = data[[mu_nm]],
    d = data[[d_nm]],
    w = data[[weight]]
  )
  
  trt_ref <- data %>%
    group_by(trttype) %>%
    reframe(
      pool_mu_d(
        mu = .data[[mu_nm]],
        d = .data[[d_nm]],
        w = .data[[weight]]
      )
    ) %>%
    rename(
      "{mu_nm}_trt" := mu_ref,
      "{d_nm}_trt" := d_ref
    )
  
  data %>%
    mutate(
      "{flag_nm}" := !is.finite(.data[[mu_nm]]) |
        !is.finite(.data[[d_nm]])
    ) %>%
    left_join(trt_ref, by = "trttype") %>%
    mutate(
      "{mu_nm}" := coalesce(
        .data[[mu_nm]],
        .data[[paste0(mu_nm, "_trt")]],
        overall_ref$mu_ref
      ),
      "{d_nm}" := coalesce(
        .data[[d_nm]],
        .data[[paste0(d_nm, "_trt")]],
        overall_ref$d_ref
      )
    ) %>%
    select(-all_of(c(paste0(mu_nm, "_trt"), paste0(d_nm, "_trt"))))
}

# Fill a proportion from the treatment-type reference, then the overall
# reference, and constrain the completed value to the probability scale.
fill_prop_var <- function(data, var, weight = "n") {
  flag_nm <- paste0(var, "_imputed")
  
  overall_ref <- pool_prop(
    p = data[[var]],
    w = data[[weight]]
  )
  
  trt_ref <- data %>%
    group_by(trttype) %>%
    summarise(
      "{var}_trt" := pool_prop(
        p = .data[[var]],
        w = .data[[weight]]
      ),
      .groups = "drop"
    )
  
  data %>%
    mutate("{flag_nm}" := !is.finite(.data[[var]])) %>%
    left_join(trt_ref, by = "trttype") %>%
    mutate(
      "{var}" := coalesce(
        .data[[var]],
        .data[[paste0(var, "_trt")]],
        overall_ref
      ),
      "{var}" := pmin(pmax(.data[[var]], 0), 1)
    ) %>%
    select(-all_of(paste0(var, "_trt")))
}

complete_bl_summaries <- function(data, num_vars, prop_vars, weight = "n") {
  out <- data
  
  for (v in num_vars) {
    out <- fill_num_var(out, v, weight = weight)
  }
  
  for (v in prop_vars) {
    out <- fill_prop_var(out, v, weight = weight)
  }
  
  # Ensure the mutually exclusive race proportions sum to one
  out %>%
    mutate(
      race_total = white + black_aa + asian + other,
      across(
        c(white, black_aa, asian, other),
        ~ if_else(race_total > 0, .x / race_total, .x)
      )
    ) %>%
    select(-race_total)
}

compl_bl <- complete_bl_summaries(
  data = agd_tidy,
  num_vars = num_vars,
  prop_vars = prop_vars,
  weight = "n"
)

# Comorbidity and concomitant-medication counts were not available for AgD.
# Use the same independent count distributions used in the IPD summary.
count_vars <- c("comorbs_mu", "conmeds_mu", "comorbs_d", "conmeds_d")
count_ref <- imp_bl_ipd %>%
  select(all_of(count_vars)) %>%
  distinct()

stopifnot(nrow(count_ref) == 1L, !anyNA(count_ref))

compl_all <- compl_bl %>%
  mutate(
    comorbs_mu = count_ref$comorbs_mu,
    conmeds_mu = count_ref$conmeds_mu,
    comorbs_d = count_ref$comorbs_d,
    conmeds_d = count_ref$conmeds_d
  )

# Match the IPD column names and order exactly
agd_bl_sum <- compl_all %>%
  rename(nct_id = trial_id) %>%
  select(all_of(ipd_schema)) %>%
  arrange(match(nct_id, imp_agd$trial_id))

# Final integrity checks
race_total <- agd_bl_sum %>%
  transmute(total = white + black_aa + asian + other) %>%
  pull(total)

stopifnot(
  identical(names(agd_bl_sum), ipd_schema),
  nrow(agd_bl_sum) == nrow(imp_agd),
  n_distinct(agd_bl_sum$nct_id) == nrow(imp_agd),
  setequal(agd_bl_sum$nct_id, imp_agd$trial_id),
  !anyNA(agd_bl_sum),
  all(abs(race_total - 1) < 1e-8)
)

check_nas(agd_bl_sum) %>% print(n = Inf)
glimpse(agd_bl_sum)

# Save for manuscript summaries
write_csv(
  agd_bl_sum,
  "outputs/manuscript_summaries/agd_bl_sum.csv"
)
