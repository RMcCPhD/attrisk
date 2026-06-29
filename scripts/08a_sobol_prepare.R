
# Prepare IPD baseline covariate information for generating Sobol sequences
# IPD baseline characteristic distributions pulled from previous export of
# materials for network meta-analysis of treatment effect heterogeneity
# Export materials had distribution data for all covariates except comorbidites
# and concomitant medications
# Sobol generation requires distribution summaries and correlation matrices
# Since correlation matrices could not be exported (risk of deidentification),
# ecological approximation is used to generate (using covariate means and 
# arm size)
# Data prepared as one overall IPD-derived reference distribution so Sobol
# generation produces 64 synthetic participants from a single row

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import comorbidity and conmed summaries
imp_count_sum <- read_csv("vivli/v2/Summaries/overall_summary.csv") %>% 
  filter(name %in% c("comorbs", "conmeds"))

# Import Vivli arm-level summary
imp_arm <- read_csv("vivli/v2/Summaries/trial_level_summary.csv")

# Get covariate names
imp_covars <- read_csv("vivli/v2/Models/export_res.csv") %>% 
  distinct(term)

# Numeric covariates except age
# Distribution-dependent, kept nested to obtain mean and d
files_bl <- list.files(
  "data/ipd_covariates",
  pattern = "covariates", 
  full.names = TRUE
)

imp_bl <- tibble(
  dist = gsub("covariates_|\\.csv", "", basename(files_bl)),
  imp = map(files_bl, read_csv, show_col_types = FALSE)
)

# Age and sex distributions
imp_bl2 <- read_csv("data/ipd_covariates/ipd_age_sex.csv")

# Categorical covariates
imp_bl3 <- read_csv("data/ipd_covariates/categorical_variables_percentage.csv")

# HbA1c
imp_bl4 <- read_csv("data/ipd_covariates/hba1c_base_summary_ipd.csv")

# Import simulated IPD to obtain arm-level counts
imp_n <- readRDS("data/ipd_covariates/agg_ipd_hba1c.Rds") %>% 
  select(ipd) %>% 
  unnest(ipd) %>% 
  ungroup() %>% 
  group_by(nct_id, arm_f, arm_lvl, trtcls5) %>% 
  reframe(total = n())

# Prepare numeric covariates ---------------------------------------------------

# Determine the appropriate distribution to pass to add_integration() for
# each covariate - based on the distribution that most frequently best
# described a covariate
imp_bl %>% 
  unnest(imp) %>% 
  ungroup() %>% 
  count(dist, var) %>% 
  group_by(var) %>% 
  arrange(desc(n), .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup() %>% 
  write_csv("scratch_data/sobol_covariate_dist.csv")

# Get mean and sd for numeric covariates
# Converted depending on the distribution that best described covariates
bl_num <- imp_bl %>% 
  mutate(
    est = map2(imp, dist, ~ {
      if (.y == "gamma") {
        .x %>% 
          mutate(
            mu = est_shape / est_rate,
            d = sqrt(est_shape) / est_rate
          )
      } else if (.y == "lnorm") {
        .x %>% 
          mutate(
            mu = exp(est_meanlog + 0.5 * est_sdlog^2),
            d = sqrt(
              (exp(est_sdlog^2) - 1) * exp(2 * est_meanlog + est_sdlog^2)
            )
          )
      } else if (.y == "norm") {
        .x %>% 
          mutate(
            mu = est_mean,
            d = est_sd
          )
      } else if (.y == "weibull") {
        .x %>% 
          mutate(
            mu = est_scale * gamma(1 + 1 / est_shape),
            d = est_scale * sqrt(
              gamma(1 + 2 / est_shape) - gamma(1 + 1 / est_shape)^2
            )
          )
      }
    })
  ) %>% 
  select(dist, est) %>% 
  unnest(est) %>% 
  ungroup() %>% 
  select(dist, nct_id, arm_f, var, mu, d)

# Prepare age distribution
# Trial-level so assumed to be applicable to arms within trials
bl_age <- bl_num %>% 
  distinct(nct_id, arm_f) %>% 
  left_join(
    imp_bl2 %>% 
      select(nct_id, age_m, age_s) %>% 
      mutate(var = "age", dist = "norm") %>% 
      rename(mu = age_m, d = age_s)
  )

# Combine numeric distributions
bl_num_all <- bl_num %>%
  full_join(bl_age) %>% 
  arrange(nct_id, arm_f, var)

# No missingness
check_nas(bl_num_all)

# Join arm counts from simulated IPD
bl_num_n <- bl_num_all %>% 
  mutate(arm_f = gsub("_dcmplx", "", arm_f)) %>% 
  left_join(
    imp_n %>% 
      mutate(
        arm_f_old = arm_f,
        arm_f = paste0(
          substr(arm_f_old, 2, 3), arm_lvl
        )
      ) %>% 
      select(-arm_f_old) %>% 
      distinct()
  )

check_nas(bl_num_n)

# Get AttRisk IPD trials
bl_num_set <- bl_num_n %>% inner_join(imp_arm %>% distinct(nct_id))
n_distinct(bl_num_set$nct_id)
check_nas(bl_num_set)

# Drop experimental arms
# Done by harmonising covariate distribution arm variable to match
# AttRisk arm-level attrition counts metadata, then inner-join to drop arms
# Manually add arm count for NCT00800683 placebo
bl_num_rm_exp <- bl_num_set %>% 
  mutate(trtcls5 = substr(arm_f, 3, 7)) %>% 
  inner_join(imp_arm %>% distinct(nct_id, trtcls5, trttype, n)) %>% 
  select(dist, nct_id, trtcls5, trttype, total, var, mu, d) %>% 
  distinct() %>% 
  mutate(
    total = if_else(is.na(total) & nct_id == "NCT00800683", 63, total)
  )

# Compare number of arms per trial to verify join was done correctly
# 83 arms in each, join was done correctly
bl_num_set %>% 
  mutate(trtcls5 = substr(arm_f, 3, 7)) %>%
  distinct(nct_id, trtcls5) %>% 
  inner_join(imp_arm %>% distinct(nct_id, trtcls5))

imp_arm %>% distinct(nct_id, trtcls5)

# Address arms with >1 distribution (arms of different doses)
# These were collapsed for AttRisk
bl_num_rm_exp %>% 
  group_by(nct_id, trtcls5, var) %>% 
  filter(n() > 1)

# Join age
bl_num_age <- bl_num_rm_exp %>% 
  full_join(
    bl_num_rm_exp %>% 
      distinct(nct_id, trtcls5, trttype) %>% 
      left_join(
        imp_bl2 %>% 
          mutate(var = "age", dist = "norm") %>% 
          rename(mu = age_m, d = age_s) %>% 
          select(-male_p)
      )
  )

# Join HbA1c
bl_num_hba1c <- imp_bl4 %>% 
  rename(arm_lvl = arm_f) %>% 
  left_join(
    imp_n %>% 
      mutate(
        arm_f_old = arm_f,
        arm_f = paste0(
          substr(arm_f_old, 2, 3), arm_lvl
        )
      ) %>% 
      select(-arm_f_old) %>% 
      distinct()
  ) %>% 
  mutate(
    var = "hba1c", 
    dist = "norm"
  ) %>% 
  inner_join(
    bl_num_age %>% distinct(nct_id, trtcls5, total)
  )

bl_num_hba1c <- bl_num_age %>% 
  full_join(
    bl_num_hba1c %>% 
      rename(mu = hba1cbase_m, d = hba1cbase_s) %>% 
      select(-c(arm_f, arm_lvl))
  ) %>% 
  mutate(
    trttype = case_when(
      is.na(trttype) & trtcls5 == "place" ~ "plc",
      is.na(trttype) & trtcls5 != "place" ~ "ac",
      TRUE ~ trttype
    )
  )

# Multiple arms exist in some trials (e.g. metformin of different doses)
# Arms were collapsed for AttRisk set
# Collapsing distribution summary across these trials
# Use sample size-weighting to collapse
bl_num_clp <- bl_num_hba1c %>% 
  group_by(nct_id) %>% 
  mutate(
    trttype = case_when(
      n_distinct(trttype) == 1 ~ trttype,
      TRUE ~ "plc_ac"
    )
  ) %>% 
  ungroup() %>% 
  group_by(dist, nct_id, trttype, var) %>% 
  reframe(
    arm_total = sum(total),
    mu_pool = sum(total * mu) / arm_total,
    var_pool = sum((total - 1) * d^2 + total * (mu - mu_pool)^2) / (arm_total - 1),
    d_pool = sqrt(var_pool)
  ) %>% 
  select(dist, nct_id, trttype, var, mu_pool, d_pool)

# Join AttRisk arm size
bl_num_rdy <- bl_num_clp %>% 
  left_join(
    imp_arm %>% 
      distinct(nct_id, trtcls5, trttype, n) %>% 
      group_by(nct_id) %>% 
      mutate(
        trttype = case_when(
          n_distinct(trttype) == 1 ~ trttype,
          TRUE ~ "plc_ac"
        )
      ) %>% 
      group_by(nct_id, trttype) %>% 
      reframe(n = sum(n))
  ) %>%
  arrange(nct_id, var)

# No missingness
check_nas(bl_num_rdy)

# Prepare categorical covariates -----------------------------------------------

# Collapse other race categories into one
# Convert percentages to proportions
bl_cat <- imp_bl3 %>% 
  mutate(
    across(-c(nct_id, arm_f, cursmoke), ~ if_else(is.na(.), 0, .)),
    other = ai_an + other
  ) %>% 
  select(
    nct_id, arm_f,
    cursmoke,
    white, black_aa, asian, other
  ) %>% 
  mutate(
    across(-c(nct_id, arm_f), ~ . / 100)
  )

# Verify that race categories sum to 100%
# All do, rounding errors expected (either 0.99 or 1.01)
bl_cat %>% 
  mutate(race_total = white + black_aa + asian + other) %>% 
  filter(race_total != 1) %>% 
  print(n = Inf)

# Join groups for merging with AttRisk metadata
bl_cat_grps <- bl_cat %>% 
  left_join(
    bl_num_n %>% distinct(nct_id, arm_f, trtcls5)
  ) %>% 
  mutate(
    trtcls5 = if_else(
      is.na(trtcls5) 
      & nct_id == "NCT00800683"
      & arm_f == "1_placebo",
      "place",
      trtcls5
    )
  )

# Get AttRisk set
bl_cat_set <- bl_cat_grps %>% 
  inner_join(
    imp_arm %>% distinct(nct_id, trtcls5, trttype)
  )

# Add sex proportion
bl_cat_rdy <- bl_cat_set %>% 
  left_join(
    bl_cat_set %>% 
      distinct(nct_id, trtcls5, trttype) %>% 
      left_join(imp_bl2 %>% select(nct_id, male_p))
  )

# Collapse percentages across trttype within trials
bl_cat_clp <- bl_cat_rdy %>% 
  group_by(nct_id) %>% 
  mutate(
    trttype = case_when(
      n_distinct(trttype) == 1 ~ trttype,
      TRUE ~ "plc_ac"
    )
  ) %>% 
  left_join(bl_num_rdy %>% distinct(nct_id, trttype, n)) %>% 
  group_by(nct_id, trttype) %>% 
  reframe(
    across(
      c(male_p, cursmoke, white, black_aa, asian, other),
      ~ weighted.mean(.x, w = n, na.rm = TRUE)
    )
  )

# Current smoker missing in 7 trials
check_nas(bl_cat_clp)

# Combine covariate distribution summaries and tidy ----------------------------

# Widen numeric covariates summaries, one row per IPD trial arm
cmbn_num_wide <- bl_num_rdy %>% 
  rename(mu = mu_pool, d = d_pool) %>% 
  select(-dist) %>% 
  pivot_wider(
    names_from = "var",
    values_from = c("mu", "d"),
    names_glue = "{var}_{.value}"
  )

# Combine
# Swap male p for female p
# Drop cholesterol ratio, not in AttRisk covariates
cmbn_all <- cmbn_num_wide %>% 
  left_join(bl_cat_clp) %>% 
  mutate(female_p = 1 - male_p) %>% 
  select(-c(male_p, contains("tot_hdl")))

# Confirm all AttRisk covariates are present
# HbA1c, comorbs and conmeds need added
imp_covars
glimpse(cmbn_all)

# Check missingness
# Missingness in durn, cholesterol, egfr, cursmoke, bmi, sbp
# Mostly durn
check_nas(cmbn_all) %>% arrange(desc(value)) %>% print(n = Inf)

# Define covariate sets
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

# Variables as arm-level proportions
prop_vars <- c(
  "cursmoke",
  "female_p",
  "white",
  "black_aa",
  "asian",
  "other"
)

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
