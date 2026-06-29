
# Conduct a random-effects meta-analysis to obtain pooled estimates of
# participant-level effects on attrition risk
# Obtained from the IPD trials, model outputs from Vivli prepared here

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import Vivli outputs
# Model estimates on log-scale
# Covariate correlation, upper triangle of VCOV converted to correlations
imp_vivli_res <- read_csv("vivli/v2/Models/export_res.csv")
imp_vivli_cor <- read_csv("vivli/v2/Models/export_covcor.csv")

# Prepare meta-analysis inputs -------------------------------------------------

# Get unique terms
res_unique_terms <- unique(imp_vivli_res$term)

# Create grid of terms
res_grid_terms <- expand.grid(
  nct_id = unique(imp_vivli_res$nct_id),
  term = res_unique_terms
) %>% as_tibble()

# Fill in missing terms (i.e. covariate or covariate level not fitted in
# trial-level IPD models due to convergence issues)
# Add zero estimates and large standard errors
res_compl_terms <- imp_vivli_res %>%
  full_join(res_grid_terms) %>% 
  mutate(
    estimate = if_else(is.na(estimate), 0, estimate),
    std.error = if_else(is.na(std.error), 99999, std.error)
  ) %>% 
  arrange(nct_id, term)

# Sort terms for VCOV matrices
vcov_terms <- res_unique_terms %>% sort()

# Helper: Complete correlation matrix in long-format for each trial
complete_cor_mat <- function(df, vcov_terms) {
  
  # Canonical pairs
  canon <- df %>% 
    transmute(nct_id, i = pmin(row, col), j = pmax(row, col), value) %>% 
    distinct()
  
  # Get current terms
  obs_vars <- sort(unique(c(df$row, df$col)))
  
  # Include missing terms
  vars <- sort(union(obs_vars, vcov_terms))
  
  # Build grid
  if (length(vars) < 2) {
    return(
      tibble(
        nct_id = unique(df$nct_id),
        row = character(),
        col = character(),
        value = numeric()
      )
    )
  }
  grid <- as_tibble(t(combn(vars, 2)))
  names(grid) <- c("i", "j")
  
  # Left-join observed terms, rename as row/col
  grid %>% 
    mutate(nct_id = unique(df$nct_id)) %>% 
    left_join(canon, by = c("nct_id", "i", "j")) %>% 
    mutate(value = if_else(is.na(value), 0, value)) %>% 
    transmute(nct_id, row = i, col = j, value)
  
}

# Helper: Apply above function to each trial
apply_complete_cor <- function(imp_vivli_cor, vcov_terms) {
  imp_vivli_cor %>% 
    group_by(nct_id) %>% 
    group_split() %>% 
    map_dfr(complete_cor_mat, vcov_terms)
}

# Get complete correlation matrix values
vcov_complete_cor <- apply_complete_cor(imp_vivli_cor, vcov_terms)

# List unique covariate names per trial
vcov_unique_terms <- vcov_complete_cor %>% 
  pivot_longer(row:col, values_to = "variable") %>% 
  distinct(nct_id, variable) %>% 
  group_by(nct_id) %>% 
  reframe(vars = list(unique(variable)))

# Bind correlation and rename components for upper and lower triangles
# Create diagonals (i.e. cor = 1)
vcov_construct_prep <- vcov_complete_cor %>% 
  bind_rows(vcov_complete_cor %>% rename(col = row, row = col)) %>% 
  group_by(nct_id) %>% 
  group_split() %>% 
  map_dfr(function(data_part) {
    key <- data_part %>% select(nct_id) %>% slice(1)
    vars <- unique(c(data_part$row, data_part$col))
    
    # Get diagonals
    diags <- expand.grid(row = vars, col = vars) %>% 
      filter(row == col) %>% 
      mutate(value = 1) %>% 
      bind_cols(key[rep(1, nrow(.)), ])
    
    # Bind diagonals and triangles
    bind_rows(data_part, diags)
  })

# Construct correlation matrices
vcov_construct_cor <- vcov_construct_prep %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  mutate(
    cor_matrix = map(data, ~ {
      vars <- union(.x$row, .x$col)
      cor_tbl <- xtabs(value ~ row + col, data = .x)[vars, vars]
      matrix(
        as.numeric(cor_tbl),
        nrow = length(vars),
        ncol = length(vars),
        dimnames = list(vars, vars)
      )
    })
  ) %>% 
  select(-data)

# Convert correlation to VCOV
vcov_construct <- vcov_construct_cor %>% 
  mutate(
    vcov_matrix = pmap(list(nct_id, cor_matrix), function(id, cor_mat) {
      se_df <- res_compl_terms %>% filter(nct_id == id)
      se_vec <- setNames(se_df$std.error, se_df$term)
      se_order <- se_vec[rownames(cor_mat)]
      cov_mat <- outer(se_order, se_order) * cor_mat
      cov_mat
    })
  )

# Quick check that dimensions are correct for each (16x16)
for (i in seq_along(vcov_construct$vcov_matrix)) {
  print(dim(vcov_construct$vcov_matrix[[i]]))
}

# Save prepared outputs
saveRDS(res_compl_terms, "outputs/conduct/ma/data/ipd_trial_res.rds")
saveRDS(vcov_construct, "outputs/conduct/ma/data/ipd_trial_vcov.rds")

# Meta-analysis ----------------------------------------------------------------

# Define canonical trial and term order
trial_levels <- sort(unique(res_compl_terms$nct_id))
term_levels  <- sort(unique(res_compl_terms$term))

# Ensure coefficients are ordered consistently with block-diagonal VCOV
ma_dat <- res_compl_terms %>%
  mutate(
    ma_id = factor(nct_id, levels = trial_levels),
    ma_param = factor(term, levels = term_levels)
  ) %>%
  arrange(ma_id, ma_param)

# Ensure VCOV matrices are ordered consistently with coefficient table
ma_vcov <- vcov_construct %>%
  arrange(match(nct_id, trial_levels)) %>%
  mutate(
    vcov_matrix = map(
      vcov_matrix, 
      ~ .x[term_levels, term_levels, drop = FALSE]
    )
  )

# Check alignment of estimates and vcov matrices
# No warnings = aligned
stopifnot(
  identical(as.character(ma_dat$ma_id), rep(trial_levels, each = length(term_levels))),
  identical(as.character(ma_dat$ma_param), rep(term_levels, times = length(trial_levels))),
  all(map_lgl(vcov_construct$vcov_matrix, ~ identical(rownames(.x), term_levels))),
  all(map_lgl(vcov_construct$vcov_matrix, ~ identical(colnames(.x), term_levels)))
)

# Construct block-diagonal VCOV
ma_vcov_bdiag <- as.matrix(bdiag(ma_vcov$vcov_matrix))

# Fit meta-analysis model
# Pass trial-specific, participant-level effect estimates on log-HR scale
# Pass block-diagonal VCOV to account for within-trial uncertainty
# Terms are fixed effects
# Random participant-level effects varying across trials
# Define between-trial covariance structure using struct = "DIAG"
# As prepared above, missing terms retained with estimate = 0 and very large
# SEs, contributes negligible information while preserving a common beta and
# VCOV structure across trials
# ma_fit <- rma.mv(
#   yi = estimate,
#   V = ma_vcov_bdiag,
#   mods = ~ ma_param - 1,
#   random = ~ ma_param | ma_id,
#   struct = "DIAG",
#   data = ma_dat,
#   method = "REML",
#   verbose = TRUE
# )

# Save to skip refitting unstructured model
ma_fit <- readRDS("outputs/conduct/ma/out/fit.rds")
# saveRDS(ma_fit, "outputs/conduct/ma/out/fit.rds")

# Extract term-specific heterogeneity
ma_het <- tibble(
  ma_param = levels(ma_dat$ma_param),
  tau2 = ma_fit$tau2,
  tau = sqrt(ma_fit$tau2),
  hr_lower_approx = exp(-sqrt(ma_fit$tau2)),
  hr_upper_approx = exp(sqrt(ma_fit$tau2))
) %>%
  reframe(
    mean_tau = mean(tau, na.rm = TRUE),
    sd_tau = sd(tau, na.rm = TRUE),
    median_tau = median(tau, na.rm = TRUE),
    min_tau = min(tau, na.rm = TRUE),
    max_tau = max(tau, na.rm = TRUE),
    mean_hr_upper_approx = exp(mean_tau),
    median_hr_upper_approx = exp(median_tau)
  ) %>%
  pivot_longer(everything())

# Generate outputs -------------------------------------------------------------

# Tidy meta-analysis results
# - Get hazard ratio estimate summaries
# - Round estimates
# - Join tidy term names
out_tidy <- tidy(ma_fit) %>% 
  mutate(term = gsub("ma_param", "", term)) %>% 
  mutate(
    lo = estimate - 1.96 * std.error,
    hi = estimate + 1.96 * std.error,
    hr = exp(estimate),
    lo_hr = exp(lo),
    hi_hr = exp(hi)
  ) %>% 
  select(term:estimate, lo, hi, hr, lo_hr, hi_hr, everything()) %>% 
  mutate(
    across(estimate:p.value, ~ round(., 3)),
    term_tidy = case_match(
      term,
      "age10" ~ "Age (decades)",
      "bmi" ~ "BMI (kg/m2)",
      "chol_hdl" ~ "HDL cholesterol (mmol/L)",
      "chol_total" ~ "Total cholesterol (mmol/L)",
      "comorbs" ~ "Number of comorbidities",
      "durn" ~ "Duration of type 2 diabetes (years)",
      "egfr" ~ "eGFR (ml/min/1.73m2)",
      "hba1c" ~ "HbA1c (%)",
      "conmeds" ~ "Number of concomitant medications",
      "raceasian" ~ "Asian",
      "raceblack_aa" ~ "Black/African American",
      "raceother_na_unk" ~ "Other race/ethnicity",
      "sbp" ~ "Systolic blood pressure (mmHg)",
      "sex" ~ "Male",
      "smoke" ~ "Current smoker",
      "trttype" ~ "Treatment type: Active comparator",
      .default = term
    ),
    term_tidy = factor(
      term_tidy,
      levels = c(
        "Treatment type: Active comparator", 
        "Age (decades)", "Male", "Asian", "Black/African American",
        "Other race/ethnicity", "Current smoker", 
        "Duration of type 2 diabetes (years)", "BMI (kg/m2)",
        "Number of comorbidities", "Number of concomitant medications",
        "HbA1c (%)", "eGFR (ml/min/1.73m2)", "Systolic blood pressure (mmHg)",
        "HDL cholesterol (mmol/L)", "Total cholesterol (mmol/L)"
      )
    )
  ) %>% 
  filter(!term == "trttype") %>% 
  select(term_tidy, hr, std.error, lo_hr, hi_hr, p.value, everything())

saveRDS(out_tidy, "outputs/conduct/ma/out/ma_tidy.rds")
write_csv(out_tidy, "outputs/conduct/ma/out/ma_tidy.csv")

# Plot hazard ratio
out_hr_plot <- readRDS("outputs/conduct/ma/out/ma_tidy.rds") %>% 
  ggplot(aes(x = hr, y = fct_rev(term_tidy))) +
  geom_point() + 
  geom_linerange(aes(xmin = lo_hr, xmax = hi_hr)) +
  scale_x_continuous(limits = c(0.4, 1.6), n.breaks = 7) +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  theme_classic() +
  labs(x = "Pooled mean hazard ratio (95% confidence intervals)", y = NULL)

ggsave(
  "outputs/conduct/ma/out/ma_plot_hr.png",
  out_hr_plot,
  height = 6,
  width = 6,
  units = "in"
)
