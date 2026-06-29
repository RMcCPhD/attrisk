
# Generate Sobol sequences representing synthetic IPD participant sets
# Sequences generated for each AgD trial

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import prepared Sobol materials
imp_dist <- readRDS("processed_data/sobol_dist.rds")
imp_cor <- readRDS("processed_data/sobol_cor_mat.rds")
imp_best <- read_csv("scratch_data/sobol_covariate_dist.csv")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import pooled log-hazard estimates from IPD models
imp_loghr <- readRDS("outputs/conduct/ma/out/ma_tidy.rds") %>%
  rename(loghr = estimate)

# Import covariate scaling
imp_scale <- read_csv("vivli/v2/Summaries/covariate_mu_d.csv")

# Helper: Centre variables on mean and scale by standard deviation
scale_centre_cov <- function(dat, imp_scale) {
  
  scale_lkp <- imp_scale %>% filter(name %in% names(dat))
  mu <- setNames(scale_lkp$mu, scale_lkp$name)
  d <- setNames(scale_lkp$d, scale_lkp$name)
  
  dat %>% 
    mutate(
      across(
        all_of(scale_lkp$name),
        ~ (.x - mu[cur_column()]) / d[cur_column()]
      )
    )
}

# Helper: assign mutually exclusive race indicators.
# Race is generated as one categorical variable and then dummy-coded into
# non-reference race indicators. This avoids impossible profiles where more
# than one race/ethnicity indicator is active.
assign_race_indicators <- function(dat, p_asian, p_black_aa, p_other, shift = 0) {
  
  p_asian <- as.numeric(p_asian)
  p_black_aa <- as.numeric(p_black_aa)
  p_other <- as.numeric(p_other)
  
  p_asian <- max(p_asian, 0, na.rm = TRUE)
  p_black_aa <- max(p_black_aa, 0, na.rm = TRUE)
  p_other <- max(p_other, 0, na.rm = TRUE)
  
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
  
  race <- case_when(
    race_u < cuts[["asian"]] ~ "asian",
    race_u < cuts[["black_aa"]] ~ "black_aa",
    race_u < cuts[["other"]] ~ "other",
    TRUE ~ "white"
  )
  
  dat %>%
    mutate(
      asian = as.integer(race == "asian"),
      black_aa = as.integer(race == "black_aa"),
      other = as.integer(race == "other")
    )
}

# Prepare distribution parameters ----------------------------------------------

# Check counts
imp_best

# Helper functions
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

# Add parameters
prep_dist <- imp_dist %>% 
  mutate(
    across(ends_with("_d"), ~ pmax(.x, .Machine$double.eps)),
    
    # log-normal: bmi, chol_hdl, sbp, tot_hdl_ratio
    bmi_meanlog = get_lnorm_meanlog(bmi_mu, bmi_d),
    bmi_sdlog   = get_lnorm_sdlog(bmi_mu, bmi_d),
    
    chol_hdl_meanlog = get_lnorm_meanlog(chol_hdl_mu, chol_hdl_d),
    chol_hdl_sdlog   = get_lnorm_sdlog(chol_hdl_mu, chol_hdl_d),
    
    sbp_meanlog = get_lnorm_meanlog(sbp_mu, sbp_d),
    sbp_sdlog   = get_lnorm_sdlog(sbp_mu, sbp_d),
    
    # gamma: chol_total
    chol_total_shape = get_gamma_shape(chol_total_mu, chol_total_d),
    chol_total_rate  = get_gamma_rate(chol_total_mu, chol_total_d),
    
    # Weibull: durn, egfr
    durn_shape = get_weibull_shape(durn_mu, durn_d),
    durn_scale = get_weibull_scale(durn_mu, durn_shape),
    
    egfr_shape = get_weibull_shape(egfr_mu, egfr_d),
    egfr_scale = get_weibull_scale(egfr_mu, egfr_shape)
  )

# Generate Sobol sequences -----------------------------------------------------

# Nest AgD trials, generated 64 Sobol sequences per trial
# Join sobol data (distribution summaries, correlation matrices)
sobol_agd <- imp_agd %>% 
  group_by(trial_id) %>% 
  nest(.key = "agd") %>% 
  cbind(prep_dist %>% nest(.key = "dist")) %>% 
  cbind(imp_cor %>% select(cor = out))

# Generate 64 Sobol sequences per AgD trial
# Race indicators are generated separately as one mutually exclusive
# categorical variable, so remove race indicators from the ecological
# correlation matrix used by add_integration().
sobol_gen <- sobol_agd %>% 
  mutate(
    race_shift = ((row_number() - 1) * 0.6180339887498949) %% 1,
    sobol = map2(dist, cor, ~ {
      
      race_cols <- c("asian", "black_aa", "other")
      cor_no_race <- .y
      
      if (!is.null(rownames(cor_no_race)) && all(race_cols %in% rownames(cor_no_race))) {
        keep_cols <- setdiff(rownames(cor_no_race), race_cols)
        cor_no_race <- cor_no_race[keep_cols, keep_cols, drop = FALSE]
      }
      
      add_integration(
        .x,
        age           = distr(qnorm,    mean = .x$age_mu, sd = .x$age_d),
        male_p        = distr(qbern,    prob = 1 - .x$female_p),
        bmi           = distr(qlnorm,   meanlog = .x$bmi_meanlog, sdlog = .x$bmi_sdlog),
        durn          = distr(qweibull, shape = .x$durn_shape, scale = .x$durn_scale),
        sbp           = distr(qlnorm,   meanlog = .x$sbp_meanlog, sdlog = .x$sbp_sdlog),
        chol_hdl      = distr(qlnorm,   meanlog = .x$chol_hdl_meanlog, sdlog = .x$chol_hdl_sdlog),
        chol_total    = distr(qgamma,   shape = .x$chol_total_shape, rate = .x$chol_total_rate),
        hba1c         = distr(qnorm,    mean = .x$hba1c_mu, sd = .x$hba1c_d),
        egfr          = distr(qweibull, shape = .x$egfr_shape, scale = .x$egfr_scale),
        cursmoke      = distr(qbern,    prob = .x$cursmoke),
        comorbs       = distr(qpois,    lambda = .x$comorbs_mu),
        conmeds       = distr(qpois,    lambda = .x$conmeds_mu),
        cor = cor_no_race,
        n_int = 64
      )
    }),
    sobol_unnest = map2(sobol, race_shift, ~ {
      .x %>% 
        select(starts_with(".int")) %>% 
        rename_with(., ~ gsub(".int_", "", .)) %>% 
        unnest(cols = everything())
    })
  )

# Add mutually exclusive race indicators before downstream transformation
sobol_gen <- sobol_gen %>%
  mutate(
    sobol_unnest = pmap(
      list(sobol_unnest, dist, race_shift),
      ~ assign_race_indicators(
        dat = ..1,
        p_asian = ..2$asian,
        p_black_aa = ..2$black_aa,
        p_other = ..2$other,
        shift = ..3
      )
    )
  )

# Save outputs
saveRDS(sobol_gen, "processed_data/sobol_out.rds")

# Transform sobol points - centre on mean, scale by standard deviation
# Summaries for transformation exported from Vivli, across all trials
sobol_unnest <- sobol_gen %>% 
  mutate(
    sobol_rdy = map(sobol_unnest, ~ {
      .x %>% 
        rename(
          age10 = age,
          sex = male_p
        ) %>%
        mutate(age10 = age10 / 10) %>%
        scale_centre_cov(imp_scale)
    })
  )

# Inspect distributions
sobol_plot <- sobol_unnest %>% 
  select(trial_id, sobol_rdy) %>% 
  unnest(sobol_rdy) %>% 
  ungroup() %>% 
  pivot_longer(-trial_id) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  theme_classic()

ggsave(
  "outputs/conduct/sobol/sobol_distributions.png",
  width = 8,
  height = 6,
  units = "in"
)

# Prepare CITL inputs ----------------------------------------------------------

citl_ma <- imp_loghr %>% 
  select(term, loghr) %>% 
  pivot_wider(names_from = "term", values_from = "loghr") %>% 
  rename_with(., ~ paste0(., "_loghr"))

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

# Parameters
I <- nrow(imp_agd)      # Number of trial arms
J <- 64L                # Number of Sobol points
Pnm <- names(citl_beta) # Predictor names
pK <- length(Pnm)       # Number of predictors

# Construct a matrix of sobol sequences per trial (IxJ) x pK
Z_blocks <- vector("list", I)

for (i in seq_len(I)) {
  trial <- imp_agd$trial_id[i]
  idx <- match(trial, sobol_unnest$trial_id)
  M <- sobol_unnest$sobol_rdy[[idx]]
  Z_blocks[[i]] <- as.matrix(M)
}

Z <- do.call(rbind, Z_blocks)

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

# Supplementary Sobol summaries ------------------------------------------------

summ_vec <- function(x, metric) {
  tibble(
    metric = metric,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE)
  )
}

summ_long <- function(dat, id_cols = "trial_id") {
  dat %>%
    select(all_of(id_cols), where(is.numeric)) %>%
    pivot_longer(-all_of(id_cols)) %>%
    group_by(name) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      q025 = quantile(value, 0.025, na.rm = TRUE),
      q975 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

sobol_raw <- sobol_unnest %>%
  select(trial_id, sobol_unnest) %>%
  unnest(sobol_unnest) %>%
  ungroup()

sobol_scaled <- sobol_unnest %>%
  select(trial_id, sobol_rdy) %>%
  unnest(sobol_rdy) %>%
  ungroup()

sobol_raw_summary <- summ_long(sobol_raw)
sobol_scaled_summary <- summ_long(sobol_scaled)

sobol_target_summary <- bind_rows(
  tibble(
    name = c("age", "bmi", "chol_hdl", "chol_total", "hba1c", "durn", "egfr", "sbp"),
    target_mean = c(
      imp_dist$age_mu,
      imp_dist$bmi_mu,
      imp_dist$chol_hdl_mu,
      imp_dist$chol_total_mu,
      imp_dist$hba1c_mu,
      imp_dist$durn_mu,
      imp_dist$egfr_mu,
      imp_dist$sbp_mu
    ),
    target_sd = c(
      imp_dist$age_d,
      imp_dist$bmi_d,
      imp_dist$chol_hdl_d,
      imp_dist$chol_total_d,
      imp_dist$hba1c_d,
      imp_dist$durn_d,
      imp_dist$egfr_d,
      imp_dist$sbp_d
    ),
    target_proportion = NA_real_
  ),
  tibble(
    name = c("male_p", "cursmoke", "asian", "black_aa", "other"),
    target_mean = NA_real_,
    target_sd = NA_real_,
    target_proportion = c(
      1 - imp_dist$female_p,
      imp_dist$cursmoke,
      imp_dist$asian,
      imp_dist$black_aa,
      imp_dist$other
    )
  ),
  tibble(
    name = c("comorbs", "conmeds"),
    target_mean = c(imp_dist$comorbs_mu, imp_dist$conmeds_mu),
    target_sd = c(imp_dist$comorbs_d, imp_dist$conmeds_d),
    target_proportion = NA_real_
  )
)

sobol_generation_check <- sobol_raw_summary %>%
  left_join(sobol_target_summary, by = "name") %>%
  mutate(
    mean_minus_target = case_when(
      !is.na(target_mean) ~ mean - target_mean,
      !is.na(target_proportion) ~ mean - target_proportion,
      TRUE ~ NA_real_
    ),
    sd_minus_target = if_else(!is.na(target_sd), sd - target_sd, NA_real_)
  )

sobol_profile_check_dat <- sobol_raw %>%
  mutate(
    n_race_indicators = rowSums(across(c(asian, black_aa, other)), na.rm = TRUE),
    multiple_race_indicators = n_race_indicators > 1,
    no_nonwhite_indicator = n_race_indicators == 0,
    male = male_p,
    female = 1 - male_p
  )

sobol_profile_checks <- sobol_profile_check_dat %>%
  summarise(
    n_profiles = n(),
    n_trials = n_distinct(trial_id),
    prop_multiple_race_indicators = mean(multiple_race_indicators, na.rm = TRUE),
    prop_no_nonwhite_indicator = mean(no_nonwhite_indicator, na.rm = TRUE),
    prop_asian = mean(asian, na.rm = TRUE),
    prop_black_aa = mean(black_aa, na.rm = TRUE),
    prop_other = mean(other, na.rm = TRUE),
    prop_current_smoker = mean(cursmoke, na.rm = TRUE),
    prop_female_generated = mean(female, na.rm = TRUE),
    prop_male_generated = mean(male, na.rm = TRUE),
    n_missing_values = sum(is.na(as.matrix(
      sobol_profile_check_dat %>% select(where(is.numeric))
    )))
  )

write_csv(sobol_raw_summary, "outputs/conduct/sobol/sobol_raw_summary.csv")
write_csv(sobol_scaled_summary, "outputs/conduct/sobol/sobol_scaled_summary.csv")
write_csv(sobol_generation_check, "outputs/conduct/sobol/sobol_generation_check.csv")
write_csv(sobol_profile_checks, "outputs/conduct/sobol/sobol_profile_checks.csv")

saveRDS(sobol_unnest, "processed_data/sobol_unnest.rds")

# Summarise participant-level linear predictors and hazard multipliers ----------

Z <- Z[, names(citl_beta), drop = FALSE]

sobol_lp <- as.vector(Z %*% citl_beta)
sobol_hr_multiplier <- exp(sobol_lp)

sobol_lp_dat <- tibble(
  arm = rep(seq_len(I), each = J),
  trial_id = rep(imp_agd$trial_id, each = J),
  profile = rep(seq_len(J), times = I),
  lp = sobol_lp,
  hr_multiplier = sobol_hr_multiplier
)

sobol_lp_arm <- sobol_lp_dat %>%
  group_by(trial_id) %>%
  reframe(
    mean_lp = mean(lp, na.rm = TRUE),
    median_lp = median(lp, na.rm = TRUE),
    mean_hr_multiplier = mean(hr_multiplier, na.rm = TRUE),
    median_hr_multiplier = median(hr_multiplier, na.rm = TRUE)
  )

sobol_lp_summary <- bind_rows(
  summ_vec(sobol_lp_dat$lp, "participant_linear_predictor"),
  summ_vec(sobol_lp_dat$hr_multiplier, "participant_hazard_multiplier"),
  summ_vec(sobol_lp_arm$mean_lp, "arm_mean_linear_predictor"),
  summ_vec(sobol_lp_arm$mean_hr_multiplier, "arm_mean_hazard_multiplier")
)

write_csv(sobol_lp_dat, "outputs/conduct/sobol/sobol_linear_predictors.csv")
write_csv(sobol_lp_arm, "outputs/conduct/sobol/sobol_linear_predictors_by_arm.csv")
write_csv(sobol_lp_summary, "outputs/conduct/sobol/sobol_linear_predictor_summary.csv")

