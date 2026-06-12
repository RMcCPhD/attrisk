
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import characteristics
a_imp_bl <- readRDS("processed_data/agg_bl_pooled_arms.rds")

# Heatmap of missingness before imputation
plot_hmap_preimp <- a_imp_bl %>% 
  mutate(across(trttype:one_pom, as.character)) %>% 
  pivot_longer(age_mu:former_never) %>% 
  mutate(na = if_else(is.na(value), 0L, 1L)) %>% 
  ggplot(aes(x = trial_id, y = fct_rev(name), fill = na)) +
  geom_tile() +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank()
  )

ggsave(
  "outputs/conduct_model/agg_bl/hmap_pre_imputation.png",
  plot_hmap_preimp,
  width = 10,
  height = 10,
  units = "in"
)

# Plot before imputation
plot_dens <- a_imp_bl %>% 
  pivot_longer(age_mu:former_never) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.01) +
  facet_wrap(~name, scales = "free")

plot_dens

ggsave(
  "outputs/conduct_model/agg_bl/dist_before_imp.png",
  plot_dens,
  height = 10,
  width = 10,
  units = "in"
)

# Helper functions for setting up proportion imputations
softmax <- function(x) {
  ex <- exp(x - max(x, na.rm = TRUE))
  ex / sum(ex)
}

renorm_simplex <- function(p, eps = 1e-3) {
  p <- pmax(eps, pmin(1 - eps, p))
  s <- sum(p)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(p)))
  p / s
}

logit <- function(p) log(p/(1-p))
inv_logit <- function(z) 1/(1+exp(-z))
alr  <- function(p) log(p[-length(p)] / p[length(p)])
ialr <- function(z) softmax(c(z, 0))

# Format trttype, remove refs for categorical variables
# Remove old variables
# Split out blinding
b_format <- a_imp_bl %>% 
  mutate(
    hba1c_d = case_when(
      trial_id == "NCT00294723" & hba1c_d == 12 ~ 1.2,
      TRUE ~ hba1c_d
    ),
    trttype = factor(trttype),
    blind_ol = if_else(blinding == "open_label", 1L, 0L),
    blind_dbl = if_else(blinding == "double", 1L, 0L),
    blind_trip = if_else(blinding == "triple", 1L, 0L),
    blind_quad = if_else(blinding == "quadruple", 1L, 0L),
  ) %>% 
  select(trial_id:phase, blind_ol:blind_quad, everything()) %>% 
  select(-blinding, -source)

race_cols <- c("asian","black","other","white")

# Function to build dataset with race and smoking status prepped for imputation
prep_df <- function(df) {
  
  # Race log-transformed
  p_race <- as.matrix(df[, race_cols])
  mask_race <- is.finite(p_race)
  p_race[mask_race & p_race <= 0] <- 1e-3
  p_race[mask_race & p_race >= 1] <- 1 - 1e-3
  
  race_mat <- t(apply(p_race, 1, function(p) {
    if (all(is.finite(p))) alr(p) else rep(NA_real_, length(p) - 1)
  }))
  
  race_df <- as_tibble(
    race_mat, 
    .name_repair = ~ paste0("race_alr", seq_len(ncol(race_mat)))
  )
  
  # Smoke - logit
  smoke_cur <- df$current
  smoke_cur[is.finite(smoke_cur) & smoke_cur <= 0] <- 1e-3
  smoke_cur[is.finite(smoke_cur) & smoke_cur >= 1] <- 1 - 1e-3
  smoke_logit <- if_else(is.finite(smoke_cur), logit(smoke_cur), NA_real_)
  
  # Add transformed variables
  df %>% 
    select(-all_of(c(race_cols, "current", "former_never"))) %>% 
    bind_cols(race_df, tibble(smoke_logit = smoke_logit))
}

# Apply function to df
b_amelia_df <- prep_df(b_format)

# Bounds for male proportion
prop_male <- intersect(c("male"), names(b_amelia_df))

prop_bounds <- if (length(prop_male)) {
  do.call(
    rbind,
    lapply(prop_male, function(v) c(match(v, names(b_amelia_df)), 0, 1))
  )
} else {
  matrix(numeric(0), nrow = 0, ncol = 3)
}

# Bounds for numeric variables
b_amelia_df %>% 
  # mutate(across(ends_with("_d"), ~ log(.))) %>% 
  reframe(across(ends_with(c("_mu", "_d")), ~ unlist(
    list(min = min(., na.rm = TRUE), max = max(., na.rm = TRUE))
  ))) %>% 
  pivot_longer(everything()) %>% 
  arrange(name) %>% 
  print(n = 32)

# Limits pre-log SDs
marg_limits <- rbind(
  c(match("age_mu", names(b_format)), 39.5, 74.9),
  c(match("age_d", names(b_format)), 3.67, 13.4),
  c(match("bmi_mu", names(b_format)), 24.2, 37.4),
  c(match("bmi_d", names(b_format)), 1.83, 10.6),
  c(match("egfr_mu", names(b_format)), 24.1, 154),
  c(match("egfr_d", names(b_format)), 4.4, 70),
  c(match("hba1c_mu", names(b_format)), 6.7, 10.2),
  c(match("hba1c_d", names(b_format)), 0.4, 6.3),
  c(match("hdl_mu", names(b_format)), 0.94, 2.7),
  c(match("hdl_d", names(b_format)), 0.02, 0.7),
  c(match("sbp_mu", names(b_format)), 124, 151),
  c(match("sbp_d", names(b_format)), 6.7, 17.4),
  c(match("t2dur_mu", names(b_format)), 0.72, 20.7),
  c(match("t2dur_d", names(b_format)), 0.87, 18),
  c(match("total_chol_mu", names(b_format)), 4, 9.8),
  c(match("total_chol_d", names(b_format)), 0.07, 1.9)
)

# Limits with log-SDs - leads to extreme values
# marg_limits <- rbind(
#   c(match("age_mu", names(b_format)), 39.5, 74.9),
#   c(match("age_d", names(b_format)), 1.3, 2.6),
#   c(match("bmi_mu", names(b_format)), 24.2, 37.4),
#   c(match("bmi_d", names(b_format)), 0.6, 2.36),
#   c(match("egfr_mu", names(b_format)), 24.1, 154),
#   c(match("egfr_d", names(b_format)), 1.48, 4.25),
#   c(match("hba1c_mu", names(b_format)), 6.7, 10.2),
#   c(match("hba1c_d", names(b_format)), -0.92, 0.59),
#   c(match("hdl_mu", names(b_format)), 0.94, 2.7),
#   c(match("hdl_d", names(b_format)), -3.91, -0.36),
#   c(match("sbp_mu", names(b_format)), 124, 151),
#   c(match("sbp_d", names(b_format)), 1.9, 2.86),
#   c(match("t2dur_mu", names(b_format)), 0.72, 20.7),
#   c(match("t2dur_d", names(b_format)), -0.14, 2.89),
#   c(match("total_chol_mu", names(b_format)), 4, 9.8),
#   c(match("total_chol_d", names(b_format)), -2.66, 0.64)
# )

# Combine bounds for male and numeric covariates
bounds <- rbind(
  prop_bounds, 
  marg_limits[!is.na(marg_limits[,1]), , drop = FALSE]
)

# Helper function to obtain mean/sd ratio
# Log-ratio used to constrain sd proportional to means
# eps <- 1e-6
# get_logcv <- function(df, pairs) {
#   for (p in pairs) {
#     mu <- p[1]
#     sd <- p[2]
#     lmu <- paste0(mu, "_logmu")
#     lcv <- paste0(sd, "_logcv")
#     df[[lmu]] <- if_else(
#       is.finite(df[[mu]]), 
#       log(pmax(df[[mu]], 0) + eps),
#       NA_real_
#     )
#   }
#   df
# }
# 
# mu_sd_pairs <- list(
#   c("age_mu", "age_d"),
#   c("bmi_mu", "bmi_d"),
#   c("egfr_mu", "egfr_d"),
#   c("hba1c_mu", "hba1c_d"),
#   c("hdl_mu", "hdl_d"),
#   c("sbp_mu", "sbp_d"),
#   c("t2dur_mu", "t2dur_d"),
#   c("total_chol_mu", "total_chol_d")
# )
# 
# # Get log ratios
# b_logratio <- get_logcv(b_amelia_df, mu_sd_pairs)
# 
# # Add credible ranges for logratios
# logr_limits <- lapply(mu_sd_pairs, function(p) {
#   mu <- p[1]
#   sd <- p[2]
#   lcv <- paste0(sd, "_logcv")
#   ok <- is.finite(b_format[[mu]]) & is.finite(b_format[[sd]]) &
#     b_format[[mu]] > 0 & b_format[[sd]] > 0
#   if (!any(ok)) return(c(-Inf, Inf))
#   cv <- b_format[[sd]][ok] / b_format[[mu]][ok]
#   lo <- max(log(quantile(cv, 0.01, na.rm = TRUE)), log(0.05))
#   hi <- min(log(quantile(cv, 0.99, na.rm = TRUE)), log(2.0))
#   c(lo, hi)
# })
# 
# # Build into boundaries matrix
# logr_bounds <- do.call(rbind, lapply(seq_along(mu_sd_pairs), function(i) {
#   p <- mu_sd_pairs[[i]]
#   lcv <- paste0(p[2], "_logcv")
#   idx <- match(lcv, names(b_format))
#   if (is.na(idx)) return(NULL)
#   c(idx, cv_limits[[i]][1], cv_limits[[i]][2])
# }))
# 
# bounds <- rbind(
#   prop_bounds,
#   marg_limits[!is.na(marg_limits[,1]), , drop = FALSE],
#   logr_bounds
# )

# Run amelia
set.seed(123)

b_impute <- amelia(
  b_amelia_df,
  m = 20,
  idvars  = "trial_id",
  noms    = "trttype",
  empri   = 0.01 * nrow(b_amelia_df),
  bounds  = bounds,
  # logs = c(
  #   "age_d", "bmi_d", "egfr_d", "hba1c_d", "hdl_d",
  #   "sbp_d", "t2dur_d", "total_chol_d"
  # ),
  incheck = TRUE
)

# Function to inverse transform race and smoking to original format
conv_transf <- function(df) {
  alr_cols <- grep("^race_alr", names(df), value = TRUE)
  stopifnot(length(alr_cols) == length(race_cols) - 1)
  
  Z <- as.matrix(df[, alr_cols, drop = FALSE])
  race_rec <- t(apply(Z, 1, function(z) {
    if (all(is.finite(z))) ialr(z) else rep(NA_real_, length(race_cols))
  }))
  race_rec <- t(apply(race_rec, 1, renorm_simplex, eps = 1e-3))
  colnames(race_rec) <- race_cols
  
  # Smoking back from sm_logit (current determines both)
  cur <- inv_logit(df$smoke_logit)
  cur[!is.finite(cur)] <- NA_real_
  cur <- pmin(1 - 1e-3, pmax(1e-3, cur))
  fn  <- 1 - cur
  
  df %>%
    select(-all_of(c(alr_cols, "smoke_logit"))) %>%
    bind_cols(
      as.data.frame(race_rec),
      tibble::tibble(current = cur, former_never = fn)
    )
}

# Function to reconstruct SDs and drop helpers (i.e. logX)
# reconstruct_sds <- function(df, pairs, eps = 1e-6) {
#   for (p in pairs) {
#     mu <- p[1]
#     sd <- p[2]
#     lmu <- paste0(mu, "_logmu")
#     lcv <- paste0(sd, "_logcv")
#     has_mu <- is.finite(df[[mu]])
#     has_lcv <- is.finite(df[[lcv]])
#     rebuilt <- has_mu & has_lcv
#     if (any(rebuilt)) {
#       df[[sd]][[result]] <- (pax(df[[mu]][rebuilt], 0) + eps) * exp(df[[lcv]][rebuilt])
#     }
#     df[[lmu]] <- NULL
#     df[[lcv]] <- NULL
#   }
#   df
# }

# Apply inverse transformation
b_imp_transf <- lapply(b_impute$imputations, conv_transf)

# Add iteration to each imputed dataset
c_iter <- lapply(seq_along(b_imp_transf), \(j) mutate(b_imp_transf[[j]], imp = j))

# Stack, check race and smoking sum to 1
# Exponentiate log-SDs
c_bind <- bind_rows(c_iter) %>% 
  # mutate(across(ends_with("_d"), ~ exp(.))) %>% 
  select(imp, trial_id:male, asian:former_never, everything())

# Impute con meds and comorbs --------------------------------------------------

# Function to obtain base values for negative binomial parameters
# Setup to match the spread in IPD marginals for conmeds and comorbs
base_nb <- function(minx, medx, maxx, p_lo = 0.05, p_hi = 0.95) {
  stopifnot(minx <= medx, medx <= maxx)
  loss <- function(par) {
    mu   <- exp(par[1]); size <- exp(par[2])
    q_lo <- qnbinom(p_lo, size = size, mu = mu)
    q_md <- qnbinom(0.50, size = size, mu = mu)
    q_hi <- qnbinom(p_hi, size = size, mu = mu)
    (q_lo - minx)^2 + 2*(q_md - medx)^2 + (q_hi - maxx)^2
  }
  opt <- optim(
    par   = c(log(medx + 1e-6), log(1.5)),
    fn    = loss,
    method= "L-BFGS-B",
    lower = c(log(1e-8), log(1e-3)),
    upper = c(log(1e6),  log(1e3))
  )
  list(mu = exp(opt$par[1]), size = exp(opt$par[2]))
}

# Draw 20 sets of parameters with log-normal random effects
draw_nb_params <- function(mu_base, size_base, n_draws = 20, 
                           sigma_log_mu = 0.30, sigma_log_size = 0.15) {
  tibble(imp = seq_len(n_draws)) |>
    mutate(
      com_mu   = exp(rnorm(n_draws, log(mu_base$com$mu),   sigma_log_mu)),
      com_size = exp(rnorm(n_draws, log(mu_base$com$size), sigma_log_size)),
      con_mu   = exp(rnorm(n_draws, log(mu_base$con$mu),   sigma_log_mu)),
      con_size = exp(rnorm(n_draws, log(mu_base$con$size), sigma_log_size))
    )
}

# Get global marginals
d_base_comorbs <- base_nb(minx = 0, medx = 1, maxx = 9)
d_base_conmeds <- base_nb(minx = 0, medx = 2, maxx = 48)
d_mu_base <- list(com = d_base_comorbs, con = d_base_conmeds)

# Build tibble of 20 parameter draws per trial (corresponding to each imputation)
d_comorbs_conmeds <- tibble(trial_id = unique(c_bind$trial_id)) |>
  group_by(trial_id) |>
  group_modify(
    ~ draw_nb_params(mu_base = d_mu_base, size_base = NULL, n_draws = 20)
  ) |>
  ungroup()

# Join with rest of variables
e_join <- c_bind %>% 
  left_join(d_comorbs_conmeds) %>% 
  select(imp:former_never, com_mu:con_size, everything()) %>% 
  rename_with(., ~ gsub("com_", "comorbs_", .)) %>% 
  rename_with(., ~ gsub("con_", "conmeds_", .))

# Inspect
e_plot_stack <- e_join %>% 
  pivot_longer(age_mu:conmeds_size) %>% 
  select(trial_id, name, value) %>% 
  mutate(type = "imp") %>% 
  full_join(
    a_imp_bl %>% 
      pivot_longer(age_mu:former_never) %>% 
      mutate(
        name = case_when(
          grepl("asian|black|other|white", name) ~ paste0("race_", name),
          grepl("curr|form", name) ~ paste0("smoke_", name),
          TRUE ~ name
        )
      ) %>% 
      mutate(type = "pre")
  ) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value, colour = type)) +
  geom_density(alpha = 0.01) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~name, scales = "free") +
  theme_classic()

e_plot_stack

ggsave(
  "outputs/conduct_model/agg_bl/dist_before_clp.png",
  e_plot_stack,
  height = 10,
  width = 10,
  units = "in"
)

# Save
saveRDS(e_join, "processed_data/agg_bl_imputed_stack.rds")

# Ensure positive values
# Truncated normal for almost symmetrical variables
# Lognormal or Gamma for skewed
# Inverse CDF to get positive values
pos_bases <- c("egfr", "hdl", "t2dur", "total_chol")
eps <- -1e-8
`%||%` <- function(x, y) if (!is.null(x)) x else y

# Moment matching helper functions
moments_to_lnorm <- function(m, s) {
  m <- pmax(m, eps)
  s <- pmax(s, eps)
  sigma2 <- log(1 + (s^2 / m^2))
  list(meanlog = log(m) - 0.5*sigma2, sdlog = sqrt(sigma2))
}

moments_to_gamma <- function(m, s) {
  m <- pmax(m, eps)
  s <- pmax(s, eps)
  list(shape = (m^2)/(s^2), scale = (s^2)/m)
}

# Choose positive family based on coef of variation
# cv <= 0.3 - truncated normal at 0
# 0.3 - < 1.5 - log-normal
# >= 1.5 - gamma
choose_pos_family <- function(m, s) {
  if (!is.finite(m) || !is.finite(s) || m <= 0 || s <= 0) {
    m <- max(m, 1e-4); s <- max(s, 1e-4)
  }
  cv <- s / m
  if (cv <= 0.30) {
    list(family = "tnorm", pars = list(mean = m, sd = s, a = 0, b = Inf))
  } else if (cv >= 1.5) {
    list(family = "gamma", pars = moments_to_gamma(m, s))
  } else {
    list(family = "lnorm", pars = moments_to_lnorm(m, s))
  }
}

# Build table of mu/sd pairs and positive parameters for vars
pos_param_long <- c_bind %>%
  select(imp, trial_id, ends_with("_mu"), ends_with("_d")) %>%
  pivot_longer(
    cols = matches("(_mu|_d)$"),
    names_to   = c("base", "stat"),
    names_pattern = "(.*)_(mu|d)$",
    values_to  = "val"
  ) %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  filter(base %in% pos_bases) %>%
  mutate(
    mu = pmax(mu, eps),
    d  = pmax(d,  eps)
  ) %>%
  rowwise() %>%
  mutate(
    choice = list(choose_pos_family(mu, d)),
    family = choice$family,
    pars   = list(choice$pars)
  ) %>%
  ungroup() %>%
  mutate(
    meanlog = map_dbl(pars, ~ .x$meanlog %||% NA_real_),
    sdlog   = map_dbl(pars, ~ .x$sdlog   %||% NA_real_),
    shape   = map_dbl(pars, ~ .x$shape   %||% NA_real_),
    scale   = map_dbl(pars, ~ .x$scale   %||% NA_real_),
    tmean   = map_dbl(pars, ~ .x$mean    %||% NA_real_),
    tsd     = map_dbl(pars, ~ .x$sd      %||% NA_real_),
    ta      = map_dbl(pars, ~ .x$a       %||% NA_real_),
    tb      = map_dbl(pars, ~ .x$b       %||% NA_real_)
  ) %>%
  select(
    imp, trial_id, base, mu, d, family,
    meanlog, sdlog, shape, scale, tmean, tsd, ta, tb
  )

# Save
saveRDS(pos_param_long, "processed_data/agg_bl_pos_sup.rds")

# Inverse CDF
inv_cdf_positive <- function(u, family, meanlog, sdlog, shape, scale, tmean, tsd, a = 0, b = Inf) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  if (family == "lnorm") return(qlnorm(u, meanlog = meanlog, sdlog = sdlog))
  if (family == "gamma") return(qgamma(u, shape = shape, scale = scale))
  if (family == "tnorm") return(truncnorm::qtruncnorm(u, a = a, b = b, mean = tmean, sd = tsd))
  stop("Unknown family: ", family)
}

# Check values generated by each variable marginals
c_check_age <- c_bind %>% 
  group_by(trial_id) %>% 
  nest() %>% 
  mutate(
    age = map(data, function(df) {
      rnorm(100, df$age_mu, df$age_d)
    }),
    sbp = map(data, function(df) {
      rnorm(100, df$sbp_mu, df$sbp_d)
    }),
    bmi = map(data, function(df) {
      rnorm(100, df$bmi_mu, df$bmi_d)
    }),
    t2dur = map(data, function(df) {
      rnorm(100, df$t2dur_mu, df$t2dur_d)
    }),
    egfr = map(data, function(df) {
      rnorm(100, df$egfr_mu, df$egfr_d)
    }),
    hba1c = map(data, function(df) {
      rnorm(100, df$hba1c_mu, df$hba1c_d)
    }),
    hdl = map(data, function(df) {
      rnorm(100, df$hdl_mu, df$hdl_d)
    }),
    total_chol = map(data, function(df) {
      rnorm(100, df$total_chol_mu, df$total_chol_d)
    })
  )

c_check_age %>% 
  select(trial_id, age:total_chol) %>% 
  unnest(c(age, sbp, bmi, t2dur, egfr, hba1c, hdl, total_chol)) %>% 
  pivot_longer(age:total_chol) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  scale_x_continuous(n.breaks = 10) +
  facet_wrap(~name, scales = "free") +
  theme_classic()
