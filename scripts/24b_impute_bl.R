
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import bl chars
a_imp_bl <- readRDS("processed_data/bl_chars_noimp.rds")

# Check missingness
a_imp_bl %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  print(n = 25)

# Setup for amelia -------------------------------------------------------------

# Format trttype, remove refs for proportions
b_format <- a_imp_bl %>% 
  mutate(trttype = factor(trttype)) %>% 
  select(-c(white, former_never))

# Set bounds for proportion imputations
prop_vars <- intersect(
  c("female_prop", "asian", "black", "other", "current"), 
  names(b_format)
)

bounds <- do.call(
  rbind, 
  lapply(prop_vars, function(v) c(match(v, names(b_format)), 0, 1))
)

# Set bounds for mean and sd imputations
b_format %>% 
  reframe(across(contains(c("mean", "sd")), ~ unlist(
    list(min = min(., na.rm = TRUE), max = max(., na.rm = TRUE))
  ))) %>% 
  pivot_longer(everything()) %>% 
  print(n = 32)

marg_limits <- rbind(
  c(match("mean_age", names(b_format)), 39.5, 74.9),
  c(match("sd_age", names(b_format)), 3.8, 12.8),
  c(match("mean_sbp", names(b_format)), 124, 151),
  c(match("sd_sbp", names(b_format)), 6.7, 17.4),
  c(match("mean_bmi", names(b_format)), 24.2, 37.4),
  c(match("sd_bmi", names(b_format)), 2.6, 10.6),
  c(match("mean_t2dur", names(b_format)), 0.7, 20.7),
  c(match("sd_t2dur", names(b_format)), 1.4, 10.1),
  c(match("mean_egfr", names(b_format)), 24.1, 154),
  c(match("sd_egfr", names(b_format)), 0.25, 70),
  c(match("mean_hba1c", names(b_format)), 6.7, 10.2),
  c(match("sd_hba1c", names(b_format)), 0.35, 1.8),
  c(match("mean_hdl", names(b_format)), 0.94, 2.7),
  c(match("sd_hdl", names(b_format)), 0.02, 0.07),
  c(match("mean_total_chol", names(b_format)), 4, 9.8),
  c(match("sd_total_chol", names(b_format)), 0.07, 1.9)
)

bounds <- rbind(
  bounds, 
  marg_limits[!is.na(marg_limits[,1]), , drop = FALSE]
)

# Run amelia -------------------------------------------------------------------

set.seed(123)

b_impute <- amelia(
  b_format,
  m = 20,
  idvars  = "trial_id",
  noms    = "trttype",
  empri   = 0.01 * nrow(b_format),
  bounds  = bounds,
  incheck = TRUE
)

# Fix race and smoking groups imputed, sum to 1
fix_grps <- function(df) {
  
  # Smoking groups
  df$current <- pmin(pmax(df$current, 0), 1)
  df$former_never <- 1 - df$current
  
  # Race groups
  race_cols <- c("asian", "black", "other")
  race_X <- as.matrix(df[race_cols])
  race_X[!is.finite(race_X)] <- NA
  race_X[is.na(race_X)] <- 0
  race_X[race_X < 0] <- 0
  race_X[race_X > 1] <- 1
  
  race_s <- rowSums(race_X)
  zero <- !is.finite(race_s) | race_s <= 0
  if (any(zero)) race_X[zero, ] <- 1/3
  
  race_s <- rowSums(race_X)
  race_big <- race_s > 1
  if (any(race_big)) race_X[race_big, ] <- race_X[race_big, ] / race_s[race_big]
  
  df[race_cols] <- race_X
  df$white <- 1 - rowSums(race_X)
  df$white <- pmin(pmax(df$white, 0), 1)
  
  df %>%
    select(trial_id:other, white, everything()) %>% 
    mutate(
      across(mean_age:former_never, ~ round(., 2)),
      race_total = asian + black + other + white,
      white = case_when(
        race_total == 0.99 ~ white + 0.01,
        race_total == 1.01 ~ white - 0.01,
        TRUE ~ white
      ),
      race_total = asian + black + other + white,
      smoke_total = current + former_never
    )
}

b_fixed <- lapply(b_impute$imputations, fix_grps)

# Add imputation iteration to datasets, bind rows
c_bind <- bind_rows(
  lapply(seq_along(b_fixed), \(j) mutate(b_fixed[[j]], imp = j))
)

# Get marginals across imputations ---------------------------------------------

# Check distributions of variable marginals
c_bind %>% 
  pivot_longer(mean_age:former_never) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value, group = imp)) +
  geom_density(alpha = 0.01) +
  facet_wrap(~name, scales = "free")

# Pool variables across imputations based on distribution
# Normally distributed
pool_normal <- function(df, mu_var, delta_var) {
  mu <- do.call(cbind, lapply(df, `[[`, mu_var))
  delta <- do.call(cbind, lapply(df, `[[`, delta_var))
  mu_rows <- rowMeans(mu, na.rm = TRUE)
  var <- rowMeans(delta^2, na.rm = TRUE) + matrixStats::rowVars(mu)
  tibble(mu = mu_rows, sd = sqrt(pmax(var, 0)))
}

pooled_age <- pool_normal(b_fixed, "mean_age", "sd_age")
pooled_egfr <- pool_normal(b_fixed, "mean_egfr", "sd_egfr")
pooled_hba1c <- pool_normal(b_fixed, "mean_hba1c", "sd_hba1c")

# Gamma distributed
pool_gamma <- function(df, )

# Inspect distributions

# ----- NORMAL-LIKE (e.g., age, sbp) ------------------------------------------
pool_normal <- function(completed, mean_col, sd_col) {
  M <- do.call(cbind, lapply(completed, `[[`, mean_col))  # I x m
  S <- do.call(cbind, lapply(completed, `[[`, sd_col))    # I x m
  mu  <- rowMeans(M, na.rm = TRUE)
  sd2 <- rowMeans(S^2, na.rm = TRUE) + rowVars(M, na.rm = TRUE)
  tibble(mu = mu, sd = sqrt(pmax(sd2, 0)))
}

# ----- LOGNORMAL-LIKE (e.g., duration; positive skew) ------------------------
pool_lognormal <- function(completed, mean_col, sd_col) {
  M <- do.call(cbind, lapply(completed, `[[`, mean_col))
  S <- do.call(cbind, lapply(completed, `[[`, sd_col))
  # per-imputation log-scale params
  sigma2L <- log1p((S^2) / (M^2))
  muL     <- log(M) - 0.5 * sigma2L
  # mixture moments on original scale
  m1  <- rowMeans(exp(muL + 0.5 * sigma2L), na.rm = TRUE)
  v_w <- rowMeans(exp(2*muL + sigma2L) * (exp(sigma2L) - 1), na.rm = TRUE)  # within
  v_b <- rowVars(exp(muL + 0.5 * sigma2L), na.rm = TRUE)                     # between
  v   <- v_w + v_b
  # match a single lognormal
  sigma2L_star <- log1p(v / (m1^2))
  muL_star     <- log(m1) - 0.5 * sigma2L_star
  tibble(muL = muL_star, sdL = sqrt(pmax(sigma2L_star, 0)))
}

# ----- SINGLE PROPORTION as BETA (e.g., female_prop) -------------------------
# Combines between-imp variance with binomial within-imp variance via n_eff
pool_beta <- function(completed, prop_col, n_eff) {
  P <- do.call(cbind, lapply(completed, `[[`, prop_col))  # I x m
  pbar <- rowMeans(P, na.rm = TRUE)
  v_b  <- rowVars(P, na.rm = TRUE)
  v_w  <- pbar * (1 - pbar) / pmax(n_eff, 1)              # binomial approx
  v    <- v_b + v_w
  # method-of-moments Beta: Var = p(1-p)/(alpha0+1)  => alpha0 = p(1-p)/Var - 1
  alpha0 <- pmax(pbar*(1 - pbar) / pmax(v, 1e-12) - 1, 2) # clamp >=2 for stability
  alpha <- pmax(pbar * alpha0, 1e-6)
  beta  <- pmax((1 - pbar) * alpha0, 1e-6)
  tibble(alpha = alpha, beta = beta)
}

# ----- COMPOSITION as DIRICHLET (e.g., race/smoking) -------------------------
# cols: vector of column names (must already sum ~1 in each completed dataset)
pool_dirichlet <- function(completed, cols, n_eff) {
  mats <- lapply(completed, function(d) as.matrix(d[cols]))   # list of I x K
  Pbar <- Reduce(`+`, mats) / length(mats)                    # I x K
  # between-imputation variance per part
  Vb   <- Reduce(`+`, lapply(mats, function(M) (M - Pbar)^2)) / (length(mats) - 1)
  # within-imputation (multinomial) variance approx: Var(p_i) ~ p_i(1-p_i)/(n_eff + 1)
  Vw   <- Pbar * (1 - Pbar) / (pmax(n_eff, 1) + 1)
  Vtot <- Vb + Vw                                              # I x K
  # Solve alpha0 from Var_dirichlet(p_i) ≈ p_i(1-p_i)/(alpha0+1)
  alpha0_hat <- rowMeans( Pbar * (1 - Pbar) / pmax(Vtot, 1e-12) - 1, na.rm = TRUE )
  alpha0     <- pmax(alpha0_hat, 2)                            # clamp for stability
  Alpha      <- Pbar * alpha0                                  # I x K
  as_tibble(Alpha, .name_repair = ~cols) %>% mutate(alpha0 = alpha0)
}

# --------- BUILD ONE CONSENSUS TABLE PER ARM ---------------------------------
# completed: list of m completed, renormalised datasets
# keys: columns that define an arm (e.g., c("trial_id","trttype"))
build_consensus <- function(completed, keys = c("trial_id","trttype")) {
  base <- completed[[1]][keys]
  n_eff <- if ("n" %in% names(completed[[1]])) completed[[1]]$n else rep(200, nrow(base))  # fallback
  
  age_ps   <- pool_normal(completed,   "mean_age",   "sd_age")   %>% rename(age_mu = mu, age_sd = sd)
  sbp_ps   <- pool_normal(completed,   "mean_sbp",   "sd_sbp")   %>% rename(sbp_mu = mu, sbp_sd = sd)
  dur_ps   <- pool_lognormal(completed,"mean_t2dur", "sd_t2dur") %>% rename(t2dur_muL = muL, t2dur_sdL = sdL)
  bmi_ps   <- pool_lognormal(completed,"mean_bmi",   "sd_bmi")   %>% rename(bmi_muL = muL,  bmi_sdL  = sdL)
  
  fem_ps   <- pool_beta(completed, "female_prop", n_eff) %>% rename(female_alpha = alpha, female_beta = beta)
  
  race_ps  <- pool_dirichlet(completed, c("asian","black","other","white"), n_eff) %>%
    rename(alpha_asian = asian, alpha_black = black, alpha_other = other, alpha_white = white)
  
  smoke_ps <- pool_dirichlet(completed, c("current","former_never"), n_eff) %>%
    rename(alpha_current = current, alpha_former_never = former_never)
  
  dplyr::bind_cols(base, age_ps, sbp_ps, dur_ps, bmi_ps, fem_ps, race_ps, smoke_ps)
}





