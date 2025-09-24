
source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)

# Import data (variable names, model)
a_imp_df <- colnames(readRDS("mars/design_model/latest_stan_design.rds")$X)
a_imp_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")

# Get posterior draws
b_pos <- as_draws_df(a_imp_mdl) %>%
  as_tibble() %>% 
  pivot_longer(everything())

# Get relevant parameters for attrisk
c_params <- b_pos %>% 
  mutate(
    param = str_extract(name, "^[^\\[]+"),
    draw = str_extract(name, "\\[\\d+"),
    draw = gsub("\\[", "", draw)
  ) %>% 
  filter(param %in% c("mu1", "beta", "a")) %>% 
  select(-name)

# Join covariate names for betas
d_add_names <- c_params %>% 
  left_join(
    tibble(cov = a_imp_df) %>% 
      mutate(
        param = "beta",
        draw = as.character(row_number())
      )
  )

# Summarise draws for mu1 and shape
e_sum <- d_add_names %>% 
  filter(!param %in% c("beta")) %>% 
  group_by(param) %>% 
  reframe(
    mean_param = mean(value),
    median_param = median(value),
    q2.5_param = quantile(value, 0.025),
    q97.5_param = quantile(value, 0.975)
  )

# Summarise draws for betas
e_sum_beta <- d_add_names %>% 
  filter(param %in% c("beta")) %>% 
  group_by(cov) %>% 
  reframe(
    mean_param = mean(value),
    median_param = median(value),
    q2.5_param = quantile(value, 0.025),
    q97.5_param = quantile(value, 0.975)
  )

# Join summaries
f_join <- e_sum %>% full_join(e_sum_beta %>% rename(param = cov))

# Save
saveRDS(f_join, "outputs/design_model/sum_draws_shiny.rds")

# Sample posterior draws from multivariate normal distribution -----------------

# Get draws matrix of target parameters
g_draws_mat <- as_draws_matrix(as.array(a_imp_mdl))
g_keep <- grepl("^mu1|^beta|^a", colnames(g_draws_mat))
g_params <- g_draws_mat[, g_keep, drop = FALSE]

# Separate out betas
g_gomp <- g_params[, !grepl("beta", colnames(g_params)), drop = FALSE]
g_betas <- g_params[, grepl("beta", colnames(g_params)), drop = FALSE]

# Sample Gompertz parameters
h_mu <- colMeans(g_gomp)
h_Sigma <- cov(g_gomp)
h_mvr <- MASS::mvrnorm(1000, mu = h_mu, Sigma = h_Sigma)

# Sample betas
h_mu_beta <- colMeans(g_betas)
h_Sigma_beta <- cov(g_betas)
h_mvr_beta <- MASS::mvrnorm(1000, mu = h_mu_beta, Sigma = h_Sigma_beta)

# Tidy samples
i_tidy <- as_tibble(h_mvr) %>% 
  pivot_longer(everything()) %>% 
  separate(name, into = c("param", "trial"), sep = "\\[") %>% 
  mutate(trial = gsub("\\]", "", trial)) %>% 
  arrange(trial) %>% 
  group_by(trial, param) %>% 
  mutate(sample = row_number())

# Tidy beta samples
i_tidy_beta <- as_tibble(h_mvr_beta) %>% 
  pivot_longer(everything()) %>% 
  separate(name, into = c("beta", "cov"), sep = "\\[") %>% 
  mutate(cov = gsub("\\]", "", cov)) %>% 
  arrange(cov) %>% 
  group_by(cov, beta) %>% 
  mutate(sample = row_number())

# Save
saveRDS(i_tidy, "scratch_data/mvr_draws_post_testing.rds")
saveRDS(i_tidy_beta, "scratch_data/mvr_draws_beta.rds")
