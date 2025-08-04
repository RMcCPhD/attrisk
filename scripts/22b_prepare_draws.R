
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
