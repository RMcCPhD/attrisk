
# Interal validation
# Using random samples from posterior draws per arm

source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)
library(loo)

# Import fitted model
a_imp_df <- readRDS("processed_data/latest_dataset.rds")
a_imp_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")
a_imp_draws <- readRDS("mars/design_model/outputs/design_draws.rds")

# Get draws df and save for future use (imported above)
# b_draws_df <- as_draws_df(a_imp_mdl)
# saveRDS(b_draws_df, "mars/design_model/outputs/design_draws.rds")

# Brier score ------------------------------------------------------------------

# Get predicted attrition
b_probs <- a_imp_draws %>% 
  pivot_longer(everything()) %>% 
  filter(grepl("^p", name)) %>% 
  separate(name, into = c("var", "row"), sep = "\\[") %>% 
  mutate(
    row = gsub("\\]", "", row),
    row = as.integer(row),
    iter = rep(1:(n() / 366), each = 366)
  )

b_id <- a_imp_df %>% 
  select(trial_id, trttype, arm_n, attr, attr_prop) %>% 
  mutate(row = as.integer(row_number()))

# Join observed variables
b_join <- b_id %>%
  left_join(b_probs) %>% 
  rename(
    obs_attr = attr_prop,
    pred_attr = value
  ) %>% 
  select(-var) %>% 
  mutate(
    trttype = case_when(
      trttype == 0L ~ "Plc",
      TRUE ~ "AC"
    ),
    trttype = factor(trttype, levels = c("Plc", "AC"))
  )

# Compute Brier score
b_brier <- b_join %>% 
  group_by(trttype, iter) %>% 
  reframe(brier = mean((obs_attr - pred_attr)^2))

# Summarise Brier scores across posterior draws
b_brier_sum <- b_brier %>% 
  group_by(trttype) %>% 
  reframe(
    mean_brier = mean(brier),
    sd_brier = sd(brier),
    med_brier = median(brier),
    min_brier = min(brier),
    max_brier = max(brier)
  )

write_csv(b_brier_sum, "mars/design_model/outputs/brier_sum.csv")

# Plot Brief scores
b_brier_dist <- b_brier %>%
  mutate(trttype = as.factor(trttype)) %>% 
  ggplot(aes(x = brier, group = trttype)) +
  geom_density(aes(fill = trttype), alpha = 0.3, colour = "black") +
  scale_x_continuous(n.breaks = 8) +
  theme_classic() +
  labs(x = "Brier score", y = "Density", fill = "Treatment type")

b_brier_dist

ggsave(
  "mars/design_model/outputs/brier_dist.png",
  b_brier_dist,
  width = 8,
  height = 8,
  units = "in"
)

# Calibration plot -------------------------------------------------------------

# Summarise predicted attrition
# Create bins of predicted attrition
c_calib_df <- b_join %>% 
  group_by(row, trttype) %>% 
  reframe(
    obs_attr = first(obs_attr),
    pred_attr = mean(pred_attr)
  ) %>% 
  mutate(
    bin = ntile(pred_attr, 10)
  ) %>% 
  group_by(bin, trttype) %>% 
  reframe(
    mean_pred = mean(pred_attr),
    mean_obs = mean(obs_attr),
    n = n()
  )

# Plot
c_calib_plot <- c_calib_df %>% 
  ggplot(aes(x = mean_pred, y = mean_obs, group = trttype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point() +
  geom_line() +
  facet_wrap(~ trttype, scales = "free") +
  theme_classic() +
  labs(x = "Mean predicted attrition", y = "Observed attrition")

c_calib_plot

ggsave(
  "mars/design_model/outputs/calib_plot.png",
  c_calib_plot,
  width = 8,
  height = 8,
  units = "in"
)