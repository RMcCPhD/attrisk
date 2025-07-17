
# Interal validation: Evaluate discriminatory ability
# Computing AUC and plotting ROC curve to evaluate
# Using random samples from posterior draws per arm

source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)
library(pROC)

# Import fitted model
a_imp_df <- readRDS("processed_data/latest_dataset.rds")
a_imp_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")
a_imp_draws <- readRDS("mars/design_model/outputs/design_draws.rds")

# Get draws df and save for future use (imported above)
# b_draws_df <- as_draws_df(a_imp_mdl)
# saveRDS(b_draws_df, "mars/design_model/outputs/design_draws.rds")

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
  select(trial_id, arm_n, attr, attr_prop) %>% 
  mutate(row = as.integer(row_number()))

# Join observed variables
b_join <- b_id %>%
  left_join(b_probs) %>% 
  rename(
    obs_attr = attr_prop,
    pred_attr = value
  ) %>% 
  select(-var)

# Randomly select 1000 samples per iteration
set.seed(123)
b_sample_df <- b_join %>% group_by(iter) %>% nest()
b_samples <- sample(b_sample_df$iter, 1000)
b_sampled <- b_sample_df %>% filter(iter %in% b_samples)

# Simulate binary outcomes and compute AUC per arm
sim_attr_auc <- function(df) {
  sim_df <- df %>% mutate(response = rbinom(n(), size = 1, prob = obs_attr))
  auc(
    roc(
      response = sim_df$response,
      predictor = sim_df$pred_attr, 
      quiet = TRUE
    )
  )
}

# Compute auc
c_auc <- b_sampled %>% mutate(auc = map(data, sim_attr_auc))

# Unnest, get object and values
d_auc_dat <- c_auc %>% 
  mutate(
    auc_obj = map(auc, ~ attr(.x, "roc")),
    auc_val = map(auc, ~ .x[[1]])
  )

# Summarise values
d_auc_sum <- d_auc_dat %>% 
  select(auc_val) %>% 
  unnest(auc_val) %>% 
  ungroup() %>% 
  reframe(
    mean_auc = mean(auc_val),
    sd_auc = sd(auc_val),
    med_auc = median(auc_val),
    min_auc = min(auc_val),
    max_auc = max(auc_val)
  )

write_csv(d_auc_sum, "mars/design_model/outputs/auc_sum.csv")

# Plot AUC values across posterior draws
d_auc_dist <- d_auc_dat %>% 
  select(auc_val) %>% 
  unnest(auc_val) %>% 
  ungroup() %>% 
  ggplot(aes(x = auc_val)) +
  geom_density() +
  theme_bw() +
  labs(x = "AUC", y = "Density")

d_auc_dist

ggsave(
  "mars/design_model/outputs/auc_dist.png",
  d_auc_dist,
  width = 8,
  height = 8,
  units = "in"
)

# Get plot data
e_roc <- d_auc_dat %>% 
  mutate(
    roc_dat = map(
      auc_obj,
      ~ data.frame(
        sens = .x$sensitivities,
        spec = .x$specificities
      )
    )
  )

# Get median and credible intervals for plotting
e_roc_sum <- e_roc %>% 
  unnest(roc_dat) %>% 
  mutate(
    fpr = 1 - spec,
    fpr_bin = cut(fpr, breaks = seq(0, 1, by = 0.01))
  ) %>% 
  group_by(fpr_bin) %>% 
  reframe(
    fpr = mean(as.numeric(gsub("\\((.+),.*", "\\1", fpr_bin))),
    median_tpr = median(sens),
    q2.5 = quantile(sens, 0.025),
    q97.5 = quantile(sens, 0.975),
  )

# Plot median and credible intervals
e_roc_plot <- e_roc_sum %>% 
  ggplot(aes(x = fpr, y = median_tpr, ymin = q2.5, ymax = q97.5)) +
  geom_line() +
  geom_ribbon(alpha = 0.3, fill = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_classic() +
  labs(
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)"
  )

e_roc_plot

ggsave(
  "mars/design_model/outputs/roc_plot.png",
  e_roc_plot,
  width = 8,
  height = 8,
  units = "in"
)


