
source("scripts/00_config.R")
library(tidyverse)
library(metafor)
library(broom)

# Import prepared data
a_imp_res <- readRDS("processed_data/conduct_cox_res.rds")
a_imp_vcov <- readRDS("processed_data/conduct_cox_vcov.rds")

# Fit meta-analytic model
# Random effects for term and trial
ma_id <- a_imp_res$nct_id
ma_params <- a_imp_res$term
ma_ests <- a_imp_res$estimate
ma_vcov <- as.matrix(bdiag(a_imp_vcov$cov_matrix))

b_ma <- rma.mv(
  yi = ma_ests,
  V = ma_vcov,
  mods = ~ factor(ma_params) - 1,
  random = ~ ma_params | ma_id
)

# Tidy results
b_ma_tidy <- tidy(b_ma) %>% 
  mutate(term = gsub("factor\\(ma_params\\)", "", term)) %>% 
  mutate(
    lo = estimate - 1.96 * std.error,
    hi = estimate + 1.96 * std.error,
    hr = exp(estimate),
    lo_hr = exp(lo),
    hi_hr = exp(hi)
  ) %>% 
  select(term:estimate, lo, hi, hr, lo_hr, hi_hr, everything()) %>% 
  mutate(across(estimate:p.value, ~ round(., 3)))

saveRDS(b_ma_tidy, "outputs/conduct_model/cox_ma/ma_res.rds")

# Plot pooled log-hazard
plot_ma <- b_ma_tidy %>% 
  arrange(desc(estimate)) %>% 
  ggplot(aes(x = term, y = estimate, ymin = lo, ymax = hi)) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 15) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red") +
  theme_classic() +
  labs(x = NULL, y = "Log estimates")

plot_ma

ggsave(
  "outputs/conduct_model/cox_ma/ma_plot.png",
  plot_ma,
  height = 6,
  width = 8,
  units = "in"
)

# Plot hazard ratio
plot_ma_hr <- b_ma_tidy %>% 
  arrange(desc(estimate)) %>% 
  ggplot(aes(x = term, y = hr, ymin = lo_hr, ymax = hi_hr)) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 15) +
  coord_flip() +
  geom_hline(yintercept = 1, colour = "red") +
  theme_classic() +
  labs(x = NULL, y = "Log estimates")

plot_ma_hr

ggsave(
  "outputs/conduct_model/cox_ma/ma_plot_hr.png",
  plot_ma_hr,
  height = 6,
  width = 8,
  units = "in"
)





# Fit model (doesn't seem to take in variance-covariance to perform
# meta-analysis accounting for correlation between variables)

# Take vcovs as list
b_vcov <- a_imp_vcov$cov_matrix
names(b_vcov) <- unique(a_imp_res$nct_id)

b_fit <- brm(
  formula = estimate ~ term + (term|gr(nct_id, cov = unlist(cov_matrix))),
  data = a_imp_res %>% left_join(a_imp_vcov),
  family = gaussian(link = "identity"),
  chains = 4,
  iter = 3000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "rstan",
  control = list(adapt_delta = 0.95)
)

