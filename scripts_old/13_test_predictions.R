
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)
library(flexsurv)

# Import datasets
a_imp_raw <- readRDS("processed_data/definitive_dataset.rds")
a_imp_df <- readRDS("mars/design_model_mk4/stan_dat.rds")
a_imp_design <- readRDS("mars/design_model_mk4/outputs/design4_sum_draws.rds")

# Scaling parameters
#' Male prop: Mean=0.5535239, SD=0.09109676
#' Age: 4.044924, 0.07224443
#' T2D dur: 1.923038, 0.6038713
#' Sample: 6.090808, 0.6236182
#' Arm size: 5.109074, 0.6125837
#' Duration: -0.7143654, 0.438717
#' Facilities: 4.068796, 1.021258
scale(a_imp_raw$male_prop)
scale(log(a_imp_raw$age))
scale(log(a_imp_raw$t2dur))
scale(log(a_imp_raw$sample_size))
scale(log(a_imp_raw$arm_n))
scale(log(a_imp_raw$dur_yr))
scale(log(a_imp_raw$n_facilities))

summary(a_imp_raw[c(12:15, 31, 34:36)])
covs <- a_imp_df$X
summary(covs[, c(15:18)])
cov_names <- colnames(a_imp_df$X)

# Prepare draws
posterior_samples <- a_imp_design %>% 
  filter(grepl("^beta\\[|^mu1\\[|^a\\[", variable)) %>% 
  separate(variable, into = c("variable", "draw"), sep = "\\[") %>% 
  mutate(draw = gsub("\\]", "", draw)) %>% 
  select(variable, draw, mean) %>% 
  pivot_wider(names_from = "variable", values_from = "mean")

# User inputted design features
# Scale numeric variables in design features by means and sds of scaling
input <- data.frame(
  trttype = 0L,
  p3 = 1L,
  p4 = 0L,
  po = 0L,
  bdbl = 1L,
  btrip = 0L,
  bquad = 0L,
  bol = 0L,
  reg_am = 1L,
  reg_eur = 1L,
  reg_wpac = 1L,
  reg_afr = 0L,
  reg_sea = 0L,
  reg_em = 0L,
  sample = 500,
  arm_size = 250,
  dur = 26/52,
  facilities = 100,
  arms_2 = 1L,
  arms_3 = 0L,
  arms_abv = 0L,
  pom_1 = 1L,
  pom_2 = 0L,
  pom_abv = 0L
) %>% 
  pivot_longer(
    everything(), names_to = "variable", values_to = "value"
  ) %>% 
  mutate(
    scaling = case_when(
      variable == "sample" ~ as.numeric(scale(log(value), center = 6.090808, scale = 0.6236182)),
      variable == "arm_size" ~ as.numeric(scale(log(value), center = 5.109074, scale = 0.6125837)),
      variable == "dur" ~ as.numeric(scale(log(value), center = -0.7143654, scale = 0.438717)),
      variable == "facilities" ~ as.numeric(scale(log(value), center = 4.068796, scale = 1.021258)),
      TRUE ~ value
    )
  )

# Test with just means ---------------------------------------------------------

# Extract beta
beta <- posterior_samples %>% filter(!is.na(beta))

# Sequence trial duration
dt <- seq(1:24)

# Extract mean mu1 and a
mu1 <- mean(posterior_samples$mu1)
a <- mean(posterior_samples$a)

# Get matrix of betas and compute dot products (X[i] * beta[i])
beta <- as.matrix(beta$beta)
b <- exp(mu1 + t(input$scaling) %*% beta)

# Tibble of gompertz inputs
est <- tibble(
  dt = dt,
  shape = rep(a, 24),
  rate = rep(b, 24)
  )

# Map hazard and chazard
haz <- est %>% 
  mutate(
    haz = pmap_dbl(
      list(dt, shape, rate),
      ~ hgompertz(x = ..1, shape = ..2, rate = ..3)
    ),
    chaz = pmap_dbl(
      list(dt, shape, rate),
      ~ Hgompertz(x = ..1, shape = ..2, rate = ..3)
    )
  )

# Plot
haz %>% 
  ggplot(aes(x = seq(1:24), y = haz)) +
  geom_line() +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Hazard rate")

haz %>% 
  ggplot(aes(x = dt, y = chaz)) +
  geom_line() +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Cumulative hazard")

# Trying with uncertainty ------------------------------------------------------

# Import model to get draws
a_imp_mdl <- readRDS("mars/design_model_mk4/outputs/design4_fit.rds")

# Get draws
draws_df <- as_draws_matrix(a_imp_mdl)
draws_subset <- subset_draws(draws_df, variable = c("mu1", "a", "beta"))

# Get parameter means and variance-covariance across data points
draws_mean <- colMeans(draws_subset)
draws_vcov <- cov(draws_subset)
dim(draws_vcov)

# Sample values from a multivariate normal distribution
library(MASS)
draws_sample <- mvrnorm(100, mu = draws_mean, Sigma = draws_vcov)

# Get tidy dataframe of draws
draws_tidy <- as.data.frame(draws_sample) %>% 
  mutate(sample = row_number()) %>% 
  pivot_longer(-sample, names_to = "var_point", values_to = "sample_value") %>% 
  separate(var_point, into = c("variable", "data_point"), sep = "\\[") %>% 
  mutate(data_point = gsub("\\]", "", data_point))

# Matrix of user inputs
user_input <- t(input$scaling)

# Matrix of mu1 and betas row-means in each sample
sample_betas <- as.matrix(
  draws_tidy %>% 
    filter(grepl("^beta", var_point)) %>% 
    pivot_wider(names_from = "var_point", values_from = "sample_value") %>% 
    dplyr::select(-sample)
)
sample_betas <- rowMeans(sample_betas)

sample_mu1 <- as.matrix(
  draws_tidy %>% 
    filter(grepl("^mu1", var_point)) %>% 
    pivot_wider(names_from = "var_point", values_from = "sample_value") %>% 
    dplyr::select(-sample)
)
sample_mu1 <- rowMeans(sample_mu1)

sample_a <- as.matrix(
  draws_tidy %>% 
    filter(grepl("^a", var_point)) %>% 
    pivot_wider(names_from = "var_point", values_from = "sample_value") %>% 
    dplyr::select(-sample)
)
sample_a <- rowMeans(sample_a)

# Check distribution of sample draws for a and b
# Compare real draws and sampled draws

draws_real_df <- as.data.frame(draws_subset) %>% 
  mutate(draw = row_number()) %>% 
  pivot_longer(-draw, names_to = "parameter_point", values_to = "value")

draws_real_df %>% 
  filter(grepl("^a", parameter_point)) %>%
  ggplot(aes(x = value, group = parameter_point)) +
  geom_density()

draws_tidy %>% 
  filter(grepl("^a", var_point)) %>% 
  ggplot(aes(x = sample_value, group = var_point)) +
  geom_density()

draws_tidy %>% 
  filter(grepl("^b", var_point)) %>% 
  ggplot(aes(x = sample_value, group = var_point)) +
  geom_density()

# Vectorised computation of X * beta
user_input_matrix <- matrix(
  rep(t(input$scaling), each = 100),
  ncol = length(user_input),
  byrow = FALSE
)

# Convert to tibble
beta_multip <- rowSums(sample_betas * user_input_matrix) %>% 
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "combo", values_to = "value") %>% 
  group_by(combo) %>% 
  mutate(sample = row_number()) %>% 
  ungroup() %>% 
  select(-combo)

mu1_tibble <- draws_tidy %>% 
  filter(grepl("^mu1", var_point)) %>% 
  select(sample, mu1_point = var_point, mu1 = sample_value)

a_tibble <- draws_tidy %>% 
  filter(grepl("^a\\[", var_point)) %>% 
  select(sample, a_point = var_point, a = sample_value)

# Join and compute b
parameters_tibble <- beta_multip %>% 
  left_join(mu1_tibble) %>% 
  mutate(b = exp(mu1 + value))

# Get gompertz parameters
gomp_tibble <- parameters_tibble %>% 
  mutate(
    parameter = "b",
    data_point = as.integer(str_extract(mu1_point, "(?<=\\[)\\d+(?=\\])"))
  ) %>% 
  select(sample, parameter, data_point, value = b) %>% 
  full_join(
    a_tibble %>% 
      mutate(
        parameter = "a",
        data_point = as.integer(str_extract(a_point, "(?<=\\[)\\d+(?=\\])"))
      ) %>% 
      select(sample, parameter, data_point, value = a)
  ) %>% 
  arrange(sample, data_point, parameter)

# Summarise mean across data points
gomp_mean <- gomp_tibble %>% 
  group_by(sample, parameter) %>% 
  reframe(mean_value = mean(value)) %>% 
  pivot_wider(names_from = parameter, values_from = mean_value)

# Prepare for estimating uncertainty
haz_grid <- crossing(
  sample = gomp_mean$sample,
  dt = dt
) %>% 
  mutate(
    a = gomp_mean$a[sample],
    b = exp(gomp_mean$b[sample]),
    haz = pmap_dbl(
      list(dt, a, b),
      ~ hgompertz(x = ..1, shape = ..2, rate = ..3)
    ),
    chaz = pmap_dbl(
      list(dt, a, b),
      ~ Hgompertz(x = ..1, shape = ..2, rate = ..3)
    )
  ) %>% 
  group_by(sample) %>% 
  filter(!any(haz > 10)) %>% 
  ungroup()

haz_grid %>% 
  ggplot(aes(x = dt, y = haz, group = sample)) +
  geom_line() +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Hazard rate")

haz_grid %>% 
  ggplot(aes(x = dt, y = chaz, group = sample)) +
  geom_line() +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Cumulative hazard")

# Summarise and plot
haz_sum <- haz_grid %>% 
  group_by(dt) %>% 
  reframe(
    median_haz = median(haz),
    lower_haz = quantile(haz, 0.025),
    upper_haz = quantile(haz, 0.975),
    median_chaz = median(chaz),
    lower_chaz = quantile(chaz, 0.025),
    upper_chaz = quantile(chaz, 0.975)
  )

haz_sum %>% 
  ggplot(aes(x = dt, y = median_haz)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_haz, ymax = upper_haz), alpha = 0.2) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Hazard rate")

haz_sum %>% 
  ggplot(aes(x = dt, y = median_chaz)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_chaz, ymax = upper_chaz), alpha = 0.2) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw() +
  labs(x = "Trial duration in weeks", y = "Cumulative hazard")













predict_attrition <- function(new_X, posterior_samples, dt) {
  
  mu1 <- mean(posterior_samples$mu1)
  a <- mean(posterior_samples$a)
  beta <- beta$beta
  
  n_draws <- length(beta)
  beta <- as.matrix(beta)
  risk_estimates <- matrix(NA, n_draws, length(dt))
  
  b <- exp(posterior_samples$mu1[i] + new_X %*% posterior_samples$beta[i, ])
  estimate <- hgompertz(dt, shape = a, rate = exp(b))
  
  return(risk_estimates)
}



predict_attrition(new_X, a_imp_design, 0.5)
