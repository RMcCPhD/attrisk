
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)
library(flexsurv)

# Import attrition rates data
a_imp_rates_chars <- read_csv("data/attrition_rates_data/study_chars.csv")
a_imp_rates_pars <- read_csv("data/attrition_rates_data/parameters.csv")
a_imp_rates_vcov <- read_csv("data/attrition_rates_data/vcov_pblc.csv")

# Prepare data
# Parameters scaled to time in weeks
b_tidy_rates <- a_imp_rates_chars %>% 
  select(ctgov, condition) %>% 
  distinct() %>% 
  left_join(
    a_imp_rates_pars %>% 
      filter(dist == "gompertz") %>% 
      select(ctgov, parameter, est) %>% 
      spread(parameter, est)
  ) %>% 
  mutate(
    shape_scaled = shape * 7,
    rate_scaled = rate + log(7)
  )

b_tidy_rates %>% reframe(across(rate:rate_scaled, ~ mean(.)))

# Test scaling
plot(hgompertz(1:365, rate = exp(-6.95), shape = -0.00491))
plot(hgompertz(1:52, rate = exp(-5.01), shape = -0.0344))

b_tidy_rates %>% count(condition)

# Linear regression model
c_fit <- lm(shape_scaled ~ rate_scaled, data = b_tidy_rates)
summary(c_fit)

# Add residuals
b_tidy_rates$residual <- broom::augment(c_fit)$.resid

# Plots
b_tidy_rates %>% 
  gather("parameter", "est", rate, shape) %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw() +
  labs(x = "Estimates (rate as log-rate)", y = "Density")

b_tidy_rates %>% 
  gather("parameter", "est", rate_scaled, shape_scaled) %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw() +
  labs(x = "Estimates (rate as log-rate)", y = "Density")

b_tidy_rates %>% 
  ggplot(aes(x = rate, y = shape)) +
  geom_point() +
  theme_bw() +
  labs(x = "Log-rate", y = "Shape")

b_tidy_rates %>% 
  ggplot(aes(x = rate_scaled, y = shape_scaled)) +
  geom_point() +
  theme_bw() +
  labs(x = "Log-rate", y = "Shape")
