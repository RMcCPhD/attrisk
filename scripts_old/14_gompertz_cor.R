
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
    shape_scaled = shape * 365,
    rate_scaled = rate + log(365)
  )

b_tidy_rates %>% reframe(across(rate:rate_scaled, ~ unlist(list(mean = mean(.), sd = sd(.)))))

# Test scaling
plot(hgompertz(1:365, rate = exp(-6.95), shape = -0.00491))
plot(hgompertz(seq(0, 1, by = 0.1), rate = exp(-1.05), shape = -1.79))
cor(b_tidy_rates$shape_scaled, b_tidy_rates$rate_scaled)

b_tidy_rates %>% count(condition)

# Linear regression model
c_fit <- lm(shape_scaled ~ rate_scaled, data = b_tidy_rates)
d_fit <- lm(rate_scaled ~ shape_scaled, data = b_tidy_rates)

summary(c_fit)
summary(d_fit)

# Add residuals to tidy data
b_tidy_rates$residual_shape <- broom::augment(c_fit)$.resid
b_tidy_rates$residual_rate <- broom::augment(d_fit)$.resid

hist(b_tidy_rates$residual_shape)
hist(b_tidy_rates$residual_rate)

summary(b_tidy_rates$residual_shape)
sd(b_tidy_rates$residual_shape)

summary(b_tidy_rates$residual_rate)
sd(b_tidy_rates$residual_rate)

mad(b_tidy_rates$residual_shape)
mad(b_tidy_rates$residual_rate)

# Contour plot
plot_density <- b_tidy_rates %>% 
  ggplot(aes(x = rate_scaled, y = shape_scaled)) +
  geom_density_2d() +
  theme_bw() +
  labs(x = "Log-rate", y = "Shape")

ggsave(
  "mpes/plot_density.png",
  plot_density,
  height = 4,
  width = 4,
  units = "in"
)

# Parameter correlation
plot_cor <- b_tidy_rates %>% 
  ggplot(aes(x = rate_scaled, y = shape_scaled)) +
  geom_point() +
  theme_bw() +
  labs(x = "Log-rate", y = "Shape") +
  scale_y_continuous(n.breaks = 6)

ggsave(
  "mpes/plot_cor.png",
  plot_cor,
  height = 4,
  width = 4,
  units = "in"
)

# Parameter distribution
plot_density_gomp <- b_tidy_rates %>% 
  gather("parameter", "est", rate_scaled, shape_scaled) %>% 
  mutate(
    parameter = case_match(
      parameter, 
      "rate_scaled" ~ "rate", "shape_scaled" ~ "shape"
    )
  ) %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw() +
  labs(x = "Estimates (rate as log-rate)", y = "Density")

ggsave(
  "mpes/plot_density_gomp.png",
  plot_density_gomp,
  height = 4,
  width = 6,
  units = "in"
)

# Residual deviance for shape
plot_residuals <- b_tidy_rates %>% 
  gather("parameter", "est", residual_rate, residual_shape) %>% 
  mutate(
    parameter = case_match(
      parameter, 
      "residual_rate" ~ "rate", 
      "residual_shape" ~ "shape"
    )
  ) %>% 
  ggplot(aes(x = est)) +
  geom_histogram() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()  + 
  labs(x = "Residual deviance", y = "Density")

ggsave(
  "mpes/plot_residuals.png",
  plot_residuals,
  height = 4,
  width = 6,
  units = "in"
)

# Testing samples
set.seed(123)

test_shape_int <- data.frame(
  scale_centre = rnorm(1000, mean = -0.034, sd = 0.04),
  scale_scale = rnorm(1000, mean = -10.69, sd = 1.39),
  scale_latent = rnorm(1000, mean = 0, sd = 0.025)
  ) %>% 
  mutate(shape_int = scale_centre + (scale_scale * scale_latent))

hist(test_shape_int$shape_int)

test_rate_int <- data.frame(
  rate_centre = rnorm(1000, mean = -5.01, sd = 0.73),
  rate_rate = rnorm(1000, mean = -0.037, sd = 0.004),
  rate_latent = rnorm(1000, mean = 0, sd = 0.034)
) %>% 
  mutate(rate_int = rate_centre + (rate_rate * rate_latent))

hist(test_rate_int$rate_int)

# Dataframe of outputs needed for informed priors
shape_inf_priors <- data.frame(
  type = "shape",
  centre_mean = -0.034,
  centre_sd = 0.04,
  scale_mean = -10.69,
  scale_sd = 1.39,
  latent_mean = 0,
  latent_sd = 0.025
)

rate_inf_priors <- data.frame(
  type = "rate",
  centre_mean = -5.01,
  centre_sd = 0.73,
  scale_mean = -0.037,
  scale_sd = 0.004,
  latent_mean = 0,
  latent_sd = 0.034
)

z_inf_priors <- shape_inf_priors %>% full_join(rate_inf_priors)
write_csv(z_inf_priors, "mars/intercept_model2/gomp_prior_params.csv")
