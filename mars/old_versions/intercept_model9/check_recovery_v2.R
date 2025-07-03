
rm(list = ls())

library(tidyverse)
library(MASS)

# Settings
N <- 361
mu <- c(-1.21, -0.671)
sigma <- c(0.594, 1.5)
rho <- -0.34

# Build correlation matrix
Sigma_corr <- matrix(c(1, rho, rho, 1), nrow = 2)
L_corr <- t(chol(Sigma_corr))  # Cholesky factor of correlation matrix

# Simulate latent z ~ N(0, 1)
set.seed(123)
z <- matrix(rnorm(2 * N), nrow = 2)

# Generate theta_std = L_corr * z
theta_std <- L_corr %*% z

# Generate true mu1 and a
mu1_sim <- mu[1] + sigma[1] * theta_std[1, ]
a_sim   <- mu[2] + sigma[2] * theta_std[2, ]

posterior_sum <- readRDS("mars/intercept_model9/int_sum_draws.rds")

rec_mu1 <- posterior_sum %>%
  filter(str_detect(variable, "^mu1\\[")) %>%
  mutate(trial = row_number(),
         mu1_recovered = median) %>%
  dplyr::select(trial, mu1_recovered)

rec_a <- posterior_sum %>%
  filter(str_detect(variable, "^a\\[")) %>%
  mutate(trial = parse_number(variable),
         a_recovered = median) %>%
  dplyr::select(trial, a_recovered)

recovery_df <- tibble(
  trial = 1:N,
  mu1_sim = mu1_sim,
  a_sim = a_sim
) %>%
  left_join(rec_mu1, by = "trial") %>%
  left_join(rec_a, by = "trial")

cat("Correlation (log-rate): ", cor(recovery_df$mu1_sim, recovery_df$mu1_recovered), "\n")
cat("Correlation (shape):    ", cor(recovery_df$a_sim, recovery_df$a_recovered), "\n")

ggplot(recovery_df, aes(x = mu1_sim, y = mu1_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Log-Rate Recovery", x = "Simulated", y = "Recovered") +
  theme_bw()

ggplot(recovery_df, aes(x = a_sim, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Shape Recovery", x = "Simulated", y = "Recovered") +
  theme_bw()
