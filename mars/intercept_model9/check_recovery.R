rm(list = ls())
gc()

library(tidyverse)
library(MASS)

# Load STAN input data
a_imp_df <- readRDS("mars/intercept_model9/stan_dat.rds")
a_imp_sum <- readRDS("mars/intercept_model9/int_sum_draws.rds")

# Set constants
set.seed(123)
N <- 361

# --- Step 1: Simulate bivariate Gompertz parameters ---

# Use prior means for mu and sigma
mu_vec <- c(-1.21, -0.671)
sigma_vec <- c(0.594, 1.85)  # SDs
rho <- -0.34                 # Prior mean for correlation

# Construct covariance matrix
Sigma <- matrix(c(
  sigma_vec[1]^2, rho * sigma_vec[1] * sigma_vec[2],
  rho * sigma_vec[1] * sigma_vec[2], sigma_vec[2]^2
), nrow = 2, byrow = TRUE)

# Simulate trial-level (log_rate, shape)
theta_sim <- MASS::mvrnorm(n = N, mu = mu_vec, Sigma = Sigma)
mu1_sim <- theta_sim[, 1]
a_sim   <- theta_sim[, 2]
b_sim   <- exp(mu1_sim)

# --- Step 2: Simulate attrition counts ---

dt <- a_imp_df$dt
n <- a_imp_df$N

gompertz_prob <- function(b, a, dt) {
  ifelse(
    abs(a) < 1e-4,
    1 - exp(-b * dt),
    1 - exp(-b / a * (exp(a * dt) - 1))
  )
}

p <- gompertz_prob(b_sim, a_sim, dt)
r <- rbinom(N, size = n, prob = p)

# --- Step 3: Create simulation tibble ---

sim_data <- tibble(
  trial = 1:N,
  mu1_sim = mu1_sim,
  a_sim = a_sim,
  r = r
)

# --- Step 4: Extract posterior summaries from STAN output ---

rec_mu1 <- a_imp_sum %>%
  filter(str_detect(variable, "^mu1\\[")) %>%
  mutate(trial = row_number()) %>%
  rename(mu1_recovered = median) %>%
  dplyr::select(trial, mu1_recovered)

rec_a <- a_imp_sum %>%
  filter(str_detect(variable, "^a\\[")) %>%
  mutate(trial = parse_number(variable)) %>%
  rename(a_recovered = median) %>%
  dplyr::select(trial, a_recovered)

# --- Step 5: Merge simulated and recovered values ---

rec_data <- sim_data %>%
  left_join(rec_mu1, by = "trial") %>%
  left_join(rec_a, by = "trial")

# --- Step 6: Compute correlation between simulated and recovered ---

cor_log_rate <- cor(rec_data$mu1_sim, rec_data$mu1_recovered)
cor_shape <- cor(rec_data$a_sim, rec_data$a_recovered)

cat("Correlation (log-rate):", round(cor_log_rate, 3), "\n")
cat("Correlation (shape):   ", round(cor_shape, 3), "\n")

# --- Step 7: Plot parameter recovery ---

ggplot(rec_data, aes(x = mu1_sim, y = mu1_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Log-Rate: Simulated vs Recovered",
    x = "Simulated log-rate", y = "Recovered log-rate"
  ) +
  theme_bw()

ggplot(rec_data, aes(x = a_sim, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Shape: Simulated vs Recovered",
    x = "Simulated shape", y = "Recovered shape"
  ) +
  theme_bw()
