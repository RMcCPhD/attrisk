
source("mars/intercept_model6/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/intercept_model6/stan_dat.rds")
a_imp_sum <- readRDS("mars/intercept_model6/outputs/int_sum_draws.rds")

# Set seed and number of data points
set.seed(123)
N <- 361

# Cholesky factor for correlation between log(rate) and shape
L <- t(chol(matrix(c(1, -0.63, -0.63, 1), nrow = 2)))

# Simulate latent z
z_latent <- matrix(rnorm(N * 2), ncol = 2)

# Apply Cholesky to induce correlation
theta_raw <- z_latent %*% L 

# Empirical priors
mu <- c(-1.05, -1.79)
tau <- c(0.73, 2.26)

# Parameter matrix [log(rate), shape]
theta <- sweep(theta_raw, 2, tau, `*`) + matrix(mu, nrow = N, ncol = 2, byrow = TRUE)

# Extract simulated values
b_log_sim <- theta[, 1]
a_sim <- theta[, 2]
b_sim <- exp(b_log_sim)

# Generate simulated attrition probabilities and counts
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

sim_data <- tibble(
  trial = 1:N,
  a_sim = a_sim,
  b_sim = b_sim,
  b_log_sim = b_log_sim,
  dt = dt,
  N = n,
  r = r
)

# Extract recovered estimates from STAN summarised draws
rec_b <- a_imp_sum %>% 
  filter(str_detect(variable, "^mu1\\[")) %>%
  mutate(trial = row_number())

rec_a <- a_imp_sum %>%
  filter(str_detect(variable, "^a\\[")) %>%
  mutate(trial = as.numeric(str_extract(variable, "\\d+")))

# Join with simulated data
sim_data <- sim_data %>%
  left_join(rec_b %>% select(trial, b_log_recovered = median), by = "trial") %>%
  left_join(rec_a %>% select(trial, a_recovered = median), by = "trial")

# Correlation diagnostics
cor_shape <- cor(sim_data$a_sim, sim_data$a_recovered)
cor_lograte <- cor(sim_data$b_log_sim, sim_data$b_log_recovered)

cat("Correlation (shape):", round(cor_shape, 3), "\n")
cat("Correlation (log-rate):", round(cor_lograte, 3), "\n")

# RMSE and bias
rmse_shape <- sqrt(mean((sim_data$a_sim - sim_data$a_recovered)^2))
rmse_lograte <- sqrt(mean((sim_data$b_log_sim - sim_data$b_log_recovered)^2))

bias_shape <- mean(sim_data$a_recovered - sim_data$a_sim)
bias_lograte <- mean(sim_data$b_log_recovered - sim_data$b_log_sim)

cat("RMSE (shape):", round(rmse_shape, 3), "\n")
cat("Bias (shape):", round(bias_shape, 3), "\n")
cat("RMSE (log-rate):", round(rmse_lograte, 3), "\n")
cat("Bias (log-rate):", round(bias_lograte, 3), "\n")

# Plots
p_shape <- ggplot(sim_data, aes(x = a_sim, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "Simulated Shape", y = "Recovered Shape") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 6)

p_lograte <- ggplot(sim_data, aes(x = b_log_sim, y = b_log_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "Simulated Log-Rate", y = "Recovered Log-Rate") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 6)
