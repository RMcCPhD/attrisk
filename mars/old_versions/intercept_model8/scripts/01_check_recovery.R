# Load packages and config
source("mars/intercept_model7/scripts/00_config_functions_packages.R")

# Import dataset to build simulated data
a_imp_df <- readRDS("mars/intercept_model8/stan_dat.rds")

# Set seed and number of data points
set.seed(123)
N <- 361

# 1. Simulate parameters
s_true <- abs(rnorm(1, mean = 0.3, sd = 0.2))
mu_mu1_true <- -1.05                         
z_sim <- rnorm(N)                            
b_log_sim <- mu_mu1_true + s_true * z_sim    
a_sim <- -3.85 + (-0.21) * b_log_sim         
b_sim <- exp(b_log_sim)

# 2. Simulate attrition outcomes
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

# 3. Assemble simultated tibble
sim_data <- tibble(
  trial = 1:N,
  z_sim = z_sim,
  b_log_sim = b_log_sim,
  a_sim = a_sim,
  b_sim = b_sim,
  dt = dt,
  N = n,
  r = r
)

# 4. Extract posterior summaries of parameters
a_imp_sum <- readRDS("mars/intercept_model8/outputs/int_sum_draws.rds")

rec_z <- a_imp_sum %>%
  filter(str_detect(variable, "^z_out\\[")) %>%
  mutate(
    trial = as.numeric(str_extract(variable, "\\d+")),
    z_recovered = median
  ) %>%
  select(trial, z_recovered)

rec_mu_mu1 <- a_imp_sum %>%
  filter(variable == "mu_mu1") %>%
  pull(median)

rec_beta_0 <- a_imp_sum %>%
  filter(variable == "beta_0") %>%
  pull(median)

rec_beta_1 <- a_imp_sum %>%
  filter(variable == "beta_1") %>%
  pull(median)

# 5. Reconstruct recovered log-rate and shape
sim_data2 <- sim_data %>%
  left_join(rec_z, by = "trial") %>%
  mutate(
    b_log_recovered = rec_mu_mu1 + z_recovered,
    a_recovered = rec_beta_0 + rec_beta_1 * b_log_recovered
  )

# Compare observed to predicted counts
r_new <- a_imp_sum %>%
  filter(str_detect(variable, "^r_new\\[")) %>%
  mutate(
    trial = as.numeric(str_extract(variable, "\\d+")),
    r_new = median
  ) %>%
  select(trial, r_new)

summary(sim_data2$r)
hist(sim_data2$r)
hist(a_imp_df$r)
hist(r_new$r_new)

# 6. Get correlations
cor(sim_data2$a_sim, sim_data2$a_recovered)
cor(sim_data2$b_log_sim, sim_data2$b_log_recovered)
cor(sim_data2$z_sim, sim_data2$z_recovered)

# 7. Visualise correlations
ggplot(sim_data2, aes(x = a_sim, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Shape: Simulated vs Recovered", x = "Simulated", y = "Recovered") +
  theme_bw()

ggplot(sim_data2, aes(x = b_log_sim, y = b_log_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Log-Rate: Simulated vs Recovered", x = "Simulated", y = "Recovered") +
  theme_bw()

ggplot(sim_data2, aes(x = z_sim, y = z_recovered)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Latent z: Simulated vs Recovered", x = "Simulated", y = "Recovered") +
  theme_bw()

