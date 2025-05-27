
source("mars/intercept_model4/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/intercept_model5/stan_dat.rds")
a_imp_sum <- readRDS("mars/intercept_model5/outputs/int_sum_draws.rds")

# Build a simulated dataset
set.seed(123)

N <- 361
z_b <- rnorm(N)
b_mu1 <- rnorm(N, mean = -1.42, sd = 0.54)
tau_mu1 <- rnorm(N, mean = 0.54, sd = 1.08)
b_log_sim <- b_mu1 + tau_mu1 * z_b

z_a <- rnorm(N)
a_mu1 <- rnorm(N, mean = -3.85, sd = 1.31)
tau_mu1 <- rnorm(N, mean = 1.31, sd = 2.62)
a_sim <- a_mu1 + tau_mu1 * z_a

dt <- a_imp_df$dt
n <- a_imp_df$N

gompertz_prob <- function(b, a, dt) {
  ifelse(
    abs(a) < 1e-4,
    1 - exp(-b * dt),
    1 - exp(-b / a * (exp(a * dt) - 1))
  )
}

b_sim <- exp(b_log_sim)
p <- gompertz_prob(b_sim, a_sim, dt)
r <- rbinom(N, size = n, prob = p)

sim_data <- data.frame(
  trial = 1:N,
  a_sim = a_sim,
  b_sim = b_sim,
  b_log_sim = b_log_sim,
  dt = dt,
  N = n,
  r = r
)

# Recovered dataset
rec_data <- a_imp_sum %>% 
  filter(grepl("^mu1\\[", variable)) %>% 
  separate(variable, into = c("variable", "row"), sep = "\\[") %>% 
  mutate(
    row = gsub("\\]", "", row),
    row = as.numeric(row)
  )

rec_data_a <- a_imp_sum %>% 
  filter(grepl("^a\\[", variable)) %>% 
  separate(variable, into = c("variable", "row"), sep = "\\[") %>% 
  mutate(
    row = gsub("\\]", "", row),
    row = as.numeric(row)
  )

sim_data$a_recovered <- rec_data_a$median
sim_data$b_log_recovered <- rec_data$median

# Correlation
cor(sim_data$a_sim, sim_data$a_recovered)
cor(sim_data$b_log_sim, sim_data$b_log_recovered)

# Plot
sim_data %>% 
  ggplot(aes(x = a_sim, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "a_simulated", y = "a_recovered") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 6)

sim_data %>% 
  ggplot(aes(x = b_log_sim, y = b_log_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "b_log_simulated", y = "b_log_recovered") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 6)
