
source("mars/intercept_model2/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/intercept_model/stan_dat.rds")
# a_imp_mdl <- readRDS("mars/intercept_model/Outputs/int_fit.rds")
a_imp_sum <- readRDS("mars/intercept_model/outputs/int_sum_draws.rds")
# a_imp_gomp <- read_csv("mars/intercept_model/gomp_prior_params.csv")

# Individual data point recovery -----------------------------------------------

# Build a simulated dataset
set.seed(123)

N <- 361
z_b <- rnorm(N)
b_mu1 <- -1.42
tau_mu1 <- 0.54
b_log_sim <- b_mu1 + tau_mu1 * z_b

z_a <- rnorm(N)
a_mu1 <- 0
tau_mu1 <- 1
a_sim <- a_mu1 + tau_mu1 * z_a

dt <- a_imp_df$dt
n <- a_imp_df$N

gompertz_prob <- function(b, a, dt) {
  ifelse(abs(a) < 1e-4,
         1 - exp(-b * dt),
         1 - exp(-b / a * (exp(a * dt) - 1)))
}

b_rate_sim <- exp(b_log_sim)
p <- gompertz_prob(b_rate_sim, a_sim, dt)
r <- rbinom(N, size = n, prob = p)

sim_data <- data.frame(
  trial = 1:N,
  a_raw = a_sim,
  b_log = b_log_sim,
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
sim_data$b_log_recovered <- exp(rec_data$median)

# Correlation
cor(sim_data$a_raw, sim_data$a_recovered)
cor(sim_data$b_log, sim_data$b_log_recovered)

# Plot
sim_data %>% 
  ggplot(aes(x = a_raw, y = a_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw()

sim_data %>% 
  ggplot(aes(x = b_log, y = b_log_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw()
