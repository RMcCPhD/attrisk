
source("mars/intercept_model2/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/intercept_model2/stan_dat.rds")
# a_imp_mdl <- readRDS("mars/intercept_model2/Outputs/int2_fit.rds")
a_imp_sum <- readRDS("mars/intercept_model2/outputs/closest/int2_sum_draws.rds")
a_imp_gomp <- read_csv("mars/intercept_model2/gomp_prior_params.csv")

# Individual data point recovery -----------------------------------------------

# Build a simulated dataset
set.seed(123)

N <- 361
z_b <- rnorm(N)
b_mu1 <- -1.42
tau_mu1 <- 0.54
b_log_sim <- b_mu1 + tau_mu1 * z_b

a_fixed <- -3.85
dt <- a_imp_df$dt
n <- a_imp_df$N

gompertz_prob <- function(b, a, dt) {
  ifelse(abs(a) < 1e-4,
         1 - exp(-b * dt),
         1 - exp(-b / a * (exp(a * dt) - 1)))
}

b_rate_sim <- exp(b_log_sim)
p <- gompertz_prob(b_rate_sim, a_fixed, dt)
r <- rbinom(N, size = n, prob = p)

sim_data <- data.frame(
  trial = 1:N,
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

sim_data$b_log_recovered <- exp(rec_data$median)

# Correlation
cor(sim_data$b_log, sim_data$b_log_recovered)

# Plot
sim_data %>% 
  ggplot(aes(x = b_log, y = b_log_recovered)) +
  geom_point(alpha = 0.6) +
  theme_bw()

# Get tidy mean posterior draws for a and b
b_gomp <- a_imp_sum %>% 
  filter(grepl("^a\\[|^mu1\\[", variable)) %>% 
  separate(variable, into = c("variable", "row"), sep = "\\[") %>% 
  mutate(
    row = gsub("\\]", "", row),
    row = as.numeric(row)
  )

b_gomp %>% group_by(variable) %>% reframe(mean = mean(mean))

set.seed(123)

b_shape <- data.frame(
  centre = rnorm(100000, mean = -3.85, sd = 1.31),
  scale = rnorm(100000, mean = 1.31, sd = 2.62),
  latent = rnorm(100000, mean = 0, sd = 1)
  ) %>% 
  mutate(a_intercept = centre + (scale * latent))

# hist(b_shape$a_intercept)

b_rate <- data.frame(
  centre = rnorm(100000, mean = -1.42, sd = 0.54),
  scale = rnorm(100000, mean = 0.54, sd = 1.08),
  latent = rnorm(100000, mean = 0, sd = 1)
  ) %>% 
  mutate(b_intercept = centre + (scale * latent))

# hist(b_rate$b_intercept)

gomp_raw <- data.frame(
  a_raw = b_shape$a_intercept,
  b_raw = b_rate$b_intercept
  ) %>% 
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# gomp_raw %>% 
#   ggplot(aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free") +
#   scale_x_continuous(n.breaks = 6)

gomp_dat <- b_gomp %>%
  select(variable, value = mean) %>% 
  full_join(gomp_raw) %>% 
  mutate(
    variable = case_match(
      variable,
      "a" ~ "a_recovered",
      "mu1" ~ "b_recovered",
      "a_raw" ~ "a_simulated",
      "b_raw" ~ "b_simulated"
    ),
    parameter = if_else(grepl("^a", variable), "a", "b")
  )

gomp_dat %>% 
  ggplot(aes(x = value, colour = variable)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()
