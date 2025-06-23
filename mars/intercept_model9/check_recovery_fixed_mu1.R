library(tidyverse)

# Define gompertz_prob function
gompertz_prob <- function(b, a, dt) {
  ifelse(abs(a) < 1e-4,
         1 - exp(-b * dt),
         1 - exp(-b / a * (exp(a * dt) - 1)))
}

# Grid settings
dt_grid <- c(1, 2, 5)
n_grid <- c(100, 500, 1000)

# Simulation constants
N <- 361
mu1_fixed <- -1.2
b_fixed <- exp(mu1_fixed)
rho <- -0.34

# Placeholder for results
results <- tibble()

# Loop over conditions
for (dt_val in dt_grid) {
  for (n_val in n_grid) {
    
    # Simulate data
    set.seed(123)
    dt <- rep(dt_val, N)
    n <- rep(n_val, N)
    a_sim <- seq(-3, 3, length.out = N)
    p <- gompertz_prob(b_fixed, a_sim, dt)
    r <- rbinom(N, size = n, prob = p)
    
    # Simulated trial data
    sim_data <- tibble(trial = 1:N, a_sim = a_sim, r = r, N = n, dt = dt)
    
    # --- Fit model or load precomputed summaries ---
    # Replace this with model execution or summary load
    # Here we assume summary is saved with pattern "summaries/dt{dt}_n{n}.rds"
    summary_path <- glue::glue("summaries/summary_dt{dt_val}_n{n_val}.rds")
    
    if (!file.exists(summary_path)) {
      message(glue::glue("Skipping: {summary_path} not found."))
      next
    }
    
    posterior_sum <- readRDS(summary_path)
    
    rec_a <- posterior_sum %>%
      filter(str_detect(variable, "^a\\[")) %>%
      mutate(trial = parse_number(variable)) %>%
      rename(a_recovered = median) %>%
      arrange(trial)
    
    # Merge with simulation
    merged <- sim_data %>%
      left_join(rec_a, by = "trial")
    
    # Compute correlation
    cor_val <- cor(merged$a_sim, merged$a_recovered)
    
    results <- results %>%
      add_row(
        dt = dt_val,
        n = n_val,
        correlation = round(cor_val, 3)
      )
  }
}

# --- Plot results ---
ggplot(results, aes(x = factor(dt), y = correlation, fill = factor(n))) +
  geom_col(position = position_dodge(width = 0.6)) +
  labs(
    title = "Shape Recovery Correlation by dt and n",
    x = "Follow-up time (dt)", y = "Correlation (shape: sim vs recovered)",
    fill = "Sample size (n)"
  ) +
  theme_bw()
