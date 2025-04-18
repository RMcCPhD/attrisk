
# Inspect outputs
source("scripts/00_functions_packages_etc.R")

library(shinystan)

# Import model fit from MARS
a_fit <- readRDS("outputs/fit2.rds")

# Extract posterior samples
posterior_samples <- as_draws_df(a_fit)
posterior_summary <- summarise_draws(posterior_samples, default_summary_measures())

posterior_summary <- summarise_draws(
  posterior_samples,
  measures = list(
    mean = mean,
    median = median,
    q5 = ~quantile(.x, 0.05),
    q97.5 = ~quantile(.x, 0.975)
  )
)

# Summarise
tidy_summary <- posterior_summary %>% 
  rename(
    mean = measures.mean,
    median = measures.median,
    q5 = `5%`,
    q975 = `97.5%`
  ) %>% 
  separate(
    variable,
    into = c("variable", "trial"),
    sep = "\\["
  ) %>% 
  mutate(
    trial = gsub("\\]", "", trial),
    variable = if_else(variable == "beta", paste0(variable, trial), variable)
  )

# Join trial ids
trial_summary <- readRDS("processed_data/df_fu_n348.rds") %>% 
  select(trial_id) %>% 
  mutate(trial = row_number()) %>% 
  left_join(tidy_summary %>% mutate(trial = as.numeric(trial)))

# Get global parameter summaries
global_summary <- tidy_summary %>% filter(is.na(trial)) %>% select(-trial)

saveRDS(trial_summary, "outputs/posterior_sum_trials.rds")
saveRDS(global_summary, "outputs/posterior_sum_global.rds")



