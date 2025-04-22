
# Inspect outputs
source("scripts/00_functions_packages_etc.R")

library(shinystan)

# Import model fit from MARS
a_fit <- readRDS("outputs/fit.rds")

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
    variable = if_else(variable %in% c("beta", "lambda", "z"), paste0(variable, trial), variable)
  ) %>% 
  group_by(variable) %>% 
  mutate(trial = case_when(n() != 369 ~ NA, TRUE ~ trial)) %>% 
  ungroup()

# Get trial summaries (gompertz, binomial)
trial_summary <- readRDS("processed_data/df_fu_n348.rds") %>% 
  select(trial_id) %>% 
  mutate(trial = row_number()) %>% 
  left_join(
    tidy_summary %>%
      filter(!is.na(trial)) %>% 
      mutate(trial = as.numeric(trial))
  )

# Get global summaries
global_summary <- tidy_summary %>% filter(is.na(trial)) %>% select(-trial)

saveRDS(trial_summary, "outputs/posterior_sum_trials.rds")
saveRDS(global_summary, "outputs/posterior_sum_global.rds")

# Summarise Rhat for collab protocol
rhat_summary <- summary(a_fit)$summary[, "Rhat"] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter") %>% 
  rename(rhat = ".") %>% 
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("[0-9]", "", parameter)
  ) %>% 
  group_by(parameter) %>% 
  reframe(
    n = n(),
    min_rhat = min(rhat, na.rm = T),
    max_rhat = max(rhat, na.rm = T),
    mean_rhat = mean(rhat, na.rm = T),
    med_rhat = median(rhat, na.rm = T)
  )

rhat_summary %>% 
  select(parameter, med_rhat, min_rhat, max_rhat) %>% 
  filter(
    parameter %in% c("a", "a_raw", "b", "beta", "beta_arm", "mu", "rate")
  ) %>% 
  write_csv("outputs/rhat_sum_protocol.csv")
