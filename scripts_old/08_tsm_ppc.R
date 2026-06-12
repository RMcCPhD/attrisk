
# Get model outputs after validating with posterior checks
source("scripts/00_functions_packages_etc.R")

# Import models
a_design <- readRDS("outputs/design_fit.rds")
a_trial <- readRDS("outputs/trial_submodel_fit.rds")

# Get summaries
b_design_summary <- summarise_draws(
  as_draws_df(a_design),
  mean, median, sd,
  ~ quantile(.x, probs = c(0.025, 0.975))
)

b_trial_summary <- summarise_draws(
  as_draws_df(a_trial),
  mean, median, sd,
  ~ quantile(.x, probs = c(0.025, 0.975))
)

# Get draws
b_design_draws_pos <- as_draws_matrix(a_design)
b_trial_draws_pos <- as_draws_matrix(a_trial)

# Save
write_csv(b_design_summary, "outputs/design_model/design_summary.csv")
write_csv(b_trial_summary, "outputs/trial_submodel/trial_summary.csv")

saveRDS(b_design_draws_pos, "outputs/design_model/design_matrix.rds")
saveRDS(b_trial_draws_pos, "outputs/trial_submodel/trial_matrix.rds")
