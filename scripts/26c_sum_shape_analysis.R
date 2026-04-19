
# Summaries for the external analysis and modelling of shape on rate

source("scripts/00_config.R")
library(tidyverse)

# Code below copied from 20_model_shape.R to setup tidy results ----------------

# Import dataset
a_imp_df <- readRDS("processed_data/latest_dataset.rds")

# Trials of interest
z_ipd_trials <- a_imp_df %>% filter(source == "ipd") %>% distinct(trial_id)
z_ipd_trials <- as.vector(z_ipd_trials$trial_id)

# Import vivli intercept-only model results
a_imp_mdl <- read_csv("data/attrition_intercept_models/tidy_res.csv")

# Import intercept-only models from attrition rates paper
a_imp_mdl_rates <- read_csv("data/attrition_rates_data/parameters.csv")

# Prepare newer model results
# Was modelled in years, no scaling needed
b_mdl_res <- a_imp_mdl %>% 
  rename(
    trial_id = nct_id,
    shape = shape_est,
    rate = rate_est
  ) %>% 
  select(-c(shape_se, rate_se))

# All present as expected
sum(b_mdl_res$trial_id %in% z_ipd_trials)

# Prepare rates model data
# Remove trials in newer set
# Scale to duration in years
c_tidy_rates <- a_imp_mdl_rates %>% 
  filter(dist == "gompertz" & parameter %in% c("rate", "shape")) %>% 
  select(-c(dist, se)) %>% 
  pivot_wider(names_from = "parameter", values_from = "est") %>% 
  rename(trial_id = ctgov) %>% 
  filter(
    !(trial_id %in% b_mdl_res$trial_id),
    trial_id %in% z_ipd_trials
  ) %>% 
  mutate(
    shape = shape * 365,
    rate = rate + log(365)
  )

# 4 additional trials present
sum(c_tidy_rates$trial_id %in% z_ipd_trials)

# Join model datasets, 74 trials
d_join_mdl <- b_mdl_res %>% 
  full_join(c_tidy_rates) %>% 
  inner_join(
    read_csv("vivli/01_res.csv") %>% 
      distinct(nct_id) %>% 
      rename(trial_id = nct_id)
  )

# Join characteristics
# Transform rate
e_join_chars <- d_join_mdl %>% left_join(a_imp_df)

# Correlation
cor.test(e_join_chars$shape, e_join_chars$rate, method = "pearson")

# Plot distributions -----------------------------------------------------------

# Parameter distribution
plot_gomp <- d_join_mdl %>% 
  gather("parameter", "est", shape, rate) %>% 
  mutate(
    parameter = case_match(
      parameter,
      "shape" ~ "Shape",
      "rate" ~ "Log-rate"
    )
  ) %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_classic() +
  labs(x = "Model estimates", y = "Density")

ggsave(
  "outputs/plots/design_expl_shape_rate_distr.png",
  plot_gomp,
  width = 6,
  height = 4,
  units = "in"
)

# Plot correlation -------------------------------------------------------------

plot_cor <- d_join_mdl %>% 
  ggplot(aes(x = shape, y = rate)) +
  geom_point() +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "Pearson's r = ", "-0.40",
      "\np = ", "<0.01"
    ),
    hjust = 1.1, vjust = 1.5
  ) +
  theme_classic() +
  labs(x = "Shape", y = "Log-rate")

ggsave(
  "outputs/plots/design_expl_shape_rate_cor.png",
  plot_cor,
  width = 6,
  height = 4,
  units = "in"
)


# Fit shape model --------------------------------------------------------------

# Data - duration removed
f_design_df <- e_join_chars %>% 
  select(-c(source, bmi_lower:bmi_upper, contains("dbl"), n_log, duration))
  # mutate(rate = exp(rate))

# Prepare model formula
z_covars_names <- names(f_design_df[c(3, 8:26)])
z_covars <- paste0(z_covars_names, collapse = " + ")
z_fml <- as.formula(paste0("shape ~ ", z_covars))

# Construct glmnet formula and design matrices
z_cv_fml <- as.formula(paste0("~ ", z_covars))
design_X <- model.matrix(z_cv_fml, data = f_design_df)[, -1]
design_y <- f_design_df$shape

# Fit linear model with shrinkage via glmnet
cv_fit <- glmnet::cv.glmnet(design_X, design_y)
plot(cv_fit)

# Get betas with lowest lambda
design_betas <- as.numeric(coef(cv_fit, s = "lambda.1se")[-1, ])

# Get residual standard deviation
cv_pred <- predict(cv_fit, newx = design_X, s = "lambda.1se")
cv_resid <- design_y - cv_pred
cv_n <- length(cv_resid)
cv_p <- length(coef(cv_fit))
cv_resid_sd <- sqrt(sum(cv_resid^2) / (cv_n - cv_p))

# Fit linear model without shrinkage to explore
# Rate, maximum age (eligibility), maximium bmi (eligibility), phase and 
# quadruple blinding (compared to open-label) predict shape
f_mdl <- lm(z_fml, data = f_design_df)
summary(f_mdl)

# Compute residual standard dev to compare to RMSE from glmnet residuals
f_resid <- residuals(f_mdl)
f_n <- length(f_resid)
f_p <- length(coef(f_mdl))
f_resid_sd <- sqrt(sum(f_resid^2) / (f_n - f_p))

# Diagnostics for linear model
f_diag <- broom::glance(f_mdl)

# Tidy outputs
f_res <- broom::tidy(f_mdl) %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  )

# Generate tidy outputs --------------------------------------------------------

# Tidy covariate names
var_names <- read_csv("scratch_data/tidy_covar_names.csv") %>% 
  mutate(
    covar_name = factor(
      covar_name,
      levels = c(
        "Treatment type",
        "Eligibility: Lower age limit", "Eligibility: Upper age limit", 
        "Eligibility: Lower BMI limit", "Eligibility: Upper BMI limit", 
        "Eligibility: Lower eGFR limit", "Eligibility: Upper eGFR limit",
        "Trial size", "Trial phase", "Double-blinded", "Triple-blinded", 
        "Quadruple-blinded", "Number of arms", "Number of facilities", 
        "Number of POMs",
        "African regions", "Regions of the Americas", 
        "Eastern Mediterranean regions", "European regions", 
        "Western Pacific regions", "Southeast Asian regions"
      )
    )
  )

var_levels <- c(
  "Log-rate",
  "Treatment type", 
  "Eligibility: Lower age limit", "Eligibility: Upper age limit", 
  "Eligibility: Lower BMI limit", "Eligibility: Upper BMI limit", 
  "Eligibility: Lower eGFR limit", "Eligibility: Upper eGFR limit", 
  "Trial size", "Trial phase", "Double-blinded", "Triple-blinded", 
  "Quadruple-blinded", "Number of arms", "Number of facilities", 
  "Number of POMs", "African regions", "Regions of the Americas", 
  "Eastern Mediterranean regions", "European regions", "Western Pacific regions", 
  "Southeast Asian regions"
)

# Plot coefficients
plot_coefs <- f_res %>% 
  filter(!term == "(Intercept)") %>% 
  left_join(
    var_names %>% 
      rename(term = covar) %>% 
      distinct(term, covar_name)
  ) %>% 
  mutate(
    covar_name = if_else(is.na(covar_name), "Log-rate", covar_name),
    covar_name = factor(covar_name, levels = var_levels)
  ) %>% 
  ggplot(aes(x = fct_rev(covar_name), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(limits = c(-8, 8), n.breaks = 9) +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = 0, colour = "red") +
  labs(x = NULL, y = "Mean effect on Gompertz shape (95% confidence intervals)")

ggsave(
  "outputs/plots/design_expl_shape_coefs.png",
  plot_coefs,
  width = 6,
  height = 6,
  units = "in"
)

# Get predicted shape and residuals
f_resid <- broom::augment(f_mdl) %>% 
  as_tibble() %>% 
  rename(
    shape_pred = .fitted,
    resid = .resid
  ) %>% 
  select(shape, rate, shape_pred, resid)

# Plot residuals
plot_resid <- f_resid %>% 
  ggplot(aes(x = resid)) + 
  geom_density() + 
  theme_classic() +
  labs(x = "Residual variation", y = "Density")

ggsave(
  "outputs/plots/design_expl_mdl_residuals.png",
  plot_resid,
  width = 4,
  height = 4,
  units = "in"
)

# Predicted vs observed
plot_pred_obs <- f_resid %>% 
  ggplot(aes(x = shape_pred, y = shape)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  theme_classic() +
  labs(x = "Predicted shape", y = "Observed shape")

ggsave(
  "outputs/plots/design_expl_mdl_pred_vs_obs.png",
  plot_pred_obs,
  width = 4,
  height = 4,
  units = "in"
)

# Residuals vs predicted
plot_pred_resid <- f_resid %>% 
  ggplot(aes(x = shape_pred, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "Predicted shape", y = "Residuals")

# Save results
write_csv(f_res, "outputs/model_shape/design/lm_res.csv")
write_csv(f_diag, "outputs/model_shape/design/lm_diags.csv")
write_csv(f_resid, "outputs/model_shape/design/lm_pred_resid.csv")

ggsave(
  "outputs/model_shape/design/01_lm_res_coefs.png",
  plot_coefs,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  "outputs/model_shape/design/02_lm_residuals.png",
  plot_resid,
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/design/03_lm_obs_pred_shape.png",
  plot_pred_obs,
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/design/03b_lm_resid_obs_pred.png",
  cowplot::plot_grid(
    plot_resid, plot_pred_obs, ncol = 2, labels = c("A", "B")
  ),
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/design/04_lm_pred_resid.png",
  plot_pred_resid,
  width = 6,
  height = 6,
  units = "in"
)

# Fit model without blinding
z_covars_anova <- gsub("\\+ blinding ", "", z_covars)
z_fml_anova <- as.formula(paste0("shape ~ rate + ", z_covars_anova))
f_mdl2 <- lm(z_fml_anova, data = f_design_df)
summary(f_mdl2)

# ANOVA test - significance
# Including blinding does significantly improves model fit (quadruple)
f_anova <- anova(f_mdl2, f_mdl)
summary(f_anova)
