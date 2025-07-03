
# Model shape on trial covariates for design and conduct models
source("scripts/00_config.R")
library(tidyverse)

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

# Join model datasets, 76 trials
d_join_mdl <- b_mdl_res %>% full_join(c_tidy_rates)

# Join characteristics
e_join_chars <- d_join_mdl %>% left_join(a_imp_df)

# Plot distributions
plot_gomp <- d_join_mdl %>% 
  gather("parameter", "est", shape, rate) %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  theme_bw() +
  labs(x = "Estimates (rate as log-rate)", y = "Density")

# Correlation between rate and shape
# -0.4
cor(d_join_mdl$shape, d_join_mdl$rate, method = "pearson")

# Contour plot with scatter of shape and rate
plot_contour <- d_join_mdl %>% 
  ggplot(aes(x = rate, y = shape)) +
  geom_density_2d() +
  geom_point(aes(x = rate, y = shape)) +
  theme_bw() +
  labs(x = "Log-rate", y = "Shape")

# Save plots
ggsave(
  "outputs/model_shape/gomp_params_dist.png",
  plot_gomp,
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/gomp_params_contour.png",
  plot_contour,
  width = 6,
  height = 6,
  units = "in"
)

# Plot shape against covariates ------------------------------------------------

# Add shape and rate
f_df <- a_imp_df %>% inner_join(d_join_mdl)

# Continuous variables
# Age limits clearly banded
# Lower - 18 or >=18
# Upper - 80 or no upper limit
plot_shape_dbl <- f_df %>% 
  select(shape, where(is.double)) %>% 
  select(-c(n:attr, rate)) %>% 
  pivot_longer(duration:egfr_mean) %>% 
  ggplot(aes(x = value, y = shape)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(n.breaks = 10) +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "Shape")

# Binary
plot_shape_bin <- f_df %>% 
  select(shape, where(is.integer)) %>% 
  pivot_longer(trttype:one_pom) %>% 
  ggplot(aes(x = value, y = shape)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(n.breaks = 2) +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "Shape")

# Factors/ordinal
plot_shape_fct <- f_df %>% 
  select(shape, where(is.factor)) %>% 
  pivot_longer(bmi_lower:blinding) %>% 
  ggplot(aes(x = value, y = shape)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line", group = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "Shape")

# Save
ggsave(
  "outputs/model_shape/shape_vs_dbl.png",
  plot_shape_dbl,
  width = 14,
  height = 8,
  units = "in"
)

ggsave(
  "outputs/model_shape/shape_vs_bin.png",
  plot_shape_bin,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  "outputs/model_shape/shape_vs_fct.png",
  plot_shape_fct,
  width = 12,
  height = 4,
  units = "in"
)

# Design model -----------------------------------------------------------------

# Data
f_design_df <- e_join_chars %>% 
  select(-c(
    source, male_prop, bmi_lower:bmi_upper, contains(c("mean", "dbl"))
  ))

# Prepare model formula
z_covars_names <- names(f_design_df[c(10:26)])
z_covars <- paste0(z_covars_names, collapse = " + ")
z_fml <- as.formula(paste0("shape ~ rate + ", z_covars))

# Fit linear model
# Rate, maximum age (eligibility), maximium bmi (eligibility), phase and 
# quadruple blinding (compared to open-label) predict shape
f_mdl <- lm(z_fml, data = f_design_df)
summary(f_mdl)

# Diagnostics
f_diag <- broom::glance(f_mdl)

# Tidy outputs
f_res <- broom::tidy(f_mdl) %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  )

# Plot coefficients
plot_coefs <- f_res %>% 
  filter(!term == "(Intercept)") %>% 
  mutate(
    term = as.factor(term),
    term = fct_relevel(term, "rate")
  ) %>% 
  ggplot(aes(x = fct_rev(term), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 10) +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = 0, colour = "red") +
  labs(x = NULL, y = "Coefficient (95% CI)")

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
  theme_bw() +
  labs(x = "Model residuals", y = "Density")

# Predicted vs observed
plot_pred_obs <- f_resid %>% 
  ggplot(aes(x = shape_pred, y = shape)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "Predicted shape", y = "Observed shape")

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
# Including blinding significantly improves model fit (quadruple)
f_anova <- anova(f_mdl2, f_mdl)
summary(f_anova)

# Conduct model ----------------------------------------------------------------

# Data
g_conduct_df <- e_join_chars %>% 
  select(-c(source, bmi_lower:bmi_upper, contains("dbl")))

# Prepare model formula
z_covars2_names <- names(g_conduct_df[c(10:32)])
z_covars2 <- paste0(z_covars2_names, collapse = " + ")
z_fml2 <- as.formula(paste0("shape ~ rate + ", z_covars2))

# Fit linear model
# Same design predictors
# No baseline characteristic predicts shape
g_mdl <- lm(z_fml2, data = g_conduct_df)
summary(g_mdl)

# Diagnostics
g_diag <- broom::glance(g_mdl)

# Tidy outputs
g_res <- broom::tidy(g_mdl) %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  )

# Plot coefficients
plot_coefs2 <- g_res %>% 
  filter(!term == "(Intercept)") %>% 
  mutate(
    term = as.factor(term),
    term = fct_relevel(term, "rate")
  ) %>% 
  ggplot(aes(x = fct_rev(term), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 10) +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = 0, colour = "red") +
  labs(x = NULL, y = "Coefficient (95% CI)")

# Get predicted shape and residuals
g_resid <- broom::augment(g_mdl) %>% 
  as_tibble() %>% 
  rename(
    shape_pred = .fitted,
    resid = .resid
  ) %>% 
  select(shape, rate, shape_pred, resid)

# Plot residuals
plot_resid <- g_resid %>% 
  ggplot(aes(x = resid)) + 
  geom_density() + 
  theme_bw() +
  labs(x = "Model residuals", y = "Density")

# Predicted vs observed
plot_pred_obs <- g_resid %>% 
  ggplot(aes(x = shape_pred, y = shape)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "Predicted shape", y = "Observed shape")

# Residuals vs predicted
plot_pred_resid <- g_resid %>% 
  ggplot(aes(x = shape_pred, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "Predicted shape", y = "Residuals")

# Save results
write_csv(g_res, "outputs/model_shape/conduct/lm_res.csv")
write_csv(g_diag, "outputs/model_shape/conduct/lm_diags.csv")
write_csv(g_resid, "outputs/model_shape/conduct/lm_pred_resid.csv")

ggsave(
  "outputs/model_shape/conduct/01_lm_res_coefs.png",
  plot_coefs,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  "outputs/model_shape/conduct/02_lm_residuals.png",
  plot_resid,
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/conduct/03_lm_obs_pred_shape.png",
  plot_pred_obs,
  width = 6,
  height = 6,
  units = "in"
)

ggsave(
  "outputs/model_shape/conduct/04_lm_pred_resid.png",
  plot_pred_resid,
  width = 6,
  height = 6,
  units = "in"
)

# Fit model without blinding
z_covars2_anova <- gsub("\\+ blinding ", "", z_covars2)
z_fml2_anova <- as.formula(paste0("shape ~ rate + ", z_covars2_anova))
g_mdl2 <- lm(z_fml2_anova, data = g_conduct_df)
summary(g_mdl2)

# ANOVA test - significance
# Including blinding significantly improves model fit (quadruple)
g_anova <- anova(g_mdl2, g_mdl)
summary(g_anova)
