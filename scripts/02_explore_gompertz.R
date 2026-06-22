
# Exploratory analysis to determine the relationship between Gompertz
# parameters and trial characteristics among IPD trials
# Outputs of this will inform model specification, particularly the
# relationship between shape and rate with the addition of covariates
# Underlying mechanism is that Gompertz CDF acts as a link function supplying
# attrition probability parameter to binomial likelihood
# Parameters and covariates determine this so important to explore

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import tidied AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import intercept-only Gompertz models from attrition rates analyses
imp_fit <- read_csv("data/rates_tidy_res.csv")
imp_fit2 <- read_csv("data/rates_paper_res.csv")

# Prepare TTE model outputs
# Scale of estimates must match to the AgD data, i.e. in years
# This was the case for the remaining trial models, but not for the rates paper
# So rates model outputs are scaled to estimate TTE in years
prep_ipd_res <- imp_fit %>%
  select(-contains("se")) %>%
  rename(trial_id = nct_id) %>%
  rename_with(~gsub("_est", "", .)) %>%
  full_join(
    imp_fit2 %>%
      filter(dist == "gompertz") %>%
      select(-c(dist, se)) %>%
      pivot_wider(names_from = "parameter", values_from = "est") %>% 
      rename(trial_id = ctgov) %>%
      mutate(
        shape = shape * 365,
        rate = rate + log(365)
      )
  )

# 135 trials in total
n_distinct(prep_ipd_res$trial_id)

# Get TTE outputs for trials in AgD set
prep_ipd_agd <- prep_ipd_res %>% inner_join(imp_agd %>% distinct(trial_id))
n_distinct(prep_ipd_agd$trial_id)

# Get trial characteristics and nest join
prep_rdy <- prep_ipd_agd %>%
  left_join(
    imp_agd %>%
      group_by(trial_id) %>%
      nest(.key = "chars")
  ) %>% 
  rename(log_rate = rate) %>% 
  mutate(rate = exp(log_rate))

# Empirical data ---------------------------------------------------------------

# Parameter distributions
emp_distr <- prep_rdy %>%
  pivot_longer(c(shape, log_rate, rate)) %>%
  mutate(
    name = case_match(
      name, 
      "log_rate" ~ "Log-rate",
      "rate" ~ "Rate", 
      "shape" ~ "Shape"
    )
  ) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  theme_classic() +
  labs(x = "Parameter estimates", y = "Density")

# Parameter correlation
# Annotate plot with correlation test
emp_cor_tst <- cor.test(prep_rdy$log_rate, prep_rdy$shape, method = "pearson")
emp_cor_lbl <- paste0("r = ", round(emp_cor_tst$estimate, 3), "; p <0.01")
emp_cor_plot <- prep_rdy %>%
  ggplot(aes(x = log_rate, y = shape)) +
  geom_point(aes(x = log_rate, y = shape), alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) +
  theme_classic() +
  labs(x = "Log-rate", y = "Shape") +
  annotate(
    "text",
    label = emp_cor_lbl,
    x = -3.5, y = -3,
    fontface = "bold"
  )

# Combine plots
emp_plots <- cowplot::plot_grid(
  plotlist = list(emp_distr, emp_cor_plot),
  labels = c("A", "B")
)

ggsave(
  "outputs/design/gompertz_explore/01_empirical.png",
  emp_plots,
  width = 10,
  height = 4,
  units = "in"
)

# Parameters and covariates ----------------------------------------------------

# Explore the relationship between each parameter and covariates
# Informing design model specification
# - Covariates will be placed on rate by default, but exploring whether it is
#   informative to also place covariates on shape, if feasible

# Pivot parameters and characteristics
# Nest data by parameter for plotting
rel_dat <- prep_rdy %>% 
  unnest(chars) %>% 
  ungroup() %>%
  pivot_longer(c(shape, log_rate)) %>%
  mutate(
    name = case_match(
      name,
      "shape" ~ "Shape",
      "log_rate" ~ "Log-rate"
    )
  )

# Binary characteristics
rel_binary <- rel_dat %>%
  select(name, value, where(is.integer)) %>%
  pivot_longer(
    -c(name, value),
    names_to = "covar",
    values_to = "covar_value"
  ) %>%
  group_by(name) %>%
  nest() %>%
  arrange(name) %>%
  mutate(
    plot = map2(data, name, ~ {
      .x %>%
        ggplot(aes(x = covar_value, y = value)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE) +
        facet_wrap(~covar) +
        theme_classic() +
        scale_x_continuous(n.breaks = 2) +
        labs(x = "Binary covariates", y = .y)
    })
  )

rel_binary_plots <- cowplot::plot_grid(
  plotlist = rel_binary$plot,
  labels = c("A", "B")
)

ggsave(
  "outputs/design/gompertz_explore/02a_relationship_parameters_binary.png",
  rel_binary_plots,
  width = 12,
  height = 4,
  units = "in"
)

# Factors
rel_factor <- rel_dat %>%
  select(name, value, where(is.factor)) %>%
  pivot_longer(
    -c(name, value),
    names_to = "covar",
    values_to = "covar_value"
  ) %>%
  group_by(name) %>%
  nest() %>%
  arrange(name) %>%
  mutate(
    plot = map2(data, name, ~ {
      .x %>%
        ggplot(aes(x = covar_value, y = value)) +
        geom_jitter(width = 0.15, alpha = 0.1) +
        stat_summary(
          fun = mean,
          geom = "point",
          size = 3,
          shape = 18,
          # colour = "orange2"
        ) +
        stat_summary(
          aes(group = 1),
          fun = mean,
          geom = "line",
          # colour = "orange2"
        ) +
        facet_wrap(~covar, scales = "free") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Categorical covariates", y = .y)
    })
  )

rel_factor_plots <- cowplot::plot_grid(
  plotlist = rel_factor$plot,
  labels = c("A", "B")
)

ggsave(
  "outputs/design/gompertz_explore/02b_relationship_parameters_factor.png",
  rel_factor_plots,
  width = 14,
  height = 8,
  units = "in"
)

# Model parameters and covariates ----------------------------------------------

# Exploration suggests no or weak relationship between shape and covariates
# Covariates will be placed on log_rate only
# Parsimonious structure: log-rate explained by covariates, shape partially
# explained by log_rate with some degree of uncertainty
# Additionally avoids double-adjustment (i.e. covariates predict log-rate,
# log-rate predicts shape; same covariates involved twice)

# Fit model of shape on rate
mdl_shape <- lm(shape ~ log_rate, data = prep_rdy)
mdl_shape_sum <- tidy(mdl_shape)
mdl_shape_diag <- glance(mdl_shape)

# Inspect
# - Plot predicted against observed shape
# - Plot residuals distribution
mdl_shape_dat <- broom::augment(mdl_shape) %>% as_tibble()

mdl_shape_plot <- mdl_shape_dat %>%
  ggplot(aes(x = .fitted, y = shape)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Predicted shape", y = "Shape")

mdl_shape_plot2 <- mdl_shape_dat %>%
  ggplot(aes(x = .resid)) +
  geom_density() +
  theme_classic() +
  labs(x = "Model residuals", y = "Density")

mdl_shape_plots <- cowplot::plot_grid(
  plotlist = list(mdl_shape_plot, mdl_shape_plot2),
  labels = c("A", "B")
)

ggsave(
  "outputs/design/gompertz_explore/03_model_shape_rate.png",
  mdl_shape_plots,
  width = 8,
  height = 4,
  units = "in"
)

# Get residual standard deviation
# Acts as an uncertainty parameter (Uncertainty ~ Normal(0, RSD)) when
# specifying relationship between shape and rate in model specification
# i.e. Rate predicts shape to some degree of uncertainty draw from a normal
#      distribution (since both parameters are roughly regular)
mdl_shape_rsd <- sigma(mdl_shape)

# Extract outputs for model specification
mdl_out <- tibble(
  shape_alpha = mdl_shape_sum$estimate[[1]],
  lograte_beta = mdl_shape_sum$estimate[[2]],
  uncertainty = mdl_shape_rsd
)

# Save
saveRDS(mdl_out, "processed_data/shape_rate.rds")

# Save model outputs
write_csv(mdl_shape_sum, "outputs/design/gompertz_explore/03a_mdl_out.csv")
write_csv(mdl_shape_diag, "outputs/design/gompertz_explore/03b_mdl_diag.csv")
