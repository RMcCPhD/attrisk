
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)

# Import
a_imp_df <- read_csv("processed_data/data_both_models.csv")
glimpse(a_imp_df)

# Format, incl treatment type swapped (0 - placebo, 1 - active)
a_imp_fmt <- a_imp_df %>% 
  mutate(
    trttype = if_else(trttype == 1, 0, 1),
    across(male:hba1c, ~ as.numeric(scale(.)) * 0.1),
    across(n_250:reg_em, ~ as.integer(.))
  )

# Split into datasets
b_trialsub <- a_imp_df

# Design-stage -----------------------------------------------------------------

# Remove perfectly collinear variables
b_design <- a_imp_df %>%
  select(-c(male:hba1c, n_abv, arm_100, dur_abv, arm_4, p4, bol, f_abv, pom_abv))

# Vector for formula and covariates (25 covariates)
b_rhs <- paste0(names(b_design[c(9:34)]), collapse = " + ")
b_formula <- as.formula(paste("cbind(arm_n_attr, arm_n_compl) ~ trttype + ", b_rhs))

# Include adjustment by treatment type
b_fit <- glm(b_formula, data = b_design, family = "binomial")
summary(b_fit)

# Plot coefficients and limits
b_res <- broom::tidy(b_fit) %>% 
  mutate(
    p.value = round(p.value, 3),
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  )

b_res %>% 
  filter(!grepl("Intercept", term) & p.value < 0.05) %>% 
  ggplot(aes(x = fct_rev(term), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  theme_bw() +
  coord_flip()

# Nest by treatment type, fit models
b_nest <- b_design %>% 
  group_by(trttype) %>% 
  nest() %>% 
  mutate(
    mdl = map(data, ~ glm(b_formula, data = .x, family = "binomial"))
  )

# Get outputs
summary(b_nest$mdl[[1]])
summary(b_nest$mdl[[2]])

c_outputs <- b_nest %>% 
  mutate(
    res = map(mdl, ~ broom::tidy(.x)),
    diag = map(mdl, ~ broom::glance(.x)),
    vcov = map(mdl, ~ vcov(.x))
  )

c_res <- c_outputs %>% 
  select(trttype, res) %>% 
  unnest(res) %>% 
  mutate(
    p.value = round(p.value, 3),
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  )

# Plot
c_plot_res <- c_res %>% 
  filter(!grepl("Intercept", term) & p.value < 0.05) %>% 
  mutate(trttype = as.factor(trttype)) %>% 
  ggplot(aes(x = trttype, y = estimate, ymin = lci, ymax = uci)) +
  geom_point() +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~ term, scales = "free") +
  theme_bw() +
  coord_flip()

c_plot_res


# Trial submodel ---------------------------------------------------------------