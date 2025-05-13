
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

plot_b_res <- b_res %>% 
  filter(!grepl("Intercept", term)) %>% 
  mutate(
    term = if_else(term == "trttype", "trt_arm", term),
    term = factor(
      term,
      levels = c("trt_arm", "n_250", "n_500", "n_750", "n_1000",
                 "arm_200", "arm_300", "arm_abv", 
                 "dur_6m", "dur_1yr", "arm_2", "arm_3",
                 "p3", "bdbl", "btrip", "bquad",
                 "f_100", "f_200", "pom_1", "pom_2", "pom_3",
                 "reg_am", "reg_eur", "reg_wpac", "reg_afr",
                 "reg_sea", "reg_em")
    )
  ) %>% 
  ggplot(aes(x = fct_rev(term), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_hline(yintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_hline(yintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  coord_flip() +
  scale_y_continuous(n.breaks = 20) +
  labs(y = "log-odds ratio", x = "covariate")

ggsave(
  "outputs/logreg_validation.png",
  plot_b_res,
  width = 8,
  height = 8,
  units = "in"
)

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