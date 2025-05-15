
rm(list = ls())
gc()

library(tidyverse)
library(brms)
library(bayesplot)
library(posterior)

# Fit a logreg model of the dataset to validate
c_imp_df <- readRDS("definitive_dataset.rds")
c_log_df <- c_imp_df %>% select(arm_n_compl:arm_n_attr, trttype:pom_abv)
c_rhs <- paste0(names(c_log_df[c(4:26)]), collapse = " + ")
c_formula <- as.formula(paste("arm_n_attr | trials(arm_n_attr + arm_n_compl) ~ trttype + ", c_rhs))

# Bayesian fit
c_fit <- brm(
  formula = c_formula,
  data = c_log_df,
  family = binomial(),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123
)

summary(c_fit)
saveRDS(c_fit, "brms_fit.rds")