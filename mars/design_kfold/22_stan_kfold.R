
# Prepare design data for K-fold cross-validation
# Trial-level partitions to conserve hierarchical structure
# Split into 10 folds, sets of 34 or 35 (even split not possible)
rm(list = ls())
gc()

library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Import data and model coefficients
a_imp_df <- readRDS("latest_dataset.rds")
a_imp_mdl_design_nm <- read_csv("lm_res.csv")
a_imp_mdl_design <- readRDS("lm_res_shrinkage.rds")

# Get variables
b_design_df <- a_imp_df %>% 
  mutate(attr_prop = attr/arm_n) %>% 
  select(-c(contains("dbl"), bmi_lower:bmi_upper, male_prop:egfr_mean))

# Create model matrix for predictors
b_vars <- names(b_design_df[9:25])
b_fml <- as.formula(paste0("~ ", paste0(b_vars, collapse = " + ")))

# Prepare design coefficients
b_design_coefs <- a_imp_mdl_design_nm %>% 
  filter(!(term %in% c("(Intercept)", "rate"))) %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

b_coefs_names <- names(b_design_coefs)
b_coefs_fml <- as.formula(paste0("~ ", paste0(b_coefs_names, collapse = " + ")))

b_coefs_beta <- tibble(term = b_coefs_names, estimate = a_imp_mdl_design) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

b_coefs_mat <- model.matrix(b_coefs_fml, data = b_coefs_beta)[-1]

# Assign folds to trials
set.seed(123)
b_folds <- b_design_df %>% 
  select(trial_id) %>% 
  distinct() %>% 
  cbind(tibble(folds = sample(rep(1:10, length.out = 344))))

# Join folds to dataset
b_df_folds <- b_design_df %>% left_join(b_folds) %>% select(folds, everything())

# Config
options(scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile stan model
b_compile <- stan_model("mdl_design.stan")

# Loop over folds
c_cv_res <- vector("list", 10)

for (k in 1:10) {
  cat("Fold:", k, " | Time:", Sys.time(), "\n")
  
  # Define training and testing data
  train <- b_df_folds %>% filter(folds != k)
  test <- b_df_folds %>% filter(folds == k)
  
  # Design matrices
  train_mat <- model.matrix(b_fml, data = train)[, -1]
  test_mat <- model.matrix(b_fml, data = test)[, -1]
  
  # Shape betas matrix
  shape_mat <- b_coefs_mat
  
  # Stan data
  stan_data <- list(
    I = nrow(train),
    K = ncol(train_mat),
    r = train$attr,
    N = train$arm_n,
    dt = train$duration,
    B = shape_mat,
    X = train_mat
  )
  
  # Fit model
  fit <- sampling(
    b_compile,
    data = stan_data,
    warmup = 1000,
    iter = 4000,
    chains = 4,
    seed = 123,
    control = list(adapt_delta = 0.9999, max_treedepth = 15)
  )
  
  # Extract posterior means
  mu1_post <- rstan::extract(fit, pars = "mu1")$mu1 %>% colMeans()
  beta_post <- rstan::extract(fit, pars = "beta")$beta %>% colMeans()
  
  # Predict b, a, and p for test set
  b_pred <- exp(mu1_post[1] + test_mat %*% beta_post)
  a_pred <- b_pred + test_mat %*% shape_mat
  
  dt_test <- test$duration
  N_test <- test$arm_n
  r_test <- test$attr
  
  # Gompertz predicted probabilities
  p_pred <- ifelse(abs(a_pred) < 1e-4,
                   1 - exp(-b_pred * dt_test),
                   1 - exp(-b_pred / a_pred * (exp(a_pred * dt_test) - 1)))
  
  # Compute predictive quantities
  r_pred <- N_test * p_pred
  log_lik <- dbinom(r_test, size = N_test, prob = p_pred, log = TRUE)
  
  c_cv_res[[k]] <- tibble(
    fold = k,
    trial_id = test$trial_id,
    r_obs = r_test,
    r_pred = r_pred,
    log_lik = log_lik
  )
}

# Bind lists
c_cv_df <- bind_rows(c_cv_res) %>% 
  mutate(
    residual = r_obs - r_pred,
    abs_residual = abs(residual)
  )

# Summarise
c_cv_sum <- c_cv_df %>% 
  summarise(
    RMSE = sqrt(mean((r_obs - r_pred)^2)),
    mean_log_score = mean(log_lik),
    R2 = 1 - sum((r_obs - r_pred)^2) / sum((r_obs - mean(r_obs))^2)
  )

# Summarise by fold
c_cv_sum_fold <- c_cv_df %>% 
  group_by(fold) %>% 
  summarise(
    RMSE = sqrt(mean((r_obs - r_pred)^2)),
    mean_log_score = mean(log_lik),
    R2 = 1 - sum((r_obs - r_pred)^2) / sum((r_obs - mean(r_obs))^2)
  )


# Save
saveRDS(c_cv_res, "outputs/cv_lst.rds")
saveRDS(c_cv_df, "outputs/cv_df.rds")
write_csv(c_cv_sum, "outputs/cv_sum.csv")
write_csv(c_cv_sum_fold, "outputs/cv_sum_folds.csv")
