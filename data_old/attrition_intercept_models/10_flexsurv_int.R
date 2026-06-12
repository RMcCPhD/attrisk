
source("Scripts/00_config.R")
library(tidyverse)
library(flexsurv)

# Import dataset
# Remove NCT02025907 consistent with cox set (no variables converged)
# 71 trials
a_imp_df <- readRDS("Processed_data/n72_cox.rds")

# Nest
# Fit flexsurv intercept-only models
# Get tidy outputs
b_nest <- a_imp_df %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  mutate(
    flexsurv = map(
      data, 
      ~ flexsurvreg(Surv(weeks/52, event) ~ 1, dist = "gompertz", data = .x)
    ),
    tidy = map(
      flexsurv,
      ~ t(.x$res.t) %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "outputs") %>%
        filter(outputs %in% c("est", "se"))
    ),
    vcov = map(flexsurv, ~ vcov(.x))
  )

# Extract tidy outputs
c_tidy <- b_nest %>% 
  select(nct_id, tidy) %>% 
  unnest(cols = c("tidy")) %>% 
  pivot_wider(names_from = "outputs", values_from = c("shape", "rate")) %>% 
  mutate(across(shape_est:rate_se, ~ round(., 3)))

# Save
write_csv(c_tidy, "Outputs/flexsurv_int/tidy_res.csv")

# Extract vcov
c_vcov <- b_nest %>% select(nct_id, vcov)

# Convert vcov to correlation matrices
# Keep upper triangle (removing lower and diagonals), rest as NA
c_vcov <- c_vcov %>% 
  mutate(
    cor_mat = map(vcov, function(v) {
      cor <- cov2cor(v)
      cor[!upper.tri(cor)] <- NA_real_
      cor
    })
  )

# Convert to long-format data
# Long-format and filter out NAs (lower triangle, diagonals)
c_cor <- c_vcov %>% 
  mutate(
    cor_long = map(vcov, function(v) {
      df <- as.data.frame(v) %>% 
        as_tibble(rownames = "row") %>% 
        gather("col", "value", -row) %>% 
        filter(!is.na(value))
      df
    })
  ) %>% 
  select(nct_id, cor_long) %>% 
  unnest(cor_long) %>% 
  ungroup()

c_cor %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save
write_csv(c_cor, "Outputs/flexsurv_int/tidy_cov2cor.csv")

# Review -----------------------------------------------------------------------

# Summarise
c_tidy %>%
  ungroup() %>% 
  reframe(
    mean_shape = mean(shape_est),
    sd_shape = sd(shape_est),
    mean_rate = mean(rate_est),
    sd_rate = sd(rate_est)
  )

c_tidy %>% 
  ungroup() %>% 
  select(contains("est")) %>% 
  pivot_longer(contains("est"), names_to = "param", values_to = "est") %>% 
  ggplot(aes(x = est)) +
  geom_density() +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

# Correlation (-0.34)
cor(c_tidy$shape_est, c_tidy$rate_est)

# Fit linear models
d_fit1 <- lm(shape_est ~ rate_est, data = c_tidy)
summary(d_fit1)

d_fit2 <- lm(rate_est ~ shape_est, data = c_tidy)
summary(d_fit2)

