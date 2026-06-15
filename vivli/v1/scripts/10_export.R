
# Prepare model results and variance-covariance for export
source("Scripts/00_config.R")
source("Scripts/00_functions.R")
library(tidyverse)

# Import results and vcov
a_imp_res <- readRDS("processed_data/cox_outputs/cox_outputs.rds")
a_imp_vcov <- readRDS("processed_data/cox_outputs/cox_vcov.rds")

# Check sponsor frequency
a_imp_lkp <- read_csv("Created_metadata/lookup_trials.csv")

a_imp_lkp %>%
  rename(nct_id = ctgov) %>% 
  select(sponsor, nct_id) %>% 
  inner_join(a_imp_res %>% select(nct_id)) %>% 
  count(sponsor, nct_id) %>% 
  count(sponsor) %>% 
  mutate(pcnt = n/74 * 100)

# Prepare model results
# Remove NA estimates
b_res <- a_imp_res %>% 
  select(nct_id, res) %>% 
  unnest(res) %>% 
  mutate(across(estimate:p.value, ~ round(., 3))) %>% 
  filter(!is.na(estimate))

# Prepare diagnostics
b_diag <- a_imp_res %>% 
  select(nct_id, diag) %>% 
  unnest(diag) %>% 
  mutate(across(n:nobs, ~ round(., 3)))

# Create covariate correlation data from vcov matrices for export
# Keep upper triangle (removing lower and diagonals), rest as NA
c_vcov <- a_imp_vcov %>% 
  select(nct_id, vcov) %>% 
  mutate(
    vcov = map(vcov, function(v) {
      cor <- cov2cor(v)
      cor[!upper.tri(cor)] <- NA_real_
      cor
    })
  )

# Convert to long-format data
# Long-format and filter out NAs (lower triangle, diagonals)
c_tidy_vcov <- c_vcov %>% 
  mutate(
    vcov_long = map(vcov, function(v) {
      df <- as.data.frame(v) %>% 
        as_tibble(rownames = "row") %>% 
        gather("col", "value", -row) %>% 
        filter(!is.na(value))
      df
    })
  ) %>% 
  select(nct_id, vcov_long) %>% 
  unnest(vcov_long) %>% 
  ungroup()

# Save outputs
write_csv(b_res, "Export/01_res.csv")
write_csv(b_diag, "Export/02_diag.csv")
write_csv(c_tidy_vcov, "Export/03_cor.csv")
