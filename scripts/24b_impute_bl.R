
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import bl chars
a_imp_bl <- readRDS("processed_data/bl_chars_noimp.rds")

# Check missingness
a_imp_bl %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  print(n = 25)

# Setup for amelia -------------------------------------------------------------

# Format trttype, remove refs for proportions
b_format <- a_imp_bl %>% 
  mutate(trttype = factor(trttype)) %>% 
  select(-c(white, former_never))

# Set bounds for proportion imputations
prop_vars <- intersect(
  c("female_prop", "asian", "black", "other", "current"), 
  names(b_format)
)

bounds <- do.call(
  rbind, 
  lapply(prop_vars, function(v) c(match(v, names(b_format)), 0, 1))
)

# Set bounds for mean and sd imputations
b_format %>% 
  reframe(across(contains(c("mean", "sd")), ~ unlist(
    list(min = min(., na.rm = TRUE), max = max(., na.rm = TRUE))
  ))) %>% 
  pivot_longer(everything()) %>% 
  print(n = 32)

marg_limits <- rbind(
  c(match("mean_age", names(b_format)), 39.5, 74.9),
  c(match("sd_age", names(b_format)), 3.8, 12.8),
  c(match("mean_sbp", names(b_format)), 124, 151),
  c(match("sd_sbp", names(b_format)), 6.7, 17.4),
  c(match("mean_bmi", names(b_format)), 24.2, 37.4),
  c(match("sd_bmi", names(b_format)), 2.6, 10.6),
  c(match("mean_t2dur", names(b_format)), 0.7, 20.7),
  c(match("sd_t2dur", names(b_format)), 1.4, 10.1),
  c(match("mean_egfr", names(b_format)), 24.1, 154),
  c(match("sd_egfr", names(b_format)), 0.25, 70),
  c(match("mean_hba1c", names(b_format)), 6.7, 10.2),
  c(match("sd_hba1c", names(b_format)), 0.35, 1.8),
  c(match("mean_hdl", names(b_format)), 0.94, 2.7),
  c(match("sd_hdl", names(b_format)), 0.02, 0.07),
  c(match("mean_total_chol", names(b_format)), 4, 9.8),
  c(match("sd_total_chol", names(b_format)), 0.07, 1.9)
)

bounds <- rbind(
  bounds, 
  marg_limits[!is.na(marg_limits[,1]), , drop = FALSE]
)

# Run amelia -------------------------------------------------------------------

set.seed(123)

b_impute <- amelia(
  b_format,
  m = 20,
  idvars  = "trial_id",
  noms    = "trttype",
  empri   = 0.01 * nrow(b_format),
  bounds  = bounds,
  incheck = TRUE
)

# Fix race and smoking groups imputed, sum to 1
fix_grps <- function(df) {
  
  # Smoking groups
  df$current <- pmin(pmax(df$current, 0), 1)
  df$former_never <- 1 - df$current
  
  # Race groups
  race_cols <- c("asian", "black", "other")
  race_X <- as.matrix(df[race_cols])
  race_X[!is.finite(race_X)] <- NA
  race_X[is.na(race_X)] <- 0
  race_X[race_X < 0] <- 0
  race_X[race_X > 1] <- 1
  
  race_s <- rowSums(race_X)
  zero <- !is.finite(race_s) | race_s <= 0
  if (any(zero)) race_X[zero, ] <- 1/3
  
  race_s <- rowSums(race_X)
  race_big <- race_s > 1
  if (any(race_big)) race_X[race_big, ] <- race_X[race_big, ] / race_s[race_big]
  
  df[race_cols] <- race_X
  df$white <- 1 - rowSums(race_X)
  df$white <- pmin(pmax(df$white, 0), 1)
  df %>%
    select(trial_id:other, white, everything()) %>% 
    mutate(
      across(mean_age:former_never, ~ round(., 2)),
      race_total = asian + black + other + white,
      white = case_when(
        race_total == 0.99 ~ white + 0.01,
        race_total == 1.01 ~ white - 0.01,
        TRUE ~ white
      ),
      race_total = asian + black + other + white,
      smoke_total = current + former_never
    )
}

b_fixed <- lapply(b_impute$imputations, fix_grps)

# Add imputation iteration to datasets, bind rows
c_bind <- bind_rows(
  lapply(seq_along(b_fixed), \(j) mutate(b_fixed[[j]], imp = j))
)

# Get marginals across imputations ---------------------------------------------

# Check distributions of variable marginals
c_bind %>% 
  pivot_longer(mean_age:former_never) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value, group = imp)) +
  geom_density(alpha = 0.01) +
  facet_wrap(~name, scales = "free")






