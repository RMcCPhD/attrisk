
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)

# Function for tidy tables to check frequencies
tidy_sum <- function(df, var_name) {
  var <- df[[var_name]]
  cat("Variable:", var_name, "\n")
  if (is.integer(var)) {
    tbl <- table(var) %>% as.data.frame()
    print(tbl)
  } else if (is.numeric(var)) {
    tbl <- summary(var)
    print(tbl)
  } else {NA}
}

# Import
a_imp_df <- read_csv("processed_data/data_both_models.csv")
glimpse(a_imp_df)

# Format
# Treatment: 0 = Placebo, 1 = Active comparator
b_fmt <- a_imp_df %>% 
  mutate(
    trttype = case_match(
      trttype,
      0 ~ 1L,
      1 ~ 0L
    ),
    across(male:hba1c, ~ as.numeric(scale(.)) * 0.1),
    across(n_250:reg_em, ~ as.integer(.))
  )

glimpse(b_fmt)

# Collapse dummy groups into ordinal variables
#' Sample size: 0 = <250, 1 = 250-499, 2 = 500-749, 3 = 750-999, 4 = >=1000
#' Arm size: 0 = <100, 1 = 100-199, 2 = 200-299, 3 = >=300
#' Duration: 0 = <6m, 1 = 6m - <1yr, 2 = >=1yr
#' Number of arms: 0 = 2, 1 = 3, 2 = 4, 3 = >4
#' Phase: 0 = 3, 1 = 4, 2 = Other
#' Blinding: 0 = Open-label, 1 = Double, 2 = Triple, 3 = Quadruple
#' Number of facilities: 0 = <100, 1 = 100-199, 2 = >=200
#' Primary outcomes: 0 = 1, 1 = 2, 2 = 3, 3 = >3

covs <- names(b_fmt[c(15:48)])
for (var in covs) {tidy_sum(b_fmt, var)}

c_clp <- b_fmt %>% 
  mutate(
    clp_n = case_when(
      n_250 == 1L ~ 0L,
      n_500 == 1L ~ 1L,
      n_750 == 1L ~ 2L,
      n_1000 == 1L ~ 3L,
      n_abv == 1L ~ 4L
    ),
    clp_arm_n = case_when(
      arm_100 == 1L ~ 0L,
      arm_200 == 1L ~ 1L,
      arm_300 == 1L ~ 2L,
      arm_abv == 1L ~ 3L
    ),
    clp_dur = case_when(
      dur_6m == 1L ~ 0L,
      dur_1yr == 1L ~ 1L,
      dur_abv == 1L ~ 2L
    ),
    clp_arms = case_when(
      arm_2 == 1L ~ 0L,
      arm_3 == 1L ~ 1L,
      arm_4 == 1L ~ 2L,
    ),
    # clp_phase = case_when(
    #   p3 == 1L ~ 0L,
    #   p4 == 1L ~ 1L,
    # ),
    # clp_blind = case_when(
    #   bol == 1L ~ 0L,
    #   bdbl == 1L ~ 1L,
    #   btrip == 1L ~ 2L,
    #   bquad == 1L ~ 3L
    # ),
    clp_fac = case_when(
      f_100 == 1L ~ 0L,
      f_200 == 1L ~ 1L,
      f_abv == 1L ~ 2L
    )
    # clp_pom = case_when(
    #   pom_1 == 1L ~ 0L,
    #   pom_2 == 1L ~ 1L,
    #   pom_3 == 1L ~ 2L,
    #   pom_abv == 1L ~ 3L
    # )
  )

# Select variables
# Add additional dummy variables
d_vars <- c_clp %>% 
  select(
    trial_id:duration, male_prop:hba1c, 
    p3, p4, 
    bdbl, btrip, bquad, bol, pom_1, pom_2, pom_3, pom_abv,
    contains(c("clp", "reg"))
  ) %>% 
  mutate(
    po = if_else(p3 == 0L & p4 == 0L, 1L, 0L)
  ) %>% 
  select(trial_id:p4, po, everything())

glimpse(d_vars)
d_vars %>% reframe(across(male_prop:reg_em, ~ sum(is.na(.))))

# Import original dataset for assigning NAs
# Address NAs
d_imp_orig <- read_csv("processed_data/df_stan_n348.csv")

d_fix <- d_vars %>% 
  left_join(
    d_imp_orig %>% 
      select(trial_id, n_arms) %>% 
      distinct()
  ) %>% 
  group_by(trial_id) %>% 
  fill(n_arms, .direction = "downup") %>% 
  ungroup() %>% 
  mutate(
    clp_arm_n = case_when(
      is.na(clp_arm_n) & arm_n < 250 ~ 0L,
      is.na(clp_arm_n) & arm_n >= 250 & arm_n < 500 ~ 1L,
      is.na(clp_arm_n) & arm_n >= 500 & arm_n < 750 ~ 2L,
      is.na(clp_arm_n) & arm_n >= 750 & arm_n < 1000 ~ 3L,
      is.na(clp_arm_n) & arm_n >= 1000 ~ 4L,
      TRUE ~ clp_arm_n
    ),
    clp_arms = case_when(
      is.na(clp_arms) & n_arms == 2 ~ 0L,
      is.na(clp_arms) & n_arms == 3 ~ 1L,
      is.na(clp_arms) & n_arms == 4 ~ 2L,
      is.na(clp_arms) & n_arms > 4 ~ 3L,
      TRUE ~ clp_arms
    )
    # clp_phase = case_when(
    #   is.na(clp_phase) ~ 2L,
    #   TRUE ~ clp_phase
    # )
  ) %>% 
  # Removing covariates that add unnecessary complexity or have no
  # credible effect estimate
  select(-n_arms, -p4, -po)

d_fix %>% reframe(across(male_prop:reg_em, ~ sum(is.na(.))))
glimpse(d_fix)

# write_csv(d_fix, "processed_data/new_cov_formatting.csv")

# Frequency check before saving
covs <- names(d_fix[c(7, 14:33)])
for (var in covs) {tidy_sum(d_fix, var)}

# Save as design and trial sub model sets
names(d_fix)
e_covs <- names(d_fix[c(14:33)]) # 20 covariates

e_design <- list(
  I = nrow(d_fix),
  K = length(e_covs),
  r = d_fix$arm_n_attr,
  N = d_fix$arm_n,
  dt = d_fix$duration,
  arm = as.integer(d_fix$trttype),
  X = as.matrix(d_fix[, e_covs])
)

saveRDS(e_design, "mars/design_model_mk2/stan_dat_clp.rds")

names(d_fix)
f_covs <- names(d_fix[c(9:33)]) # 25 covariates

f_trial <- list(
  I = nrow(d_fix),
  K = length(f_covs),
  r = d_fix$arm_n_attr,
  N = d_fix$arm_n,
  dt = d_fix$duration,
  arm = as.integer(d_fix$trttype),
  X = as.matrix(d_fix[, f_covs])
)

saveRDS(f_trial, "mars/trial_sm_mk2/stan_dat_clp.rds")


