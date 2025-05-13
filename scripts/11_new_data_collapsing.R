
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
#' Sample size: 1 = <250, 2 = 250-499, 3 = 500-749, 4 = 750-999, 5 = >=1000
#' Arm size: 1 = <100, 2 = 100-199, 3 = 200-299, 4 = >=300
#' Duration: 1 = <6m, 2 = 6m - <1yr, 3 = >=1yr
#' Number of arms: 1 = 2, 2 = 3, 3 = 4, 4 = >4
#' Phase: 1 = 3, 2 = 4, 3 = Other
#' Blinding: 1 = Open-label, 2 = Double, 3 = Triple, 4 = Quadruple
#' Number of facilities: 1 = <100, 2 = 100-199, 3 = >=200
#' Primary outcomes: 1 = 1, 2 = 2, 3 = 3, 4 = >3

covs <- names(b_fmt[c(15:48)])
for (var in covs) {tidy_sum(b_fmt, var)}

c_clp <- b_fmt %>% 
  mutate(
    clp_n = case_when(
      n_250 == 1L ~ 1L,
      n_500 == 1L ~ 2L,
      n_750 == 1L ~ 3L,
      n_1000 == 1L ~ 4L,
      n_abv == 1L ~ 5L
    ),
    clp_arm_n = case_when(
      arm_100 == 1L ~ 1L,
      arm_200 == 1L ~ 2L,
      arm_300 == 1L ~ 3L,
      arm_abv == 1L ~ 4L
    ),
    clp_dur = case_when(
      dur_6m == 1L ~ 1L,
      dur_1yr == 1L ~ 2L,
      dur_abv == 1L ~ 3L
    ),
    clp_arms = case_when(
      arm_2 == 1L ~ 1L,
      arm_3 == 1L ~ 2L,
      arm_4 == 1L ~ 3L,
    ),
    clp_phase = case_when(
      p3 == 1L ~ 1L,
      p4 == 1L ~ 2L,
    ),
    clp_blind = case_when(
      bol == 1L ~ 1L,
      bdbl == 1L ~ 2L,
      btrip == 1L ~ 3L,
      bquad == 1L ~ 4L
    ),
    clp_fac = case_when(
      f_100 == 1L ~ 1L,
      f_200 == 1L ~ 2L,
      f_abv == 1L ~ 3L
    ),
    clp_pom = case_when(
      pom_1 == 1L ~ 1L,
      pom_2 == 1L ~ 2L,
      pom_3 == 1L ~ 3L,
      pom_abv == 1L ~ 4L
    )
  )

# Select variables
d_vars <- c_clp %>% 
  select(
    trial_id:duration, male_prop:hba1c, 
    contains(c("clp", "reg"))
  )

glimpse(d_vars)
d_vars %>% reframe(across(male_prop:reg_em, ~ sum(is.na(.))))

# Import original dataset for assigning NAs
# Address NAs
d_imp_orig <- read_csv("processed_data/df_stan_n348.csv")

d_fix <- d_vars %>% 
  left_join(
    d_imp_orig %>% 
      select(trial_id,  contains("phase"), n_arms) %>% 
      distinct()
  ) %>% 
  group_by(trial_id) %>% 
  fill(c(phase_3, phase_4, phase_notapp, n_arms), .direction = "downup") %>% 
  ungroup() %>% 
  mutate(
    clp_arm_n = case_when(
      is.na(clp_arm_n) & arm_n < 250 ~ 1L,
      is.na(clp_arm_n) & arm_n >= 250 & arm_n < 500 ~ 2L,
      is.na(clp_arm_n) & arm_n >= 500 & arm_n < 750 ~ 3L,
      is.na(clp_arm_n) & arm_n >= 750 & arm_n < 1000 ~ 4L,
      is.na(clp_arm_n) & arm_n >= 1000 ~ 5L,
      TRUE ~ clp_arm_n
    ),
    clp_arms = case_when(
      is.na(clp_arms) & n_arms == 2 ~ 1L,
      is.na(clp_arms) & n_arms == 3 ~ 2L,
      is.na(clp_arms) & n_arms == 4 ~ 3L,
      is.na(clp_arms) & n_arms > 4 ~ 4L,
      TRUE ~ clp_arms
    ),
    clp_phase = case_when(
      is.na(clp_phase) ~ 3L,
      TRUE ~ clp_phase
    )
  ) %>% 
  select(-c(phase_3:n_arms))

d_fix %>% reframe(across(male_prop:reg_em, ~ sum(is.na(.))))
glimpse(d_fix)

# Frequency check before saving
covs <- names(d_fix[c(7, 14:27)])
for (var in covs) {tidy_sum(d_fix, var)}

# Save as design and trial sub model sets
names(d_fix)
e_covs <- names(d_fix[c(14:27)]) # 14 covariates

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
f_covs <- names(d_fix[c(9:27)]) # 19 covariates

f_trial <- list(
  I = nrow(d_fix),
  K = length(d_covs),
  r = d_fix$arm_n_attr,
  N = d_fix$arm_n,
  dt = d_fix$duration,
  arm = as.integer(d_fix$trttype),
  X = as.matrix(d_fix[, d_covs])
)

saveRDS(f_trial, "mars/trial_sm_mk2/stan_dat_clp.rds")


