
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)

# Import data before converting numeric covariates to ordinal/factors
a_imp_df <- read_csv("processed_data/df_stan_n348.csv")
glimpse(a_imp_df)

# Select covariates
b_dbl <- a_imp_df %>% 
  select(trial_id, n_arms, n_facilities, n_poms:n_regions) %>% 
  distinct()

# Import data with ordinal/factor coding
b_imp_df <- read_csv("processed_data/new_cov_formatting.csv")

# Replace ordinal/factors with numeric versions
# Replace regions with a count of regions involved
# Add in original participant characteristics
c_swap <- b_imp_df %>% 
  left_join(b_dbl) %>% 
  select(-c(contains(c("clp", "pom_"))))

glimpse(c_swap)
c_swap %>% reframe(across(everything(), ~ sum(is.na(.))))

# Import original participant characteristics
d_part_chars <- read_csv("processed_data/data_both_models.csv") %>% 
  select(trial_id, n, arm_n, trttype, male_prop:hba1c) %>% 
  mutate(trttype = case_match(trttype, 0 ~ 1, 1 ~ 0))

# Join numeric characteristics
# Remove BMI as it's bimodal, reduce complexity
# Format
d_add_chars <- c_swap %>% 
  select(-c(male_prop:hba1c)) %>% 
  left_join(d_part_chars) %>% 
  select(trial_id:arm_n_attr, male_prop:hba1c, everything()) %>% 
  mutate(
    trttype = as.integer(trttype),
    across(p3:bol, ~ as.integer(.))
  )

glimpse(d_add_chars)

# Summarise numeric variables to determine scaling
d_add_chars %>% 
  reframe(
    across(
      c(age:hba1c, duration, n_arms:n_regions),
      ~ unlist(
        list(
          mean = mean(.), sd = sd(.), 
          median = median(.), min = min(.), max = max(.)
        )
      )
    )
  )

# Male_prop, age, t2dur, hba1c are roughly normal, scale to mean 0, sd 1
# names(d_add_chars[c(7:11, 13, 22)])
# for (i in d_add_chars[c(7:11, 13, 22)]) {print(hist(i))}
# for (i in d_add_chars[c(7:11, 13, 22)]) {print(hist(as.numeric(scale(i))))}
# for (i in d_add_chars[c(7:11, 13, 22)]) {print(hist(as.numeric(scale(log(i)))))}

# BMI is bimodal, remove to reduce complexity, scaling/transform not suitable
# hist(as.numeric(scale(log(d_add_chars$bmi))))

# Duration and n_facilities are skewed, log transform and scale the same
# for (i in d_add_chars[c(13, 22)]) {print(hist(as.numeric(scale(log(i)))))}

# Sample size and arm size can be log-transformed and scaled the same
# for (i in d_add_chars[c(3, 4)]) {print(hist(as.numeric(scale(log(i)))))}

# N_arms, n_poms, n_regions as either ordinal or factors
# names(d_add_chars[c(21, 23:24)])
# for (i in d_add_chars[c(21, 23:24)]) {print(hist(i))}

# Format variables according to evaluation
e_format <- d_add_chars %>% 
  select(-bmi) %>% 
  mutate(
    sample_size = n,
    arm_size = arm_n,
    dur_yr = duration,
    scale_male_prop = as.numeric(scale(male_prop)),
    across(
      c(age, t2dur, hba1c), 
      ~ as.numeric(scale(.) * 0.1), 
      .names = "scale_{.col}"
    ),
    across(
      c(sample_size, arm_n, dur_yr, n_facilities), 
      ~ as.numeric(scale(log(.))),
      .names = "scale_{.col}"
    ),
    across(contains("reg_"), ~ as.integer(.)),
    arms_2 = if_else(n_arms == 2, 1L, 0L),
    arms_3 = if_else(n_arms == 3, 1L, 0L),
    arms_abv = if_else(n_arms > 3, 1L, 0L),
    pom_1 = if_else(n_poms == 1, 1L, 0L),
    pom_2 = if_else(n_poms == 2, 1L, 0L),
    pom_abv = if_else(n_poms > 2, 1L, 0L)
  ) %>% 
  select(-c(male_prop:hba1c, n_arms:n_regions, sample_size:dur_yr)) %>%
  select(trial_id:arm_n_attr, duration, scale_male_prop:scale_hba1c, everything())

glimpse(e_format)
e_format %>% reframe(across(everything(), ~ sum(is.na(.))))

# saveRDS(e_format, "processed_data/definitive_dataset.rds")

e_format %>% 
  # select(male_prop:hba1c) %>% 
  # mutate(across(everything(), ~ scale(log(.)))) %>% View()
  select(scale_male_prop:scale_hba1c) %>%
  pivot_longer(everything(), names_to = "var", values_to = "num") %>% 
  ggplot(aes(x = num)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free")

# Design data
names(e_format)
f_covs <- names(e_format[c(12:35)])

# 24 variables, 10 predictor groups
f_design <- list(
  I = nrow(e_format),
  K = length(f_covs),
  r = e_format$arm_n_attr,
  N = e_format$arm_n,
  dt = e_format$duration,
  X = as.matrix(e_format[, f_covs])
)

# saveRDS(f_design, "mars/design_model_mk4/stan_dat.rds")

names(e_format)
g_covs <- names(e_format[c(8:10, 12:35)]) # 28 covariates, 13 groups

# 28 variables, 14 predictor groups
g_trial <- list(
  I = nrow(e_format),
  K = length(g_covs),
  r = e_format$arm_n_attr,
  N = e_format$arm_n,
  dt = e_format$duration,
  arm = as.integer(e_format$trttype),
  X = as.matrix(e_format[, g_covs])
)

saveRDS(g_trial, "mars/trial_sm_mk3/stan_dat.rds")

