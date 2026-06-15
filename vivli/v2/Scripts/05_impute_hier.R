
# Fit hierarchical models
# Predict values for completely missing variables within trials
# Save covariates with full imputation implemented

# Workflow
# Fit hierarchical model for each covariate that is completely missing
# Remove trial(s) missing the covariate from specification
# Fit model and get outputs
# Predict covariates for that trial

# Random intercepts model - random intercept for trial
# Otherwise takes extremely long time to run, hadn't seen one fit in runtime

source("Scripts/00_config.R")
source("Scripts/00_functions.R")

set.seed(123)

library(tidyverse)
library(broom)
library(lme4)

# Import data
a_imp_df <- readRDS("Processed_data/04_within_trial_imputed.rds")

# Race as factor, with "missing_unk" for NAs
# Log-transform cholesterol for predictions
b_race <- a_imp_df %>% 
  mutate(
    race = if_else(is.na(race), "na_unk", race),
    race = factor(race, levels = c("white", "asian", "black_aa", "other", "na_unk")),
    across(contains("chol"), ~ log(.))
  )

# Check missingness
b_race %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Identify completely missing variables
b_na <- b_race %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(n = length(unique(id))) %>% 
  left_join(
    b_race %>% 
      group_by(sponsor, nct_id) %>% 
      reframe(across(-c(id:mft), ~ sum(is.na(.))
      ))
  ) %>% 
  pivot_longer(-c(sponsor, nct_id, n)) %>% 
  group_by(sponsor, nct_id, name) %>% 
  mutate(all_na = value == n) %>% 
  ungroup()

# 4 trials
# 2 missing BMI
# 2 missing cholesterol (all three - hdl, ldl, total)
# 1 missing blood pressure (both - dbp and sbp)
b_na %>% filter(all_na) %>% count(nct_id)
b_na %>% filter(all_na) %>% count(name)
b_na %>% filter(all_na) %>% filter(grepl("bmi", name))
b_na %>% filter(all_na) %>% filter(grepl("chol", name))
b_na %>% filter(all_na) %>% filter(name %in% c("dbp", "sbp"))

# BMI --------------------------------------------------------------------------

# 74/76 trials
# Add 0 for any other variables that are completely missing
c_bmi_df <- b_race %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("bmi", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains("chol"), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  )

# Fit mixed-effects model
c_bmi_vars <- c(
  "age10", "sex", "race", "durn", "smoke", "dbp", "sbp", 
  "chol_hdl", "chol_ldl", "chol_total", 
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_bmi_p <- paste0(names(c_bmi_df[c_bmi_vars]), collapse = " + ")
c_bmi_lp <- paste0("bmi ~ ", c_bmi_p, " + (1 | nct_id)")
c_bmi_fml <- as.formula(c_bmi_lp)
c_bmi_mdl <- lmer(formula = c_bmi_fml, data = c_bmi_df)
summary(c_bmi_mdl)

# Predict BMI values using predict()
c_bmi_newdat <- b_race %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("bmi", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains("chol"), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  )

c_bmi_pred <- predict(
  c_bmi_mdl, 
  newdata = c_bmi_newdat, 
  allow.new.levels = TRUE
)

# Replace NA bmi with predictions
c_bmi_add <- c_bmi_newdat %>% 
  select(-bmi) %>% 
  cbind(tibble(bmi = c_bmi_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi,
    everything()
  ) %>% 
  as_tibble()

# Incorporate imputed BMI
c_bmi_inc <- b_race %>% 
  left_join(c_bmi_add %>% select(nct_id, id, bmi_new = bmi)) %>% 
  mutate(bmi = coalesce(bmi, bmi_new)) %>% 
  select(-bmi_new)

c_bmi_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after (saved in full at the end)
b_race %>% 
  select(bmi) %>% 
  mutate(type = "before") %>% 
  full_join(c_bmi_inc %>% select(bmi) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = bmi, colour = type)) +
  geom_density() +
  theme_classic()

# HDL cholesterol --------------------------------------------------------------

# 74/76 trials
# Add 0 for any other variables that are completely missing
# Omit LDL from prediction due to collinearity
# Predict log-HDL then back-transform
c_hdl_df <- c_bmi_inc %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("chol_hdl", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains(c("ldl", "total")), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-chol_total)

# Fit mixed-effects model
c_hdl_vars <- c(
  "age10", "sex", "race", "bmi", "durn", "smoke", "dbp", "sbp", 
  "chol_ldl", 
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_hdl_p <- paste0(names(c_hdl_df[c_hdl_vars]), collapse = " + ")
c_hdl_lp <- paste0("chol_hdl ~ ", c_hdl_p, " + (1 | nct_id)")
c_hdl_fml <- as.formula(c_hdl_lp)
c_hdl_mdl <- lmer(formula = c_hdl_fml, data = c_hdl_df)
summary(c_hdl_mdl)

# Predict chol_hdl values
c_hdl_newdat <- c_bmi_inc %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("chol_hdl", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains(c("chol")), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-chol_hdl)

c_hdl_pred <- predict(
  c_hdl_mdl, 
  newdata = c_hdl_newdat, 
  allow.new.levels = TRUE
)

# Replace NA hdl with predictions
c_hdl_add <- c_hdl_newdat %>% 
  cbind(tibble(chol_hdl = c_hdl_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi, durn, smoke, dbp, sbp,
    chol_hdl,
    everything()
  ) %>%
  as_tibble()

# Incorporate imputed BMI
c_hdl_inc <- c_bmi_inc %>% 
  left_join(c_hdl_add %>% select(nct_id, id, chol_hdl_new = chol_hdl)) %>% 
  mutate(chol_hdl = coalesce(chol_hdl, chol_hdl_new)) %>% 
  select(-chol_hdl_new)

c_hdl_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after
c_bmi_inc %>% 
  select(chol_hdl) %>% 
  mutate(type = "before") %>% 
  full_join(c_hdl_inc %>% select(chol_hdl) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = chol_hdl, colour = type)) +
  geom_density() +
  theme_classic()

# LDL cholesterol --------------------------------------------------------------

# 74/76 trials
# Add 0 for any other variables that are completely missing
# Remove HDL cholesterol due to collinearity
c_ldl_df <- c_hdl_inc %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("chol_ldl", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(chol_total, dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-chol_total)

# Fit mixed-effects model
c_ldl_vars <- c(
  "age10", "sex", "race", "bmi", "durn", "smoke", "dbp", "sbp",
  "chol_hdl",
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_ldl_p <- paste0(names(c_ldl_df[c_ldl_vars]), collapse = " + ")
c_ldl_lp <- paste0("chol_ldl ~ ", c_ldl_p, " + (1 | nct_id)")
c_ldl_fml <- as.formula(c_ldl_lp)
c_ldl_mdl <- lmer(formula = c_ldl_fml, data = c_ldl_df)
summary(c_ldl_mdl)

# Predict chol_ldl values
c_ldl_newdat <- c_hdl_inc %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("chol_ldl", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains(c("chol")), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-chol_ldl)

c_ldl_pred <- predict(
  c_ldl_mdl, 
  newdata = c_ldl_newdat, 
  allow.new.levels = TRUE
)

# Replace NA ldl with predictions
c_ldl_add <- c_ldl_newdat %>% 
  cbind(tibble(chol_ldl = c_ldl_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi, durn, smoke, dbp, sbp,
    chol_ldl,
    everything()
  ) %>%
  as_tibble()

# Incorporate imputed BMI
c_ldl_inc <- c_hdl_inc %>% 
  left_join(c_ldl_add %>% select(nct_id, id, chol_ldl_new = chol_ldl)) %>% 
  mutate(chol_ldl = coalesce(chol_ldl, chol_ldl_new)) %>% 
  select(-chol_ldl_new)

c_ldl_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after
c_hdl_inc %>% 
  select(chol_ldl) %>% 
  mutate(type = "before") %>% 
  full_join(c_ldl_inc %>% select(chol_ldl) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = chol_ldl, colour = type)) +
  geom_density() +
  theme_classic()

# Total cholesterol ------------------------------------------------------------

# 74/76 trials
# Add 0 for any other variables that are completely missing
c_total_df <- c_ldl_inc %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("chol_total", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(chol_total, dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  )

# Fit mixed-effects model
c_total_vars <- c(
  "age10", "sex", "race", "bmi", "durn", "smoke", "dbp", "sbp",
  "chol_hdl", "chol_ldl",
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_total_p <- paste0(names(c_total_df[c_total_vars]), collapse = " + ")
c_total_lp <- paste0("chol_total ~ ", c_total_p, " + (1 | nct_id)")
c_total_fml <- as.formula(c_total_lp)
c_total_mdl <- lmer(formula = c_total_fml, data = c_total_df)
summary(c_total_mdl)

# Predict chol_total values
c_total_newdat <- c_ldl_inc %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("chol_total", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(contains(c("chol")), dbp, sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-chol_total)

c_total_pred <- predict(
  c_total_mdl, 
  newdata = c_total_newdat, 
  allow.new.levels = TRUE
)

# Replace NA total with predictions
c_total_add <- c_total_newdat %>% 
  cbind(tibble(chol_total = c_total_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi, durn, smoke, dbp, sbp,
    chol_total,
    everything()
  ) %>%
  as_tibble()

# Incorporate imputed BMI
c_total_inc <- c_ldl_inc %>% 
  left_join(c_total_add %>% select(nct_id, id, chol_total_new = chol_total)) %>% 
  mutate(chol_total = coalesce(chol_total, chol_total_new)) %>% 
  select(-chol_total_new)

c_total_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after
c_ldl_inc %>% 
  select(chol_total) %>% 
  mutate(type = "before") %>% 
  full_join(c_total_inc %>% select(chol_total) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = chol_total, colour = type)) +
  geom_density() +
  theme_classic()

# Diastolic blood pressure -----------------------------------------------------

# 75/76 trials
# Add 0 for any other variables that are completely missing
c_dbp_df <- c_total_inc %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("dbp", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-sbp)

# Fit mixed-effects model
c_dbp_vars <- c(
  "age10", "sex", "race", "bmi", "durn", "smoke",
  "chol_hdl", "chol_ldl", "chol_total",
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_dbp_p <- paste0(names(c_dbp_df[c_dbp_vars]), collapse = " + ")
c_dbp_lp <- paste0("dbp ~ ", c_dbp_p, " + (1 | nct_id)")
c_dbp_fml <- as.formula(c_dbp_lp)
c_dbp_mdl <- lmer(formula = c_dbp_fml, data = c_dbp_df)
summary(c_dbp_mdl)

# Predict chol_dbp values with residual noise (sigma from model)
# Fails using predict() because of the loss of race reference
# Same approach is made simpler by specifying the model matrix
c_dbp_newdat <- c_total_inc %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("dbp", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  mutate(
    across(
      c(sbp),
      ~ if_else(is.na(.), 0, .)
    )
  ) %>% 
  select(-dbp)

c_dbp_pred <- predict(
  c_dbp_mdl, 
  newdata = c_dbp_newdat, 
  allow.new.levels = TRUE
)

# Replace NA dbp with predictions
c_dbp_add <- c_dbp_newdat %>% 
  cbind(tibble(dbp = c_dbp_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi, durn, smoke, dbp,
    everything()
  ) %>%
  as_tibble()

# Incorporate imputed BMI
c_dbp_inc <- c_total_inc %>% 
  left_join(c_dbp_add %>% select(nct_id, id, dbp_new = dbp)) %>% 
  mutate(dbp = coalesce(dbp, dbp_new)) %>% 
  select(-dbp_new)

c_dbp_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after
c_total_inc %>% 
  select(dbp) %>% 
  mutate(type = "before") %>% 
  full_join(c_dbp_inc %>% select(dbp) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = dbp, colour = type)) +
  geom_density() +
  theme_classic()

# Systolic blood pressure ------------------------------------------------------

# 75/76 trials
# Add 0 for any other variables that are completely missing
c_sbp_df <- c_dbp_inc %>% 
  inner_join(
    b_na %>% 
      filter(!all_na, grepl("sbp", name)) %>% 
      distinct(sponsor, nct_id)
  )

# Fit mixed-effects model
c_sbp_vars <- c(
  "age10", "sex", "race", "bmi", "durn", "smoke", "dbp",
  "chol_hdl", "chol_ldl", "chol_total",
  "hba1c", "egfr", "comorbs", "conmeds"
)

c_sbp_p <- paste0(names(c_sbp_df[c_sbp_vars]), collapse = " + ")
c_sbp_lp <- paste0("sbp ~ ", c_sbp_p, " + (1 | nct_id)")
c_sbp_fml <- as.formula(c_sbp_lp)
c_sbp_mdl <- lmer(formula = c_sbp_fml, data = c_sbp_df)
summary(c_sbp_mdl)

# Predict chol_dbp values with residual noise (sigma from model)
# Fails using predict() because of the loss of race reference
# Same approach is made simpler by specifying the model matrix
c_sbp_newdat <- c_dbp_inc %>% 
  inner_join(
    b_na %>% 
      filter(all_na, grepl("sbp", name)) %>% 
      distinct(sponsor, nct_id)
  ) %>% 
  select(-sbp)

c_sbp_pred <- predict(
  c_sbp_mdl, 
  newdata = c_sbp_newdat, 
  allow.new.levels = TRUE
)

# Replace NA dbp with predictions
c_sbp_add <- c_sbp_newdat %>% 
  cbind(tibble(sbp = c_sbp_pred)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft,
    age10, sex, bmi, durn, smoke, dbp, sbp,
    everything()
  ) %>%
  as_tibble()

# Incorporate imputed BMI
# Back-transform cholesterol
c_sbp_inc <- c_dbp_inc %>% 
  left_join(c_sbp_add %>% select(nct_id, id, sbp_new = sbp)) %>% 
  mutate(sbp = coalesce(sbp, sbp_new)) %>% 
  select(-sbp_new) %>% 
  mutate(across(contains("chol"), ~ exp(.)))

c_sbp_inc %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Compare before and after
c_dbp_inc %>% 
  select(sbp) %>% 
  mutate(type = "before") %>% 
  full_join(c_sbp_inc %>% select(sbp) %>% mutate(type = "after")) %>% 
  ggplot(aes(x = sbp, colour = type)) +
  geom_density() +
  theme_classic()

# Overall checks ---------------------------------------------------------------

a_imp_df %>% reframe(across(age10:race,  ~ sum(is.na(.))))
c_sbp_inc %>% reframe(across(age10:race,  ~ sum(is.na(.))))

a_imp_df %>% 
  reframe(
    across(
      c(age10, bmi:durn, dbp:conmeds),
      ~ unlist(
        list(
          mean = mean(., na.rm = TRUE), 
          min = min(., na.rm = TRUE), 
          max = max(., na.rm = TRUE)
        )
      )
    )
  ) %>% 
  mutate(sum = c("median", "min", "max")) %>% 
  select(sum, everything())

c_sbp_inc %>% 
  reframe(
    across(
      c(age10, bmi:durn, dbp:conmeds),
      ~ unlist(
        list(
          mean = mean(., na.rm = TRUE), 
          min = min(., na.rm = TRUE), 
          max = max(., na.rm = TRUE)
        )
      )
    )
  ) %>% 
  mutate(sum = c("median", "min", "max")) %>% 
  select(sum, everything())

# Compare distributions --------------------------------------------------------

# Overlay before and after
d_cmp <- a_imp_df %>% 
  select(
    age10, bmi, durn, dbp, sbp, 
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr, comorbs, conmeds
  ) %>% 
  pivot_longer(everything()) %>% 
  mutate(type = "before") %>% 
  full_join(
    c_sbp_inc %>% 
      select(
        age10, bmi, durn, dbp, sbp, 
        chol_hdl, chol_ldl, chol_total,
        hba1c, egfr, comorbs, conmeds
      ) %>% 
      pivot_longer(everything()) %>% 
      mutate(type = "after")
  ) %>% 
  arrange(name)

d_plot_cmp <- d_cmp %>% 
  ggplot(aes(x = value, colour = type)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  scale_fill_manual() +
  theme_classic()

ggsave(
  "Scratch_data/imputation/04_distributions_before_after_hier_imp.png",
  d_plot_cmp,
  width = 8,
  height = 6,
  units = "in"
)

# Final tidying ----------------------------------------------------------------

# Import numeric variable mu and d
e_imp_scale <- read_csv("Scratch_data/covariate_scales.csv")
  # filter(
  #   !name %in% c("age10", "hba1c", "chol_hdl", "chol_ldl", "chol_total")
  # )

# Function to centre variables on mean and scale by standard deviation
scale_centre_cov <- function(dat, e_imp_scale) {
  
  scale_lkp <- e_imp_scale %>% filter(name %in% names(dat))
  mu <- setNames(scale_lkp$mu, scale_lkp$name)
  d <- setNames(scale_lkp$d, scale_lkp$name)
  
  dat %>% 
    mutate(
      across(
        all_of(scale_lkp$name),
        ~ (.x - mu[cur_column()]) / d[cur_column()]
      )
    )
}

# Nest data into two versions
# Original
# Variables centred on mean and scaled by standard deviation
e_tidy <- c_sbp_inc %>% 
  group_by(sponsor, nct_id) %>% 
  nest() %>% 
  mutate(
    centred_scaled_data = map(data, ~ scale_centre_cov(.x, e_imp_scale))
  )

# Plot scaled variables
e_plot <- e_tidy %>% 
  select(sponsor, nct_id, centred_scaled_data) %>% 
  unnest() %>% 
  ungroup() %>% 
  pivot_longer(c(age10, bmi:durn, dbp:conmeds)) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free") + 
  theme_classic()

ggsave(
  "Scratch_data/imputation/05_distributions_scaled_vars.png",
  e_plot,
  width = 8,
  height = 6,
  units = "in"
)

# Save
saveRDS(e_tidy, "Processed_data/05_hier_imputed.rds")
