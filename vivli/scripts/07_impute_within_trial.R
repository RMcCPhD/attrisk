
# Missing data imputed using mice within-trial
# Imputation for variables that are not completely missing
# Hierarchical models will be fit for variables that are completely missing
# Refit models from mice (outputs not saved by function)
# Save coefficients from refit models used for imputation
source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(broom)
library(mice)

# readinscales
# myscale <- read_csv( "../model_covariates/Outputs/scaling_values.csv")

# Import data
a_imp_df <- readRDS("Processed_data/n77_cox.rds")
glimpse(a_imp_df)

# Import original covariates
t2d_path <- "V:/From 8697/work_folder/Extract_data/"
t2d_path2 <- "Missing_covars/t2d_covariates_n120.rds"
t2d_path_pasted <- paste0(t2d_path, t2d_path2)
a_imp_covs <- readRDS(t2d_path_pasted)

# Get relevant covariates
b_covs <- a_imp_covs %>% 
  select(
    nct_id, id = usubjid,
    age, sex, smoke, race, bmi, durdiabyr,
    bdiase:hba1c_pcnt, egfr,
    chol_hdl:chol_total
  )

# Join covariates with 77 trial groupings
b_join <- a_imp_df %>% 
  select(nct_id, id, trttype, weeks, event, comorbs, n_atc) %>% 
  left_join(b_covs) %>% 
  select(nct_id:event, age:chol_total, everything())

b_join %>% 
  pivot_longer(c(age, bmi:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "skyblue") +
  facet_wrap(~ var, scales = "free")

plot_categs <- b_join %>% 
  pivot_longer(sex:race, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_bar(fill = "skyblue", colour = "black") +
  facet_wrap(~ var, scales = "free") +
  coord_flip()

# Prepare for imputation
# Sex as integer (1 - male, 0 - female)
# Smoking status as integer (1 - current, 0 - previous/never/unknown)
# Factorise race before cox model, not here
c_prep <- b_join %>% 
  mutate(
    sex = case_when(sex == "Male" ~ 1L, TRUE ~ 0L),
    smoke = case_when(smoke == "Current" ~ 1L, TRUE ~ 0L)
  )

c_prep %>% reframe(across(age:n_atc, ~ sum(is.na(.))))

# Add hba1c
hba1c_path <- "V:/From 8697/work_folder/model_covariates/"
hba1c_path2 <- "Scratch_data/data_for_hba1c_main_models.Rds"
hba1c_path_pasted <- paste0(hba1c_path, hba1c_path2)

# Originally took durn from here but already scaled
# Taking original extraction version
c_hba1c <- readRDS(hba1c_path_pasted) %>% 
  select(sponsor:nct_id, data) %>% 
  unnest(data) %>% 
  select(nct_id, id, age10, hba1c = value_1)

# Add durn (unscaled from covariates_with_imputation)
durn_path <- "V:/From 8697/work_folder/model_covariates/"
durn_path2 <- "Scratch_data/covariates_with_imputation.rds"
durn_path_pasted <- paste0(durn_path, durn_path2)

# Originally took durn from here but already scaled
# Taking original extraction version
c_durn <- readRDS(durn_path_pasted) %>% select(nct_id, id, durn)

d_df <- c_prep %>% 
  select(-hba1c_pcnt, -age, -durdiabyr) %>% 
  left_join(c_hba1c) %>% 
  left_join(c_durn) %>% 
  select(nct_id:event, age10, sex:bmi, hba1c, durn, everything())

d_df %>% reframe(across(age10:n_atc, ~ sum(is.na(.))))

# Remove values where there's only one value and rest are NA
# Unreliable imputation, imputed by hierarchical models
d_rm_singles <- d_df %>% 
  mutate(
    egfr = case_when(nct_id == "NCT00975286" ~ NA, TRUE ~ egfr),
    chol_hdl = case_when(nct_id == "NCT00975286" ~ NA, TRUE ~ chol_hdl),
    chol_ldl = case_when(nct_id %in% c("NCT00975286", "NCT01558271") ~ NA, TRUE ~ chol_ldl),
    chol_total = case_when(nct_id == "NCT00975286" ~ NA, TRUE ~ chol_total),
    across(comorbs:n_atc, as.integer)
  )

d_rm_singles %>% reframe(across(age10:n_atc, ~ sum(is.na(.))))

# Add bmi
bmi_path <- "V:/From 8697/work_folder/model_covariates/"
bmi_path2 <- "Scratch_data/covariates_with_imputation.Rds"
bmi_path_pasted <- paste0(bmi_path, bmi_path2)

e_bmi <- readRDS(bmi_path_pasted) %>% select(nct_id, id, bmi) 

e_df <- d_rm_singles %>% 
  select(-bmi) %>% 
  left_join(e_bmi) %>% 
  rename(
    dbp = bdiase,
    sbp = bsysse
  ) %>% 
  select(nct_id:event, age10:durn, bmi, everything())

e_df %>% reframe(across(age10:n_atc, ~ sum(is.na(.))))

# Check scaling of imputation-ready data
e_df %>% 
  pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "skyblue") +
  facet_wrap(~ var, scales = "free")

# Save my scalings
e_df %>% 
  select(bmi, durn, egfr, sbp) %>% 
  reframe(across(everything(), ~ sd(., na.rm = TRUE))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "scaling_sd") %>% 
  write_csv("Scratch_data/attrisk_var_scaling.csv")

# Scale variables scaled in other analyses (based on current SDs)
e_df_scaled <- e_df %>% 
  mutate(across(c(bmi, durn, egfr, sbp), ~ . / sd(., na.rm = TRUE)))

# Check new scaling
e_df_scaled %>% 
  pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_histogram(fill = "skyblue") +
  facet_wrap(~ var, scales = "free")
  
# Imputation for partially missing variables -----------------------------------

# Split into lists by nct_id
imp_lst <- split(e_df_scaled, e_df_scaled$nct_id)

# Remove completely NA variables (hierarchical imputation afterward)
imp_rm <- lapply(imp_lst, function(df) select(df, where(~ !(all(is.na(.x))))))

# Impute
imp_run <- lapply(imp_rm, function(df) {
  df_vars <- df %>% select(-c(nct_id:event))
  run_mic <- mice(df_vars, m = 20, method = "pmm", seed = 123)
  return(run_mic)
})

# Logged events
imp_log <- lapply(imp_run, function(df) df$loggedEvents)
  
# Extract imputed values
imp_ext <- lapply(imp_run, function(df) complete(df, 1))

# Bind rows with completely NA variables untouched
imp_df <- bind_rows(imp_ext) %>% 
  cbind(bind_rows(imp_rm) %>% select(nct_id:trttype)) %>% 
  select(nct_id:trttype, everything()) %>% 
  as_tibble()

imp_df %>% reframe(across(everything(), ~ sum(is.na(.))))

# Compare distributions between imputed and pre-imputation
plot_dists <- imp_df %>% 
  pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density() +
  geom_density(
    data = e_df_scaled %>% 
      pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val"),
    aes(x = val),
    colour = "forestgreen", alpha = 0.2
  ) +
  facet_wrap(~var, scales = "free")

# Save outputs -----------------------------------------------------------------

# Plot of categorical variables
ggsave(
  "Outputs/imputation/plot_categs.png",
  plot_categs,
  width = 8,
  height = 6,
  units = "in"
)
# Plots of distribution
ggsave(
  "Outputs/imputation/plot_before_after_mice.png",
  plot_dists,
  width = 8,
  height = 6,
  units = "in"
)

# MICE-imputed dataset
saveRDS(imp_df, "Processed_data/imputation_part1.rds")
