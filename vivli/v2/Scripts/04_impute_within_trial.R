
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

# Import data
a_imp_df <- readRDS("Processed_data/03_mt10_events_n76.rds")

# Check missingness (overall and trial-level)
a_imp_df %>% 
  reframe(across(everything(), ~sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Prepare data for within-trial imputation
# Sex already an integer (0 - female, reference)
# Convert smoking status to integer (0 - former/never, 1 - current)
# Simplify race to asian/black/white/other
# Comorb and conmed counts as integers
b_prep <- a_imp_df %>% 
  mutate(
    smoke = case_when(smoke == "Current" ~ 1L, TRUE ~ 0L),
    race = case_when(
      race %in% c("white", "black_aa", "asian") ~ race,
      race %in% c("ai_an", "other") ~ "other",
      TRUE ~ race
    ),
    # race = factor(race, levels = c("white", "asian", "black_aa", "other")),
    across(comorbs:conmeds, as.integer)
  )

# Compute variable scales (for use AFTER all imputation is complete, including
# hierarchical model-based imputation in next script)
b_scales <- b_prep %>% 
  select(
    age10, bmi, durn, dbp, sbp, 
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr
  ) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  reframe(
    mu = mean(value, na.rm = TRUE),
    d = sd(value, na.rm = TRUE)
  )

write_csv(b_scales, "Scratch_data/covariate_scales.csv")

# Perform imputation -----------------------------------------------------------

# Check variable distributions before imputation
c_plot_dist_b4 <- b_prep %>% 
  select(
    age10, bmi, durn, dbp, sbp, 
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr, comorbs, conmeds
  ) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) +
  geom_density(colour = "black", fill = "skyblue", alpha = 0.3) +
  facet_wrap(~name, scales = "free") +
  theme_classic()

ggsave(
  "Scratch_data/imputation/01_distributions_before_mice.png",
  c_plot_dist_b4,
  width = 8,
  height = 6,
  units = "in"
)

# Nest by groups
# Remove completely NA variables
c_nest <- b_prep %>% 
  group_by(sponsor, nct_id) %>% 
  nest() %>% 
  mutate(
    rm_na = map(
      data, 
      ~ .x %>% select(where(~!all(is.na(.))))
    )
  )

# Perform imputation (20 datasets, seed set, pmm method)
c_impute <- c_nest %>% 
  mutate(
    imp = map(rm_na, ~ {
      
      df <- .x %>% 
        select(-c(
          id, arm, trtcls5, trtcls, trttype, 
          days, weeks, event, category, pyr, mft
        ))
      
      mice(df, m = 20, method = "pmm", seed = 123)
      
    })
  )

# Extract logged events and first imputated dataset
c_log <- c_impute %>% 
  mutate(
    logged = map(imp, ~ .x$loggedEvents)
  ) %>% 
  select(sponsor, nct_id, logged) %>% 
  unnest(logged) %>% 
  ungroup()

# Events were:
# Collinearity amongst cholesterol variables (mostly LDL)
# Constant values for smoke, race (mostly), comorbs and conmeds
c_log %>% count(meth)
c_log %>% filter(meth == "pmm") %>% count(dep)
c_log %>% filter(meth == "constant") %>% count(out)

# Extract first imputation dataset
c_ext_imp <- c_impute %>% 
  mutate(
    ext = map(imp, ~ complete(.x, 1))
  ) %>% 
  select(sponsor, nct_id, ext) %>% 
  unnest(ext) %>% 
  ungroup()

# Join with variables excluded from imputation
c_join <- b_prep %>% 
  select(
    id, arm, trtcls5, trtcls, trttype, 
    days, weeks, event, category, pyr, mft
  ) %>% 
  cbind(c_ext_imp) %>% 
  as_tibble() %>% 
  select(sponsor, nct_id, everything())

# Plot post-imputation distributions
c_plot_dist_aft <- c_join %>% 
  select(
    age10, bmi, durn, dbp, sbp, 
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr, comorbs, conmeds
  ) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) +
  geom_density(colour = "black", fill = "skyblue", alpha = 0.3) +
  facet_wrap(~name, scales = "free") +
  theme_classic()

ggsave(
  "Scratch_data/imputation/02_distributions_after_mice.png",
  c_plot_dist_aft,
  width = 8,
  height = 6,
  units = "in"
)

# Overlay before and after
c_cmp <- b_prep %>% 
  select(
    age10, bmi, durn, dbp, sbp, 
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr, comorbs, conmeds
  ) %>% 
  pivot_longer(everything()) %>% 
  mutate(type = "before") %>% 
  full_join(
    c_join %>% 
      select(
        age10, bmi, durn, dbp, sbp, 
        chol_hdl, chol_ldl, chol_total,
        hba1c, egfr, comorbs, conmeds
      ) %>% 
      pivot_longer(everything()) %>% 
      mutate(type = "after")
  ) %>% 
  arrange(name)

c_plot_cmp <- c_cmp %>% 
  ggplot(aes(x = value, colour = type)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  scale_fill_manual() +
  theme_classic()

ggsave(
  "Scratch_data/imputation/03_distributions_before_after.png",
  c_plot_cmp,
  width = 8,
  height = 6,
  units = "in"
)

# Check missingness before and after as final check
a_imp_df %>% 
  reframe(across(everything(), ~sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

c_join %>% 
  reframe(across(everything(), ~sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)


# Check trial-level missingness
a_imp_df %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(n = length(unique(id))) %>% 
  left_join(
    a_imp_df %>% 
      group_by(sponsor, nct_id) %>% 
      reframe(
        across(
          c(
            age10, sex, bmi, durn, smoke, race, dbp, sbp,
            chol_hdl, chol_ldl, chol_total, hba1c, egfr, comorbs, conmeds
          ),
          ~ sum(is.na(.))
        )
      )
  ) %>% 
  print(n = Inf)

c_join %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(n = length(unique(id))) %>% 
  left_join(
    c_join %>% 
      group_by(sponsor, nct_id) %>% 
      reframe(
        across(
          c(
            age10, sex, bmi, durn, smoke, race, dbp, sbp,
            chol_hdl, chol_ldl, chol_total, hba1c, egfr, comorbs, conmeds
          ),
          ~ sum(is.na(.))
        )
      )
  ) %>% 
  print(n = Inf)

# Save dataset with within-trial imputed values
saveRDS(c_join, "Processed_data/04_within_trial_imputed.rds")
