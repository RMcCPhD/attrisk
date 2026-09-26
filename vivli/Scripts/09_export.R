
# Prepare outputs for export

source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)

# Model outputs ----------------------------------------------------------------

# Import model outputs
a_imp_out <- readRDS("Processed_data/06_cox_out.rds")

# NULL estimates are always race/ethnicity (predominantly other)
a_imp_out %>% 
  select(sponsor, nct_id, res) %>% 
  unnest(res) %>%
  ungroup() %>% 
  filter(is.na(estimate)) %>% 
  count(term)

# Prepare estimates
# Round to 3dp, remove NULL estimates
b_res <- a_imp_out %>% 
  select(sponsor, nct_id, res) %>% 
  unnest(res) %>%
  ungroup() %>% 
  mutate(across(-c(sponsor, nct_id, term), ~ round(., 3))) %>% 
  filter(!is.na(estimate))

# Prepare diagnostics
b_diag <- a_imp_out %>% 
  select(sponsor, nct_id, diag) %>% 
  unnest(diag) %>% 
  ungroup() %>% 
  mutate(across(-c(sponsor, nct_id), ~ round(., 3)))

# Convert VCOV to correlation matrices
# Keep upper triangle
b_cov2cor <- a_imp_out %>% 
  select(sponsor, nct_id, vcov) %>% 
  ungroup() %>% 
  mutate(
    cor = map(vcov, cov2cor),
    cor_upper = map(cor, ~ {
      .x[!upper.tri(.x)] <- NA_real_
      .x
    })
  )

# Get long-format correlation tibble (upper-triangle only)
b_cor <- b_cov2cor %>% 
  mutate(
    cor_upper_long = map(cor_upper, ~ {
      tbl <- as_tibble(as.data.frame(.x), rownames = "row")
      tbl_long <- tbl %>% gather("col", "value", -row)
      tbl_filter <- tbl_long %>% filter(!is.na(value))
      tbl_filter
    })
  ) %>% 
  select(sponsor, nct_id, cor_upper_long) %>% 
  unnest(cor_upper_long) %>% 
  ungroup() %>% 
  mutate(value = round(value, 3))

# Save exportable model outputs
write_csv(b_res, "Outputs/export_res.csv")
write_csv(b_diag, "Outputs/export_diag.csv")
write_csv(b_cor, "Outputs/export_covcor.csv")

# Summary statistics -----------------------------------------------------------

# Import analysis data
b_imp_df <- readRDS("Processed_data/05_hier_imputed.rds") %>% 
  select(sponsor, nct_id, centred_scaled_data) %>% 
  unnest(centred_scaled_data) %>% 
  ungroup()

# Sponsor count (for readme in export)
b_imp_df %>% 
  distinct(sponsor, nct_id) %>% 
  group_by(sponsor) %>% 
  reframe(
    n = n()
  ) %>% 
  reframe(
    sponsor = sponsor,
    freq = paste0(n , " (", round(n / sum(n) * 100, 2), "%)")
  )

# Import covariate scaling (reverse for baseline summaries)
b_imp_scale <- read_csv("Scratch_data/covariate_scales.csv")

# Function to revert scaling
unscale_cov <- function(dat, e_imp_scale) {
  
  scale_lkp <- e_imp_scale %>% filter(name %in% names(dat))
  mu <- setNames(scale_lkp$mu, scale_lkp$name)
  d <- setNames(scale_lkp$d, scale_lkp$name)
  
  dat %>% 
    mutate(
      across(
        all_of(scale_lkp$name),
        ~ (.x + mu[cur_column()]) * d[cur_column()]
      )
    )
}

# Unscale variables
b_unscale <- b_imp_df %>% 
  group_by(sponsor, nct_id) %>% 
  nest() %>% 
  mutate(
    unscale_data = map(data, ~ unscale_cov(.x, b_imp_scale))
  ) %>% 
  select(-data) %>% 
  unnest() %>% 
  ungroup()

# Trial-level counts
# Attrition (arm, attr (n (%)))
# Person time (years), median [min, max]
b_attr_trial <- b_imp_df %>% 
  group_by(nct_id, trtcls5, trttype) %>% 
  reframe(
    n = n_distinct(id),
    attr = sum(event),
    attr_p = round(attr / n * 100, 2)
  ) %>% 
  left_join(
    b_imp_df %>% 
      group_by(nct_id, trtcls5, trttype) %>% 
      reframe(
        pt_yr = days / 365.25
      ) %>% 
      group_by(nct_id) %>% 
      reframe(
        med_pt = round(median(pt_yr), 3),
        min_pt = round(min(pt_yr), 3),
        max_pt = round(max(pt_yr), 3)
      ) %>% 
      group_by(nct_id) %>% 
      reframe(
        pt = paste0(med_pt, " [", min_pt, ", ", max_pt, "]")
      )
  )

# Trial-level summaries
# Number of trials and participants
# Attrition risk (median [min, max])
# Person time (median [min, max])
# Follow-up duration (median [min, max])
b_trial_sum <- b_imp_df %>% 
  reframe(n_trials = n_distinct(nct_id)) %>% 
  cbind(
    b_imp_df %>% 
      group_by(nct_id) %>% 
      reframe(n = n_distinct(id)) %>% 
      reframe(n_subj = sum(n))
  ) %>% 
  cbind(
    b_attr_trial %>% 
      reframe(
        med_risk = median(attr_p),
        min_risk = min(attr_p),
        max_risk = max(attr_p)
      ) %>% 
      reframe(
        risk = paste0(med_risk, " [", min_risk, ", ", max_risk, "]")
      )
  ) %>% 
  cbind(
    b_imp_df %>% 
      group_by(nct_id) %>% 
      reframe(
        pt_yr = days / 365.25
      ) %>% 
      reframe(
        med_pt = round(median(pt_yr), 3),
        min_pt = round(min(pt_yr), 3),
        max_pt = round(max(pt_yr), 3)
      ) %>% 
      reframe(
        pt = paste0(med_pt, " [", min_pt, ", ", max_pt, "]")
      )
  ) %>% 
  cbind(
    b_imp_df %>% 
      group_by(nct_id) %>% 
      reframe(fu = max(weeks)) %>% 
      reframe(
        med_fu = median(fu),
        min_fu = min(fu),
        max_fu = max(fu)
      ) %>% 
      reframe(
        followup = paste0(med_fu, " [", min_fu, ", ", max_fu, "]")
      )
  ) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(everything())

# Baseline characteristics
# Age (years) - mean (SD)
# N (%) female
# Race/ethnicity - n (%)
# Smoking status - n (%)
# Physical measures (bmi, sbp) - mean (SD)
# Laboratory measures (hba1c, cholesterol) - mean (SD)
# Medical history (t2dm duration, comorbidities, con meds) - mean (SD), median
b_char_sum <- b_unscale %>% 
  reframe(
    age_mu = mean(age10),
    age_d = sd(age10),
    age = paste0(round(age_mu, 2), " [", round(age_d, 2), "]")
  ) %>% 
  select(age) %>% 
  cbind(
    b_unscale %>% 
      reframe(
        female_n = sum(sex == 0L),
        female_p = round(female_n / n() * 100, 2),
        female = paste0(female_n, " (", female_p, "%)")
      ) %>% 
      select(female)
  ) %>% 
  cbind(
    b_unscale %>% 
      mutate(
        race = if_else(race == "other" | race == "na_unk", "other", race),
        race = paste0("race: ", race)
      ) %>% 
      group_by(race) %>% 
      reframe(race_n = n()) %>% 
      reframe(
        race = race,
        race_n = race_n,
        race_p = round(race_n / sum(race_n) * 100, 2)
      ) %>% 
      reframe(
        race = race,
        freq = paste0(race_n, " (", race_p, "%)")
      ) %>% 
      pivot_wider(names_from = "race", values_from = "freq")
  ) %>% 
  cbind(
    b_unscale %>% 
      mutate(
        smoke = case_match(
          smoke,
          0L ~ "smoke: former/never", 
          1L ~ "smoke: current"
        ),
      ) %>% 
      group_by(smoke) %>% 
      reframe(smoke_n = n()) %>% 
      reframe(
        smoke = smoke,
        smoke_n = smoke_n,
        smoke_p = round(smoke_n / sum(smoke_n) * 100, 2)
      ) %>% 
      reframe(
        smoke = smoke,
        freq = paste0(smoke_n, " (", smoke_p, "%)")
      ) %>% 
      pivot_wider(names_from = "smoke", values_from = "freq")
  ) %>% 
  cbind(
    b_unscale %>% 
      reframe(
        bmi_mu = round(mean(bmi), 2),
        bmi_d = round(sd(bmi), 2),
        bmi = paste0(bmi_mu , " [", bmi_d, "]"),
        sbp_mu = round(mean(sbp), 2),
        sbp_d = round(sd(sbp), 2),
        sbp = paste0(sbp_mu , " [", sbp_d, "]"),
        hba1c_mu = round(mean(hba1c), 2),
        hba1c_d = round(sd(hba1c), 2),
        hba1c = paste0(hba1c_mu , " [", hba1c_d, "]"),
        chol_ldl_mu = round(mean(chol_ldl), 2),
        chol_ldl_d = round(sd(chol_ldl), 2),
        chol_ldl = paste0(chol_ldl_mu , " [", chol_ldl_d, "]"),
        chol_total_mu = round(mean(chol_total), 2),
        chol_total_d = round(sd(chol_total), 2),
        chol_total = paste0(chol_total_mu , " [", chol_total_d, "]"),
        egfr_mu = round(mean(egfr), 2),
        egfr_d = round(sd(egfr), 2),
        egfr = paste0(egfr_mu , " [", egfr_d, "]")
      ) %>% 
      select(bmi, sbp, hba1c, chol_ldl, chol_total, egfr)
  ) %>% 
  cbind(
    b_unscale %>% 
      reframe(
        across(
          c(durn, comorbs, conmeds),
          list(
            med = ~ round(median(., na.rm = TRUE), 2),
            min = ~ round(min(., na.rm = TRUE), 2),
            max = ~ round(max(., na.rm = TRUE), 2)
          )
        )
      ) %>% 
      reframe(
        durn = paste0(durn_med, " [", durn_min, ", ", durn_max, "]"),
        comorbs = paste0(comorbs_med, " [", comorbs_min, ", ", comorbs_max, "]"),
        conmeds = paste0(conmeds_med, " [", conmeds_min, ", ", conmeds_max, "]")
      )
  ) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(everything())
  
# Combine trial-level summaries
b_trial_all <- b_trial_sum %>% full_join(b_char_sum)
  
# Save summaries
write_csv(b_attr_trial, "Outputs/trial_level_summary.csv")
write_csv(b_trial_all, "Outputs/overall_summary.csv")

# Missingness/Imputation exercises ---------------------------------------------

# Import datasets from each stage
# No imputation - within trial - across trials
a_imp_hier <- readRDS("Processed_data/03_mt10_events_n76.rds")
a_imp_hier1 <- readRDS("Processed_data/04_within_trial_imputed.rds")
a_imp_hier2 <- readRDS("Processed_data/05_hier_imputed.rds")

# Summarise missingness at each stage of imputation
# na_start = Beginning (no imputation)
# na_stg1 = Within-trial imputation for partially missing covariates
# na_stg2 = Across-trial imputation for completely missing covariates
c_na_sum <- a_imp_hier %>% 
  select(any_of(unique(b_res$term)), race) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  reframe(
    n = n(),
    na = sum(is.na(value))
  ) %>% 
  reframe(
    term = name,
    na_start = paste0(na, " (", round(na / n * 100, 2), "%)")
  ) %>% 
  left_join(
    a_imp_hier1 %>% 
      select(any_of(unique(b_res$term)), race) %>% 
      mutate(across(everything(), as.character)) %>% 
      pivot_longer(everything()) %>% 
      group_by(name) %>% 
      reframe(
        n = n(),
        na = sum(is.na(value))
      ) %>% 
      reframe(
        term = name,
        na_stg1 = paste0(na, " (", round(na / n * 100, 2), "%)")
      )
  ) %>% 
  left_join(
    a_imp_hier2 %>% 
      select(sponsor, nct_id, data) %>% 
      unnest() %>% 
      ungroup() %>% 
      select(any_of(unique(b_res$term)), race) %>% 
      mutate(across(everything(), as.character)) %>% 
      pivot_longer(everything()) %>% 
      group_by(name) %>% 
      reframe(
        n = n(),
        na = sum(is.na(value))
      ) %>% 
      reframe(
        term = name,
        na_stg2 = paste0(na, " (", round(na / n * 100, 2), "%)")
      )
  )

# Save
write_csv(c_na_sum, "Outputs/missingness_imputation_stages.csv")
