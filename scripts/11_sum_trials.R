
# Trial-level summary statistics
source("scripts/00_config.R")
library(tidyverse)

# Import AgD
a_imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]] %>% 
  ungroup()

# Import IPD arm-level counts
a_imp_ipd <- read_csv("vivli/Summaries/trial_level_summary.csv")

# Identify source from original data extract
a_imp_orig <- readRDS("data/tidied_agg.rds")

# Total participants
sum(a_imp_agd$arm_n)

# Follow-up duration
sum_fu <- a_imp_agd %>% 
  mutate(duration = duration * 52) %>% 
  reframe(
    med_fu = median(duration),
    min_fu = min(duration),
    max_fu = max(duration),
    across(med_fu:max_fu, ~ round(.))
  ) %>% 
  mutate(
    value = paste0(med_fu, " [", min_fu, ", ", max_fu, "]")
  )

sum_fu_ipd <- a_imp_agd %>% 
  filter(trial_id %in% a_imp_ipd$nct_id) %>% 
  mutate(duration = duration * 52) %>% 
  reframe(
    med_fu = median(duration),
    min_fu = min(duration),
    max_fu = max(duration),
    across(med_fu:max_fu, ~ round(.))
  ) %>% 
  mutate(
    value = paste0(med_fu, " [", min_fu, ", ", max_fu, "]")
  )

# Attrition
sum_attr <- a_imp_agd %>% 
  reframe(
    total = sum(arm_n),
    attr = sum(arm_attr),
    attr_p = (attr / total) * 100,
    across(total:attr_p, ~ round(., 2))
  ) %>% 
  mutate(
    value = paste0(attr, " [", attr_p, "%]")
  )

sum_attr_risk <- a_imp_agd %>% 
  mutate(attr_p = (arm_attr / arm_n) * 100) %>% 
  reframe(
    med_attr = median(attr_p),
    min_attr = min(attr_p),
    max_attr = max(attr_p),
    across(med_attr:max_attr, ~ round(., 2))
  ) %>% 
  mutate(
    value = paste0(med_attr, " [", min_attr, ", ", max_attr, "]")
  )

sum_attr_ipd <- a_imp_ipd %>% 
  reframe(
    total = sum(n),
    attr_s = sum(attr),
    attr_ps = (attr_s / total) * 100,
    across(total:attr_ps, ~ round(., 1))
  ) %>% 
  mutate(
    value = paste0(attr_s, " [", attr_ps, "%]")
  )

sum_attr_risk_ipd <- a_imp_ipd %>% 
  mutate(attr_p = (attr / n) * 100) %>% 
  reframe(
    med_attr = median(attr_p),
    min_attr = min(attr_p),
    max_attr = max(attr_p),
    across(med_attr:max_attr, ~ round(., 2))
  ) %>% 
  mutate(
    value = paste0(med_attr, " [", min_attr, ", ", max_attr, "]")
  )


# Trial-level attrition counts
sum_attr_trial <- a_imp_agd %>% 
  select(trial_id, trttype, arm_n, arm_attr) %>% 
  mutate(attr_p = round(arm_attr / arm_n * 100, 2)) %>% 
  mutate(attr = paste0(arm_attr, " (", attr_p, "%)")) %>% 
  select(-c(attr_p, arm_attr))

sum_attr_ipd <- a_imp_ipd %>% 
  mutate(attr = paste0(attr, " (", attr_p, "%)")) %>% 
  # select(-c(trtcls5, attr_p)) %>% 
  rename(arm_n = n, total_ptime_yr = pt_yr_total) %>% 
  mutate(
    trttype = recode(
      trttype,
      "plc" = "placebo",
      "ac" = "active"
    )
  )

a_imp_ipd %>% 
  group_by(nct_id) %>% 
  reframe(total = sum(n)) %>% 
  reframe(
    med = median(total),
    min = min(total),
    max = max(total)
  )

sum(a_imp_ipd$n)

# Trial characteristics
sum_trial_chars <- a_imp_agd %>% 
  select(-c(trial_id, source, arm_n:duration)) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(everything()) %>% 
  group_by(name, value) %>% 
  reframe(n = n()) %>% 
  group_by(name) %>% 
  mutate(p = round(n / sum(n) * 100, 2))

# Save
sum_save <- mget(ls(pattern = "sum_"))
lapply(names(sum_save), function(nm) {
  save_dir <- paste0("outputs/trial_summaries/", nm, ".csv")
  write_csv(sum_save[[nm]], save_dir)
})

# IPD characteristics
imp_sum_bl <- read_csv("outputs/conduct_stage/ipd_bl_sum.csv")

# Helper to summarise as pooled distribution summaries across arms
pool_continuous <- function(data, group = NULL) {
  
  if (is.null(group)) group <- character(0)
  
  cont_vars <- names(data) %>%
    str_subset("_mu$") %>%
    str_remove("_mu$") %>%
    keep(~ paste0(.x, "_d") %in% names(data))
  
  data %>%
    group_by(across(all_of(group))) %>%
    group_modify(~ {
      
      map_dfr(cont_vars, function(v) {
        
        mu_col <- paste0(v, "_mu")
        sd_col <- paste0(v, "_d")
        
        dat_v <- .x %>%
          filter(
            !is.na(n),
            !is.na(.data[[mu_col]]),
            !is.na(.data[[sd_col]]),
            n > 0
          )
        
        N <- sum(dat_v$n)
        
        if (N <= 1) {
          return(tibble(
            variable = v,
            n = N,
            mean = NA_real_,
            sd = NA_real_
          ))
        }
        
        pooled_mean <- weighted.mean(dat_v[[mu_col]], w = dat_v$n)
        
        pooled_var <- sum(
          (dat_v$n - 1) * dat_v[[sd_col]]^2 +
            dat_v$n * (dat_v[[mu_col]] - pooled_mean)^2
        ) / (N - 1)
        
        tibble(
          variable = v,
          n = N,
          mean = pooled_mean,
          sd = sqrt(pooled_var)
        )
      })
    }) %>%
    ungroup()
}

pool_proportions <- function(data) {
  
  data %>%
    select(nct_id, n, cursmoke, white, black_aa, asian, other, female_p) %>%
    tidyr::pivot_longer(
      cols = c(cursmoke, white, black_aa, asian, other, female_p),
      names_to = "variable",
      values_to = "proportion_arm"
    ) %>%
    filter(
      !is.na(n),
      !is.na(proportion_arm),
      n > 0
    ) %>%
    group_by(nct_id, variable) %>%
    summarise(
      n = sum(n),
      count = round(sum(n * proportion_arm)),
      proportion = count / n,
      .groups = "drop"
    ) %>% 
    group_by(variable) %>% 
    reframe(
      total = sum(count),
      prop = round(total / 19023 * 100, 2)
    )
}

pool_continuous(imp_sum_bl)
pool_proportions(imp_sum_bl)

# Summarise IPD trial-specific Cox model outputs -------------------------------

# Import model results and arm summaries
imp_res <- read_csv("vivli/Models/export_res.csv")
imp_sum <- read_csv("vivli/Summaries/trial_level_summary.csv")

ipd_mdl_sum <- imp_sum %>% 
  distinct(nct_id, trttype) %>% 
  arrange(nct_id, trttype) %>%
  group_by(nct_id) %>% 
  slice_tail(n = 1) %>% 
  ungroup() %>% 
  mutate(
    trttype = recode(
      trttype,
      "plc" = "Placebo",
      "ac" = "Active"
    )
  ) %>% 
  rename(reference = trttype) %>%
  left_join(imp_res) %>% 
  filter(!term == "trttype") %>% 
  mutate(
    res = paste0(
      round(estimate, 3), " [",
      round(estimate - 1.96 * std.error, 3),
      ", ",
      round(estimate + 1.96 * std.error, 3),
      " ]"
    ),
    term = factor(
      term,
      levels = c(
        "age10", "sex", "raceasian", "raceblack_aa", "raceother_na_unk",
        "bmi", "durn", "smoke", "hba1c", "egfr", "sbp", "chol_hdl", "chol_total",
        "comorbs", "conmeds"
      )
    ),
    term = recode(
      term,
      "age10" = "Age (in decades)",
      "sex" = "Male",
      "bmi" = "BMI (kg/m2)",
      "durn" = "T2DM duration (years)",
      "smoke" = "Current smoker",
      "sbp" = "Systolic blood pressure (mmHg)",
      "chol_hdl" = "HDL cholesterol (mmol/L)",
      "chol_total" = "Total cholesterol (mmol/L)",
      "hba1c" = "HbA1c (%/DCCT units)",
      "egfr" = "eGFR (ml/min/1.73m2)",
      "comorbs" = "Comorbidities",
      "conmeds" = "Concomitant medications",
      "raceasian" = "Asian",
      "raceblack_aa" = "Black/African American",
      "raceother_na_unk" = "Other race/ethnicity"
    ),
    p.value = case_when(
      p.value < 0.01 ~ "<0.01",
      TRUE ~ as.character(round(p.value, 3))
    )
  ) %>% 
  select(
    nct_id, reference, term, res, p = p.value
  ) %>% 
  arrange(nct_id, term)

write_csv(ipd_mdl_sum, "outputs/trial_summaries/ipd_cox_mdls.csv")
