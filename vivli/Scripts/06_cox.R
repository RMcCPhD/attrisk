
# Fit Cox proportional hazard models to each trial
# Estimates the effect of baseline participant characteristics on attrition risk

# Workflow
# Double-check dataset is ready for model (e.g. variable formats, exclusions)
# Fit Cox models
# Inspect and address convergence issues and other warnings
# Summarise outputs (estimates, diagnostics, prop. hazard checks)

source("Scripts/00_config.R")
source("Scripts/00_functions.R")

set.seed(123)
library(tidyverse)
library(broom)
library(survival)

# Import prepared data
a_imp_df <- readRDS("Processed_data/05_hier_imputed.rds") %>% 
  select(sponsor, nct_id, centred_scaled_data) %>% 
  unnest(centred_scaled_data) %>% 
  ungroup()

# Data prep --------------------------------------------------------------------

# Collapse unknown/unreported race into other
# Other is a fallback category captured very rare/sparse and unreported or 
# unknkown subgroups
b_race <- a_imp_df %>% 
  mutate(
    race = case_when(
      race == "na_unk" ~ "other_na_unk",
      race == "other" ~ "other_na_unk",
      TRUE ~ race
    ),
    race = factor(
      race, 
      levels = c("white", "asian", "black_aa", "other_na_unk")
    )
  )

# Check formats
# All formatted correctly
glimpse(b_race)
levels(b_race$race)
b_race %>% count(race)

# Identify trials with comorbidity and/or conmed exclusions
# Code taken from mltc_three_ways/Scripts/03_excludeconddruglab.R
# Run once to get exclusion lookup then commented out
# library(RSQLite)
# fromdb <- "../review_quality_decide/review_decide_mhdhlh/descisions.sqlite"
# todb <- "../review_quality_decide/review_decide_mhdhlh/backups/descisions.sqlite"
# file.copy(fromdb,
#           todb)
# rm(fromdb)
# con <- dbConnect(SQLite(), todb)
# dbListTables(con)
# decisions <- dbGetQuery(con, "SELECT * from decisions") %>% 
#   as_tibble()
# decisions <- decisions %>% 
#   arrange(desc(id)) %>% 
#   distinct(trial, .keep_all = TRUE) %>% 
#   arrange(id)
# write_csv(decisions, "Created_metadata/comorb_conmed_exclusions.csv")

# Import comorbidity and conmed exclusion lookup
# Simplify
b_imp_excl <- read_csv("Created_metadata/comorb_conmed_exclusions.csv") %>% 
  distinct(sponsor, trial, cqinc, dqinc) %>% 
  rename(
    nct_id = trial,
    comorb_excl = cqinc,
    conmed_excl = dqinc
  )

# Join actions to dataset
# Set comorbs and/or conmeds to NA for respective trials
b_excl <- b_imp_excl %>% 
  filter(nct_id %in% b_race$nct_id) %>% 
  left_join(b_race) %>% 
  group_by(nct_id) %>% 
  mutate(
    comorbs = if_else(comorb_excl == "exclude", NA, comorbs),
    conmeds = if_else(conmed_excl == "exclude", NA, conmeds)
  ) %>% 
  ungroup()

# Compare before and after
b_race %>% 
  inner_join(
    b_imp_excl %>% 
      filter(comorb_excl == "exclude" | conmed_excl == "exclude")
  ) %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(
    n_subj = n_distinct(id),
    comorb_excl = unique(comorb_excl),
    na_comorb = sum(is.na(comorbs)),
    conmed_excl = unique(conmed_excl),
    na_conmed = sum(is.na(conmeds))
  )

b_excl %>% 
  inner_join(
    b_imp_excl %>% 
      filter(comorb_excl == "exclude" | conmed_excl == "exclude")
  ) %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(
    n_subj = n_distinct(id),
    comorb_excl = unique(comorb_excl),
    na_comorb = sum(is.na(comorbs)),
    conmed_excl = unique(conmed_excl),
    na_conmed = sum(is.na(conmeds))
  )

# Fit cox models ---------------------------------------------------------------

# Renest data
# Rearrange variables so event is proceeded by all covariates
# Remove collinear variables - chol_ldl, dbp
b_mdl_df <- b_excl %>% 
  select(sponsor:weeks, category, pyr, mft, event, everything()) %>% 
  select(sponsor:trtcls, days:event, trttype, everything()) %>% 
  select(-contains("excl"), -c(chol_ldl, dbp)) %>%
  mutate(
    trttype = case_when(
      trttype == "plc" ~ 0L,
      TRUE ~ 1L
    )
  ) %>% 
  group_by(sponsor, nct_id) %>% 
  nest()

# Drop uniform or missing covariates
b_mdl_df_drp <- b_mdl_df %>% 
  mutate(
    drp = map(data, ~ {
      .x %>% 
        select(
          where(
            ~ !all(is.na(.)) &&
              !(is.factor(.x) && n_distinct(.x) <= 1) &&
              !(is.integer(.x) && n_distinct(.x) <= 1)
          )
        )
    })
  )

# Function to fit Cox model and capture warnings
fit_cox_capture_warnings <- function(data) {
  
  # Capture dataset names, identify event variable and pull covariates
  vars <- names(data)
  event_index <- which(vars == "event")
  covars <- vars[(event_index + 1):length(vars)]
  
  # Create Cox model formula
  f_str1 <- paste0("Surv(weeks, event) ~ ")
  f_str2 <- paste0(covars, collapse = " + ")
  f <- as.formula(paste0(f_str1, f_str2))
  
  # Fit model with call handlers to capture warnings
  wrn_msg <- character()
  
  fit <- tryCatch(
    withCallingHandlers(
      coxph(f, data = data),
      warning = function(w) {
        wrn_msg <<- c(wrn_msg, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  list(
    fit = if (inherits(fit, "error")) NULL else fit,
    warnings = if (length(wrn_msg) == 0) NULL else wrn_msg,
    error = if (inherits(fit, "error")) conditionMessage(fit) else NULL
  )
  
}

# Fit models
b_cox <- b_mdl_df_drp %>% 
  mutate(
    mdl_output = map(drp, fit_cox_capture_warnings),
    fit = map(mdl_output, "fit"),
    warnings = map(mdl_output, "warnings"),
    error = map(mdl_output, "error")
  )

# Check warnings
# Pull variable indices
# Get covariate names
b_cox_warnings <- b_cox %>% 
  select(sponsor, nct_id, drp, fit, warnings) %>% 
  unnest(warnings) %>% 
  mutate(
    var_index = str_extract(warnings, "variable?\\s+\\d+(?:\\s*,\\s*\\d+)*") %>% 
      str_remove("^variable?\\s+") %>% 
      str_remove_all("\\s"),
    var_names = map2(fit, var_index, ~ {
      coef_names <- names(coef(.x))
      ind <- as.integer(str_split(.y, ",", simplify = TRUE))
      out <- paste0(coef_names[ind], collapse = ", ")
    })
  ) %>% 
  unnest(var_names) %>% 
  ungroup()

# Count problematic variables
# Duration of T2DM (1 trial)
# Smoking status (1 trial)
# Race (30 trials)
b_cox_warnings %>% 
  count(var_names, sort = TRUE) %>% 
  mutate(agg = sum(n[grepl("race", var_names)]))

# Resolve convergence issues ---------------------------------------------------

# Duration of T2DM
# Very few (n=2) unique values, drop
# NCT00800683
b_cox_warnings %>% 
  filter(grepl("durn", var_names)) %>% 
  select(sponsor, nct_id, drp) %>% 
  unnest(drp) %>% 
  ungroup() %>% 
  count(durn)

# Smoke
# Only 2 smokers
# NCT00996658
b_cox_warnings %>% 
  filter(grepl("smoke", var_names)) %>% 
  select(sponsor, nct_id, drp) %>% 
  unnest(drp) %>% 
  ungroup() %>% 
  count(smoke)

# Omit durn and smoke from respective trials
c_fix1 <- b_cox %>% 
  mutate(
    drp2 = map2(nct_id, drp, ~ {
      if (.x == "NCT00800683") {
        .y %>% select(-durn)
      } else if (.x == "NCT00996658") {
        .y %>% select(-smoke)
      } else {
        .y
      }
    })
  )

# Get culprit trials
c_inspect_race <- b_cox_warnings %>% 
  filter(grepl("race", var_names)) %>% 
  mutate(var_names = gsub("race", "", var_names)) %>% 
  select(sponsor, nct_id, drp, var_names)

# Count subgroups
# In order: other_na_unk, black/aa, asian
c_inspect_race %>%
  unnest(drp) %>% 
  ungroup() %>% 
  distinct(nct_id, var_names) %>% 
  separate(var_names, into = c("r1", "r2", "r3"), sep = ", ") %>% 
  pivot_longer(r1:r3) %>% 
  filter(!is.na(value)) %>% 
  count(value, sort = TRUE)

# Subgroup frequencies across trials
# Saved for record keeping
c_inspect_race %>%
  unnest(drp) %>% 
  group_by(nct_id, race) %>% 
  reframe(
    var_names = unique(var_names),
    n = n()
  ) %>% 
  pivot_wider(names_from = "race", values_from = "n", values_fill = 0) %>% 
  # write_csv("Scratch_data/Cox_checks/race_freq_convg_issues.csv")
  print(n = Inf)

# Check events among subgroups
# Vast majority of cases all participants were censored
# Few cases where there were events but very few
c_inspect_race %>%
  unnest(drp) %>% 
  count(nct_id, var_names, race, event) %>% 
  pivot_wider(
    names_from = c("race", "event"), 
    values_from = "n",
    values_fill = 0
  ) %>% 
  select(
    nct_id, var_names,
    contains("white"), contains("asian"), contains("black"),
    everything()
  ) %>% 
  print(n = Inf)

# Solution: Collapse culprit subgroups into "other" to retain information
# Pooled effects will omit culprit subgroups (beta = 0, vcov = Inf)

# Collapse subgroups within culprit trials to "other"
c_fix_race <- c_fix1 %>% 
  select(sponsor, nct_id, drp, fit, drp2) %>% 
  left_join(c_inspect_race %>% distinct(nct_id, var_names)) %>% 
  mutate(
    drp2 = map2(var_names, drp2, ~ {
      if ("race" %in% names(.y) && !is.na(.x)) {
        tgts <- str_replace(.x, ", ", "|")
        mutate(
          .y, 
          race = if_else(grepl(tgts, race), "other_na_unk", race),
          race = factor(
            race, 
            levels = c("white", "asian", "black_aa", "other_na_unk")
          )
        )
      } else {
        .y
      }
    })
  )

# Verify changes were implemented correctly
c_fix_race %>% 
  filter(!is.na(var_names)) %>% 
  unnest(drp2) %>% 
  group_by(nct_id, race) %>% 
  reframe(
    var_names = unique(var_names),
    n = n()
  ) %>% 
  pivot_wider(names_from = "race", values_from = "n", values_fill = 0) %>% 
  # write_csv("Scratch_data/Cox_checks/race_freq_convg_issues_fixed.csv")
  print(n = Inf)

# Refit models
c_cox <- c_fix_race %>% 
  select(-var_names) %>% 
  mutate(
    mdl_output = map(drp2, fit_cox_capture_warnings),
    fit2 = map(mdl_output, "fit"),
    warnings = map(mdl_output, "warnings"),
    error = map(mdl_output, "error")
  )

# Check remaining warnings
# 17 trials, predominantly caused by other_na_unk
# 1 trial black_aa
c_cox_warnings <- c_cox %>% 
  select(sponsor, nct_id, drp, fit2, warnings) %>% 
  unnest(warnings) %>% 
  mutate(
    var_index = str_extract(warnings, "variable?\\s+\\d+(?:\\s*,\\s*\\d+)*") %>% 
      str_remove("^variable?\\s+") %>% 
      str_remove_all("\\s"),
    var_names = map2(fit2, var_index, ~ {
      coef_names <- names(coef(.x))
      ind <- as.integer(str_split(.y, ",", simplify = TRUE))
      out <- paste0(coef_names[ind], collapse = ", ")
    })
  ) %>% 
  unnest(var_names) %>% 
  ungroup()

c_cox_warnings %>% count(var_names)

# Inspect trials - no events in culprit subgroups
d_inspect_race <- c_cox_warnings %>% 
  filter(grepl("race", var_names)) %>% 
  mutate(var_names = gsub("race", "", var_names)) %>% 
  select(sponsor, nct_id, drp, var_names)

d_inspect_race %>%
  unnest(drp) %>% 
  count(nct_id, var_names, race, event) %>% 
  pivot_wider(
    names_from = c("race", "event"), 
    values_from = "n",
    values_fill = 0
  ) %>% 
  select(
    nct_id, var_names,
    contains("white"), contains("asian"), contains("black"),
    everything()
  ) %>% 
  print(n = Inf)

# Solution: drop participants in subgroups
# Dropping participants risks biasing estimates
d_cox <- c_cox %>% 
  select(sponsor, nct_id, drp, fit, drp2, fit2) %>% 
  left_join(d_inspect_race %>% distinct(nct_id, var_names)) %>% 
  mutate(
    drp3 = pmap(
      list(nct_id, drp2, var_names),
      function(nct_id, drp2, var_names) {
        if (nct_id %in% d_inspect_race$nct_id && !is.na(var_names)) {
          drp2 %>% filter(!race %in% str_split(var_names, ",\\s*")[[1]])
          # .y %>% select(-race)
        } else {
          drp2
        }
      }
    ),
    mdl_output = map(drp3, fit_cox_capture_warnings),
    fit3 = map(mdl_output, "fit"),
    warnings = map(mdl_output, "warnings"),
    error = map(mdl_output, "error")
  )

# No warnings remain (and no errors)
d_cox %>% filter(map_lgl(warnings, ~ !is.null(.x)))
d_cox %>% filter(map_lgl(error, ~ !is.null(.x)))

# Get outputs ------------------------------------------------------------------

# Results, diagnostics, and vcov
d_out <- d_cox %>% 
  arrange(sponsor, nct_id) %>% 
  select(sponsor, nct_id, drp3, fit3) %>% 
  mutate(
    res = map(fit3, tidy),
    diag = map(fit3, glance),
    vcov = map(fit3, vcov),
    ph = map(fit3, ~ cox.zph(.x, term = FALSE)),
    ph_df = map(ph, ~ as_tibble(.x$table, rownames = "term"))
  )
  
# Save
saveRDS(d_out, "Processed_data/06_cox_out.rds")


