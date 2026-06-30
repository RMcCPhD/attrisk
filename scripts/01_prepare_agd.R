
# Prepare AgD - trial design features, eligibility criteria
# Processing is logged, with rationales, in "logs/data_processing.xlsx"
# Starting with the differential attrition analysis set, same eligiblity apply
# Included 364 trials: 272 with AgD, 92 with IPD:
# - Prespecified follow-up duration
# - At least 125 randomised participants
# - Registered on clinicaltrials.gov

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import differential attrition trial set
imp_dat <- readRDS("data/tidied_agg.rds")

# Import AgD trial characteristics
z_agd_chars <- "scratch_data/data_processing/aact_manual_found_missingness.csv"
imp_agd_chars <- read_csv(z_agd_chars)

# Import AgD eligibility criteria
imp_agd_elig <- read_csv("data/elig.csv")

# Import IPD model outputs
imp_vivli_models <- tibble(
  file = list.files("vivli/v2/Models/",full.names = TRUE),
  imp = map(file, read_csv)
)

# Import IPD arm-level summaries
imp_vivli_arms <- tibble(
  file = list.files("vivli/v2/Summaries/", full.names = TRUE),
  imp = map(file, read_csv)
)

# Import antidiabetic ATC lookup
imp_atc <- read_csv("created_metadata/atc.csv")

# Trial design -----------------------------------------------------------------

# Summarise trial size
# Keep only non-experimental arms
# Drops trials that compared newer antidiabetics to other newer antidiabetics
# E.g. a GLP1 vs DPP4 trial
agd_rm_exp_arms <- imp_dat %>% 
  group_by(trial_id) %>% 
  mutate(n = sum(arm_n)) %>% 
  ungroup() %>% 
  filter(!class %in% c("glp1", "dpp4", "sglt2"))

# Compare trials and sources before/after
imp_dat %>% count(trial_id, source) %>% count(source)
agd_rm_exp_arms %>% count(trial_id, source) %>% count(source)

# Add treatment type
# Where a single arm (i.e. placebo or active), label appropriately
# If multiple arms, labelled to reflect this
agd_trttype <- agd_rm_exp_arms %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      n() > 1 ~ 2L,
      n() == 1 & class != "placebo" ~ 1L,
      n() == 1 & class == "placebo" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>% 
  ungroup()

# Verify - sums to 342 (one group per trial)
agd_trttype %>% count(trial_id, trttype) %>% count(trttype)

# Collapse AgD arms into single-arm trials by treatment type
# Each trial then has one row
agd_collapse <- agd_trttype %>% 
  group_by(source, trial_id, trttype, n) %>% 
  reframe(across(c(arm_n, arm_attr), ~ sum(.)))

# Join trial characteristics
# Take distinct rows and slice first row (one trial has differing value for
# south-east asian region, taking first where one site was indicated to be here)
agd_chars <- agd_collapse %>% 
  left_join(
    imp_agd_chars %>% 
      select(-c(n, arm_n, arm_n_attr, type)) %>% 
      distinct() %>% 
      group_by(trial_id) %>% 
      slice(1) %>% 
      ungroup()
  )

# Double-check missingness
# 7 trials are missing several trial design features:
# - Number of arms, primary outcomes, and WHO regions
# - Duration, mand masking
# 23 trials are missing number of facilities
check_nas(agd_chars) %>% print(n = Inf)

# Save current version of AgD data to find remaining missing characteristics
# Commented out to prevent accidental overwrite (renamed anyway)
# write_csv(
#   agd_chars, 
#   "scratch_data/data_processing/find_missing_agd_chars.csv"
# )

# Import AgD data with findings for missing characteristics
# 11 remain that are missing number of facilities
z_agd_found <- "scratch_data/data_processing/found_missing_agd_chars.csv"
agd_found_chars <- read_csv(z_agd_found)
check_nas(agd_found_chars) %>% print(n = Inf)

# Imputing n_facilities where missing for 11 trials
# Rationale:
# - Small amount of missingness
# - All other characteristics are available and complete
# - Plausible that this is predictable from geographic spread and complexities

# Vector of numeric trial characteristics
z_numeric_agd_vars <- c(
  "arm_n",
  "arm_attr",
  "number_of_arms",
  "duration",
  "n_poms",
  "n_regions",
  "region_of_the_americas",
  "european_region",
  "western_pacific_region",
  "african_region",
  "south_east_asia_region",
  "eastern_mediterranean_region"
)

# Check distribution of n_facilities
# Skewed as expected, more regular when log-transformed
agd_found_chars %>% ggplot(aes(x = n_facilities)) + geom_density()
agd_found_chars %>% ggplot(aes(x = log1p(n_facilities))) + geom_density()

# Check correlation with other numeric predictors
# Similar strengths for normal and log-transformed n_facilities
# Strongest for trial size, regions, arm size, attrition, and duration
# Meaningfully related to other numeric covariates
agd_found_chars %>% 
  select(n_facilities, where(~is.numeric(.))) %>% 
  reframe(
    across(
      c(all_of(z_numeric_agd_vars), n),
      ~ cor(n_facilities, .x, method = "spearman", use = "complete.obs")
    )
  ) %>% 
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "normal_cor"
  ) %>%
  arrange(desc(abs(normal_cor))) %>% 
  left_join(
    agd_found_chars %>% 
      mutate(
        log_facilities = log1p(n_facilities),
        n = log1p(n)
      ) %>% 
      select(log_facilities, where(~is.numeric(.))) %>% 
      select(-n_facilities) %>% 
      reframe(
        across(
          c(all_of(z_numeric_agd_vars), n),
          ~ cor(log_facilities, .x, method = "pearson", use = "complete.obs")
        )
      ) %>% 
      pivot_longer(
        everything(),
        names_to = "variable",
        values_to = "log_cor"
      ) %>%
      arrange(desc(abs(log_cor)))
  )

# Impute missing n_facilities with a simple model
# Based on moderately correlated trial design and geographic features
# Balances use of information; avoiding risk of overfitting or enforcing
# an overly strict relationship between variables
# Excluding outcome-derived variables like attrition to avoid leakage
# Variables: log1p(n), americas, europe, african, duration, number of arms
agd_impute_fit <- lm(
  log1p(n_facilities) ~ 
    log1p(n) + 
    region_of_the_americas + 
    european_region + 
    african_region +
    duration + 
    number_of_arms,
  data = agd_found_chars %>% filter(!is.na(n_facilities))
)

# Predict values for missing trials
agd_predict_na <- agd_found_chars %>% 
  filter(is.na(n_facilities)) %>% 
  mutate(
    n_facilities = round(expm1(predict(agd_impute_fit, newdata = .)))
  )

# Compare distribution before and after imputation
# Highly similar
agd_plot_imp <- agd_found_chars %>% 
  filter(!is.na(n_facilities)) %>% 
  select(n_facilities) %>% 
  mutate(type = "before") %>% 
  full_join(
    agd_found_chars %>% 
      filter(!is.na(n_facilities)) %>% 
      full_join(agd_predict_na) %>% 
      mutate(type = "after")
  ) %>% 
  mutate(
    type = factor(
      type,
      levels = c("before", "after"),
      labels = c("Before", "After")
    )
  ) %>% 
  ggplot(aes(x = n_facilities, colour = type)) +
  geom_density() +
  theme_classic() +
  labs(x = "Number of sites", y = "Density", colour = "Imputation stage")

ggsave(
  "outputs/design/dist_before_after_sites_imputation.png",
  agd_plot_imp,
  width = 6,
  height = 4,
  units = "in"
)

# Incorporate imputed n_facilities
# Same reduction applied as when joining tidying characteristics as this was
# spotted after dealing with missing characteristics, so needed to be applied
# twice
agd_complete <- agd_found_chars %>% 
  filter(!is.na(n_facilities)) %>% 
  full_join(agd_predict_na) %>% 
  arrange(trial_id) %>% 
  distinct() %>% 
  group_by(trial_id) %>% 
  slice(1) %>% 
  ungroup()

# Verify missingness
# All trials retained
check_nas(agd_complete) %>% print(n = Inf)
agd_complete %>% count(trial_id, source) %>% count(source)

# Eligibility criteria ---------------------------------------------------------

# Get eligibility criteria for analysis set
# All present
elig_set <- imp_agd_elig %>% inner_join(agd_complete %>% distinct(trial_id))
elig_set %>% count(trial_id)

# Check missingness
# Quite a lot but much is simply the absence of a restriction (e.g. upper age)
check_nas(elig_set) %>% arrange(desc(value)) %>% print(n = Inf)

# Eligibility criteria included:
# Age limits, BMI limits, eGFR limits
# Chosen as they are core determinants of T2DM populations, absence was
# interpretable as no explicit eligibility restriction for the characteristics
# Not 100% based in reality as always potential for absence of information,
# but greater confidence in that being the case for these eligibility criteria
# compared to others (e.g. systolic blood pressure, glycaemia)
elig_keep <- elig_set %>%
  select(
    trial_id, contains(c("age", "bmi", "renal"))
  )

# Inspect units - years, kg/m2, ml/min/1.73m2 or mL/min (standard)
# Kidney function not always eGFR
elig_keep %>%
  count(renal_function_unit, renal_function_measure)
  # count(bmi_unit)
  # count(minimum_age_unit, maximum_age_unit)

# Inspect criteria
# - Age, minimum either 18 or above; maximum either 80 or above
elig_keep %>% ggplot(aes(x = minimum_age)) + geom_density()
elig_keep %>% ggplot(aes(x = maximum_age)) + geom_density()

# - BMI, mix of lower and upper limits in category
elig_keep %>% count(bmi)
elig_keep %>% count(less_than_bmi_excluded)

# - Kidney function, not always eGFR uncorrected but most are (where not NA)
elig_keep %>% 
  count(
    renal_function_measure, 
    renal_function_unit, 
    renal_function_range
  ) %>% 
  filter(grepl("m2", renal_function_unit)) %>% 
  arrange(renal_function_range) %>% print(n = Inf)

# Age straightforward index
elig_age <- elig_keep %>%
  select(trial_id, contains("age")) %>% 
  mutate(
    lower_age = case_when(
      minimum_age > 18 ~ ">18",
      TRUE ~ "18"
    ),
    upper_age = case_when(
      maximum_age <= 80 ~ "80",
      TRUE ~ ">80"
    )
  )

# Verify indexing
elig_age %>% count(maximum_age, upper_age)
elig_age %>% count(minimum_age, lower_age)
elig_age %>% count(lower_age, upper_age)

# BMI
# Isolate upper and lower limits
# Convert numeric to categories
elig_bmi <- elig_keep %>% 
  select(trial_id, contains("bmi")) %>% 
  mutate(
    lower_bmi = case_when(
      !is.na(less_than_bmi_excluded) ~ as.character(less_than_bmi_excluded),
      !is.na(bmi) & grepl("<", bmi) ~ bmi,
      TRUE ~ NA_character_
    ),
    lower_bmi = case_when(
      !is.na(lower_bmi) ~ "<27",
      TRUE ~ "none"
    ),
    upper_bmi = case_when(
      !is.na(more_than_bmi_excluded) ~ as.character(more_than_bmi_excluded),
      !is.na(bmi) & grepl(">", bmi) ~ bmi,
      TRUE ~ NA_character_
    ),
    upper_bmi = case_when(
      grepl("40|45|50|60", upper_bmi) ~ ">40",
      grepl("30|35", upper_bmi) ~ ">30",
      TRUE ~ "none"
    )
  )

# Verify indexing
elig_bmi %>% count(bmi, more_than_bmi_excluded, upper_bmi)
elig_bmi %>% count(bmi, less_than_bmi_excluded, lower_bmi)
elig_bmi %>% count(upper_bmi, lower_bmi)

# eGFR straightforward
elig_egfr <- elig_keep %>% 
  select(trial_id, contains("renal")) %>% 
  mutate(
    lower_egfr = case_when(
      grepl("1.73", renal_function_unit)
      & grepl("gfr", renal_function_measure, ignore.case = TRUE)
      & grepl("15|30", renal_function_range) ~ "<30",
      grepl("1.73", renal_function_unit)
      & grepl("gfr", renal_function_measure, ignore.case = TRUE)
      & grepl("35|45|50|55|60", renal_function_range) ~ "<60",
      TRUE ~ "none"
    )
  )

# Verify indexing
elig_egfr %>% count(renal_function_unit, renal_function_range, lower_egfr)
elig_egfr %>% count(lower_egfr)

# Join indexed eligibility criteria
elig_join <- reduce(mget(ls(pattern = "elig_(age|bmi|egfr)")), left_join) %>% 
  select(trial_id, contains(c("upper", "lower"))) %>% 
  select(trial_id, contains("age"), contains("bmi"), contains("egfr")) %>% 
  distinct()

# Inspect
# More rows than trials
elig_join %>% 
  reframe(
    n_rows = n(),
    n_trials = n_distinct(trial_id)
  )

# Duplicate rows apart from one trial, collapse (had multirow eGFR eligibility)
elig_join %>% group_by(trial_id) %>% filter(n() > 1)

elig_ready <- elig_join %>% 
  mutate(lower_egfr = if_else(trial_id == "NCT02577016", "<60", lower_egfr)) %>% 
  distinct()

# Double-check trial with >2 unique rows after collapse
elig_ready %>% filter(trial_id == "NCT02577016")

# Combine and tidy -------------------------------------------------------------

# Join eligibility criteria with design features
# Check missingness
tidy_join <- agd_complete %>% left_join(elig_ready) %>% distinct()
check_nas(tidy_join) %>% print(n = Inf)
glimpse(tidy_join)

# Inspect characteristics to inform formatting
# - Trial size: Skewed, majority <1000
# - Arm size: Skewed, majority <300
tidy_join %>% ggplot(aes(x = n)) + geom_density()
tidy_join %>% ggplot(aes(x = arm_n)) + geom_density()
# - Phase: Majority 3, 2 are labelled "Not applicable"
tidy_join %>% count(phase)
# - Number of arms: Vast majority 2 or 3
tidy_join %>% count(number_of_arms)
# - Number of facilities: Vast majority <200
tidy_join %>% ggplot(aes(x = n_facilities)) + geom_density()
# - Number of POMs: Vast majority have 1
tidy_join %>% count(n_poms)
# - Masking: Vast majority double-blinded
tidy_join %>% count(masking)
# - WHO regions: Flag as either present or absent

# Excluded variables
# - arm_attr: Outcome
# - duration: Parameter of link function for predicting outcome
# - n_regions: Collinear with individual regions

# Convert covariates to binary or categorical variables
# See "Logs/data_processing.xlsx" sheet 3
tidy_convert <- tidy_join %>%
  select(-n_regions) %>%
  mutate(
    n_cat = case_when(
      n < 250 ~ "<250",
      n >= 250 & n <= 499 ~ "250-499",
      n >= 500 & n <= 749 ~ "500-749",
      n >= 750 & n <= 999 ~ "750-999",
      TRUE ~ ">=1000"
    ),
    arm_n_cat = case_when(
      arm_n < 100 ~ "<100",
      arm_n >= 100 & arm_n <= 199 ~ "100-199",
      arm_n >= 200 & arm_n <= 299 ~ "200-299",
      arm_n >= 300 ~ ">=300"
    ),
    phase = case_when(
      grepl("3", phase) ~ 0L,
      grepl("4", phase) ~ 1L,
      TRUE ~ 0L
    ),
    n_arms_cat = case_when(
      number_of_arms == 2 ~ 0L,
      TRUE ~ 1L
    ),
    n_sites_cat = case_when(
      n_facilities < 100 ~ "<100",
      n_facilities >= 100 & n_facilities <= 199 ~ "100-199",
      n_facilities >= 200 ~ ">=200"
    ),
    n_poms_cat = case_when(
      n_poms == 1 ~ 0L,
      TRUE ~ 1L
    ),
    masking = case_when(
      grepl("open", masking, ignore.case = TRUE) ~ "open_label",
      grepl("double", masking, ignore.case = TRUE) ~ "double",
      grepl("triple", masking, ignore.case = TRUE) ~ "triple",
      grepl("quad", masking, ignore.case = TRUE) ~ "quadruple"
    ),
    across(contains("region"), ~ if_else(. == 0, 0L, 1L))
  )

# Format covariates
# Convert duration to years
tidy_format <- tidy_convert %>%
  mutate(
    trttype = factor(
      trttype,
      levels = c(0, 1, 2),
      labels = c("placebo", "active", "placebo_active")
    ),
    n_cat = factor(
      n_cat,
      levels = c("<250", "250-499", "500-749", "750-999", ">=1000")
    ),
    arm_n_cat = factor(
      arm_n_cat,
      levels = c("<100", "100-199", "200-299", ">=300")
    ),
    n_sites_cat = factor(
      n_sites_cat,
      levels = c("<100", "100-199", ">=200")
    ),
    masking = factor(
      masking,
      levels = c("open_label", "double", "triple", "quadruple")
    ),
    upper_age = factor(upper_age, levels = c("80", ">80")),
    lower_age = factor(lower_age, levels = c("18", ">18")),
    upper_bmi = factor(upper_bmi, levels = c("none", ">30", ">40")),
    lower_bmi = factor(lower_bmi, levels = c("none", "<27")),
    lower_egfr = factor(lower_egfr, levels = c("none", "<30", "<60")),
    duration = duration / 52
  ) %>%
  select(-c(n, number_of_arms, n_facilities, n_poms))

glimpse(tidy_format)

# Arrange variables
tidy_arrange <- tidy_format %>%
  select(
    source, trial_id, arm_n, arm_attr, duration,
    trttype, contains("cat"), everything()
  )

# Check missingness
check_nas(tidy_arrange) %>% print(n = Inf)

# Create model dataset ---------------------------------------------------------

# Helper: Check model covariates after converting and formatting above
mdl_check_vars <- function(data, vars, stop_on_na = FALSE) {
  
  dat <- data %>% select(all_of(vars))
  
  # 1. Missingness / structure summary
  na_check <- map_dfr(names(dat), function(v) {
    x <- dat[[v]]
    
    tibble(
      variable = v,
      class = paste(class(x), collapse = "/"),
      n_rows = length(x),
      n_na = sum(is.na(x)),
      pct_na = round(mean(is.na(x)) * 100, 2),
      n_unique_non_missing = n_distinct(x, na.rm = TRUE)
    )
  })
  
  # 2. Unique values / levels and counts
  value_check <- map_dfr(names(dat), function(v) {
    x <- dat[[v]]
    
    tibble(value = x) %>%
      count(value, name = "n", .drop = FALSE) %>%
      mutate(
        variable = v,
        value = as.character(value),
        value = if_else(is.na(value), "<NA>", value),
        pct = round(n / sum(n) * 100, 2)
      ) %>%
      select(variable, value, n, pct)
  })
  
  # 3. Variables with no variation
  no_variation <- na_check %>% filter(n_unique_non_missing <= 1)
  
  cat("\n--- Missingness / structure check ---\n")
  print(na_check, n = Inf)
  
  cat("\n--- Unique values / factor levels / integer counts ---\n")
  print(value_check, n = Inf)
  
  if (nrow(no_variation) > 0) {
    cat("\n--- Variables with <=1 non-missing unique value ---\n")
    print(no_variation, n = Inf)
  }
  
  if (stop_on_na && any(na_check$n_na > 0)) {
    stop("Missing values detected. Inspect `na_check` before generating model.matrix().")
  }
  
  invisible(
    list(
      na_check = na_check,
      value_check = value_check,
      no_variation = no_variation
    )
  )
  
}

# All looks good
mdl_vars <- names(tidy_arrange[, -c(1:5)])
mdl_data_chk <- mdl_check_vars(
  data = tidy_arrange,
  vars = mdl_vars,
  stop_on_na = TRUE
)

# Use model.matrix to create model dataset
mdl_data_expr <- paste0("~ ", paste0(mdl_vars, collapse = " + "))
mdl_data_fml <- as.formula(mdl_data_expr)
mdl_data_mat <- model.matrix(mdl_data_fml, data = tidy_arrange)[, -1]

# Outputs ----------------------------------------------------------------------

# Created nested data for each stage of AgD processing
# 342 rows in each, all trials accounted for
out_dat <- tibble(
  agd_arms_chars = list(agd_complete),
  elig = list(elig_ready),
  format_covars = list(tidy_arrange)
)

# Save outputs
saveRDS(out_dat, "processed_data/agd_steps_out.rds")
saveRDS(mdl_data_mat, "processed_data/agd_covar_mat.rds")

