
# Prepare IPD baseline characteristics for generating participant profiles
# Borrows extracts from previous network meta-analysis, need to have this repo
# pulled locally in order to run this script to link with filepaths
# Does not have comorbidities and conmeds as this was constructed specifically
# for the conduct-stage model, this was extracted from Vivli for project

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Filepath to nma repo
nma_path <- "../nma_agesex_public/"
nma_path2 <- "../nma_more_covariates/"

# Import AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import comorbidity and conmed summaries
imp_count_sum <- read_csv("vivli/Summaries/overall_summary.csv") %>%
  filter(name %in% c("comorbs", "conmeds"))

# Import Vivli arm-level summary
imp_arm <- read_csv("vivli/Summaries/trial_level_summary.csv")

# Get covariate names
imp_covars <- read_csv("vivli/Models/export_res.csv") %>% distinct(term)

# Numeric covariates except age
# Distribution-dependent, kept nested to obtain mean and d
files_bl <- list.files(
  paste0(nma_path2, "Data/"),
  pattern = "covariates",
  full.names = TRUE
)

imp_bl <- tibble(
  dist = gsub("covariates_|\\.csv", "", basename(files_bl)),
  imp = map(files_bl, read_csv, show_col_types = FALSE)
)

# Age and sex distributions
imp_bl2 <- read_csv(paste0(nma_path2, "Data/ipd_age_sex.csv"))

# Categorical covariates
imp_bl3 <- read_csv(paste0(nma_path2, "Data/categorical_variables_percentage.csv"))

# HbA1c
imp_bl4 <- read_csv(paste0(nma_path2, "Data/hba1c_base_summary_ipd.csv"))

# Import simulated IPD from previous repo to get arm counts
imp_n <- readRDS(paste0(nma_path, "Scratch_data/agg_ipd_hba1c.Rds")) %>%
  select(ipd) %>%
  unnest(ipd) %>%
  ungroup() %>%
  group_by(nct_id, arm_f, arm_lvl, trtcls5) %>%
  reframe(total = n())

# Prepare numeric covariates -------------------------------------------------

# Determine the appropriate distribution to pass to add_integration() for
# each covariate - based on the distribution that most frequently best
# described a covariate
imp_bl %>%
  unnest(imp) %>%
  ungroup() %>%
  distinct(dist, nct_id, var) %>%
  count(dist, var) %>%
  group_by(var) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  write_csv("scratch_data/sobol_covariate_dist.csv")

# Get mean and sd for numeric covariates
# Converted depending on the distribution that best described covariates
bl_num <- imp_bl %>%
  mutate(
    est = map2(
      imp,
      dist,
      ~ {

        if (.y == "gamma") {

          .x %>%
            mutate(
              mu = est_shape / est_rate,
              d = sqrt(est_shape) / est_rate
            )

        } else if (.y == "lnorm") {

          .x %>%
            mutate(
              mu = exp(est_meanlog + 0.5 * est_sdlog^2),
              d = sqrt(
                (exp(est_sdlog^2) - 1) * exp(2 * est_meanlog + est_sdlog^2)
              )
            )

        } else if (.y == "norm") {

          .x %>%
            mutate(
              mu = est_mean,
              d = est_sd
            )

        } else if (.y == "weibull") {

          .x %>%
            mutate(
              mu = est_scale * gamma(1 + 1 / est_shape),
              d = est_scale * sqrt(
                gamma(1 + 2 / est_shape) - gamma(1 + 1 / est_shape)^2
              )
            )

        }

      }
    )
  ) %>%
  select(dist, est) %>%
  unnest(est) %>%
  ungroup() %>%
  select(dist, nct_id, arm_f, var, mu, d)

# Prepare age distribution
# Trial-level so assumed to be applicable to arms within trials
bl_age <- bl_num %>%
  distinct(nct_id, arm_f) %>%
  left_join(
    imp_bl2 %>%
      select(nct_id, age_m, age_s) %>%
      mutate(var = "age", dist = "norm") %>%
      rename(mu = age_m, d = age_s)
  )

# Combine numeric distributions
bl_num_all <- bl_num %>%
  full_join(bl_age) %>%
  arrange(nct_id, arm_f, var)

# No missingness
check_nas(bl_num_all)

# Join arm counts from simulated IPD
bl_num_n <- bl_num_all %>%
  mutate(arm_f = gsub("_dcmplx", "", arm_f)) %>%
  left_join(
    imp_n %>%
      mutate(
        arm_f_old = arm_f,
        arm_f = paste0(
          substr(arm_f_old, 2, 3),
          arm_lvl
        )
      ) %>%
      select(-arm_f_old) %>%
      distinct()
  )

check_nas(bl_num_n)

# Get AttRisk IPD trials
bl_num_set <- bl_num_n %>% inner_join(imp_arm %>% distinct(nct_id))
n_distinct(bl_num_set$nct_id)
check_nas(bl_num_set)

# Drop experimental arms
# Done by harmonising covariate distribution arm variable to match
# AttRisk arm-level attrition counts metadata, then inner-join to drop arms
# Manually add arm count for NCT00800683 placebo
bl_num_rm_exp <- bl_num_set %>%
  mutate(trtcls5 = substr(arm_f, 3, 7)) %>%
  inner_join(imp_arm %>% distinct(nct_id, trtcls5, trttype, n)) %>%
  select(dist, nct_id, trtcls5, trttype, total, var, mu, d) %>%
  distinct() %>%
  mutate(
    total = if_else(is.na(total) & nct_id == "NCT00800683", 63, total)
  )

# Compare number of arms per trial to verify join was done correctly
# 83 arms in each, join was done correctly
bl_num_set %>%
  mutate(trtcls5 = substr(arm_f, 3, 7)) %>%
  distinct(nct_id, trtcls5) %>%
  inner_join(imp_arm %>% distinct(nct_id, trtcls5))

imp_arm %>% distinct(nct_id, trtcls5)

# Address arms with >1 distribution (arms of different doses)
# These were collapsed for AttRisk
bl_num_rm_exp %>%
  group_by(nct_id, trtcls5, var) %>%
  filter(n() > 1)

# Join age
bl_num_age <- bl_num_rm_exp %>%
  full_join(
    bl_num_rm_exp %>%
      distinct(nct_id, trtcls5, trttype) %>%
      left_join(
        imp_bl2 %>%
          mutate(var = "age", dist = "norm") %>%
          rename(mu = age_m, d = age_s) %>%
          select(-male_p)
      )
  )

# Join HbA1c
bl_num_hba1c <- imp_bl4 %>%
  rename(arm_lvl = arm_f) %>%
  left_join(
    imp_n %>%
      mutate(
        arm_f_old = arm_f,
        arm_f = paste0(
          substr(arm_f_old, 2, 3),
          arm_lvl
        )
      ) %>%
      select(-arm_f_old) %>%
      distinct()
  ) %>%
  mutate(
    var = "hba1c",
    dist = "norm"
  ) %>%
  inner_join(
    bl_num_age %>% distinct(nct_id, trtcls5, total)
  )

bl_num_hba1c <- bl_num_age %>%
  full_join(
    bl_num_hba1c %>%
      rename(mu = hba1cbase_m, d = hba1cbase_s) %>%
      select(-c(arm_f, arm_lvl))
  ) %>%
  mutate(
    trttype = case_when(
      is.na(trttype) & trtcls5 == "place" ~ "plc",
      is.na(trttype) & trtcls5 != "place" ~ "ac",
      TRUE ~ trttype
    )
  )

# Multiple arms exist in some trials (e.g. metformin of different doses)
# Arms were collapsed for AttRisk set
# Collapsing distribution summary across these trials
# Use sample size-weighting to collapse
bl_num_clp <- bl_num_hba1c %>%
  group_by(nct_id) %>%
  mutate(
    trttype = case_when(
      n_distinct(trttype) == 1 ~ trttype,
      TRUE ~ "plc_ac"
    )
  ) %>%
  ungroup() %>%
  group_by(dist, nct_id, trttype, var) %>%
  reframe(
    arm_total = sum(total),
    mu_pool = sum(total * mu) / arm_total,
    var_pool = sum((total - 1) * d^2 + total * (mu - mu_pool)^2) / (arm_total - 1),
    d_pool = sqrt(var_pool)
  ) %>%
  select(dist, nct_id, trttype, var, mu_pool, d_pool)

# Join AttRisk arm size
bl_num_rdy <- bl_num_clp %>%
  left_join(
    imp_arm %>%
      distinct(nct_id, trtcls5, trttype, n) %>%
      group_by(nct_id) %>%
      mutate(
        trttype = case_when(
          n_distinct(trttype) == 1 ~ trttype,
          TRUE ~ "plc_ac"
        )
      ) %>%
      group_by(nct_id, trttype) %>%
      reframe(n = sum(n))
  ) %>%
  arrange(nct_id, var)

# No missingness
check_nas(bl_num_rdy)

# Prepare categorical covariates ---------------------------------------------

# Collapse other race categories into one
# Convert percentages to proportions
bl_cat <- imp_bl3 %>%
  mutate(
    across(-c(nct_id, arm_f, cursmoke), ~ if_else(is.na(.), 0, .)),
    other = ai_an + other
  ) %>%
  select(
    nct_id,
    arm_f,
    cursmoke,
    white,
    black_aa,
    asian,
    other
  ) %>%
  mutate(
    across(-c(nct_id, arm_f), ~ . / 100)
  )

# Verify that race categories sum to 100%
# All do, rounding errors expected (either 0.99 or 1.01)
bl_cat %>%
  mutate(race_total = white + black_aa + asian + other) %>%
  filter(race_total != 1) %>%
  print(n = Inf)

# Join groups for merging with AttRisk metadata
bl_cat_grps <- bl_cat %>%
  left_join(
    bl_num_n %>% distinct(nct_id, arm_f, trtcls5)
  ) %>%
  mutate(
    trtcls5 = if_else(
      is.na(trtcls5)
      & nct_id == "NCT00800683"
      & arm_f == "1_placebo",
      "place",
      trtcls5
    )
  )

# Get AttRisk set
bl_cat_set <- bl_cat_grps %>%
  inner_join(
    imp_arm %>% distinct(nct_id, trtcls5, trttype)
  )

# Add sex proportion
bl_cat_rdy <- bl_cat_set %>%
  left_join(
    bl_cat_set %>%
      distinct(nct_id, trtcls5, trttype) %>%
      left_join(imp_bl2 %>% select(nct_id, male_p))
  )

# Collapse percentages across trttype within trials
bl_cat_clp <- bl_cat_rdy %>%
  group_by(nct_id) %>%
  mutate(
    trttype = case_when(
      n_distinct(trttype) == 1 ~ trttype,
      TRUE ~ "plc_ac"
    )
  ) %>%
  left_join(bl_num_rdy %>% distinct(nct_id, trttype, n)) %>%
  group_by(nct_id, trttype) %>%
  reframe(
    across(
      c(male_p, cursmoke, white, black_aa, asian, other),
      ~ weighted.mean(.x, w = n, na.rm = TRUE)
    )
  )

# Current smoker missing in 7 trials
check_nas(bl_cat_clp)

# Combine covariate distribution summaries and tidy --------------------------

# Widen numeric covariates summaries, one row per IPD trial arm
cmbn_num_wide <- bl_num_rdy %>%
  rename(mu = mu_pool, d = d_pool) %>%
  select(-dist) %>%
  pivot_wider(
    names_from = "var",
    values_from = c("mu", "d"),
    names_glue = "{var}_{.value}"
  )

# Combine
# Swap male p for female p
# Drop cholesterol ratio, not in AttRisk covariates
cmbn_all <- cmbn_num_wide %>%
  left_join(bl_cat_clp) %>%
  mutate(female_p = 1 - male_p) %>%
  select(-c(male_p, contains("tot_hdl")))

# Confirm all AttRisk covariates are present
# HbA1c, comorbs and conmeds need added
imp_covars
glimpse(cmbn_all)

# Check missingness
# Missingness in durn, cholesterol, egfr, cursmoke, bmi, sbp
# Mostly durn
check_nas(cmbn_all) %>% arrange(desc(value)) %>% print(n = Inf)

# Define covariate sets
# Variables with mean/SD
num_vars <- c(
  "age",
  "bmi",
  "chol_hdl",
  "chol_total",
  "hba1c",
  "durn",
  "egfr",
  "sbp"
)

# Variables as arm-level proportions
prop_vars <- c(
  "cursmoke",
  "female_p",
  "white",
  "black_aa",
  "asian",
  "other"
)

# Save for manuscript summaries
write_csv(cmbn_all, "outputs/conduct_stage/ipd_bl_sum.csv")
