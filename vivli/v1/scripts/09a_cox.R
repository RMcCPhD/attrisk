
# Fit individual-level cox models with imputed dataset
# Save outputs
source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(broom)
library(survival)
library(mfp)

# Import datasets
a_imp_df <- readRDS("Processed_data/n77_cox.rds")
a_imp_hier <- readRDS("Processed_data/imputed_variables.rds")

# Prepare ----------------------------------------------------------------------

# Harmonise race, format trttype
b_df <- a_imp_df %>% 
  distinct() %>% 
  select(nct_id, id, trttype, days:category) %>% 
  left_join(a_imp_hier %>% distinct()) %>% 
  mutate(
    race = case_when(
      race %in% c("ai_an", "multiple", "nh_pi") ~ "other",
      race == "missing_redacted_unknown" ~ NA,
      TRUE ~ race
    ),
    trttype = factor(trttype, levels = c("plc", "ac"))
  ) %>% 
  select(-c(dbp, chol_ldl))

# Prepare lookup
z_rm_race <- read_csv("Scratch_data/rm_race.csv") %>% 
  pivot_longer(race:race3, names_to = "nm", values_to = "race") %>% 
  filter(!is.na(race)) %>% 
  select(-nm) %>% 
  group_by(nct_id) %>% 
  mutate(race = paste0(unique(race), collapse = ", ")) %>% 
  distinct() %>% 
  mutate(race_to_exclude = strsplit(race, ", ")) %>% 
  select(-race)

# Remove race in respective trials by matching to levels in lookup lists
# 05.08.25 - Additional removal in two trials, convergence
# issues after fixing comorbidity extraction
b_df_race <- b_df %>% 
  left_join(z_rm_race) %>% 
  mutate(
    exclude = map2_lgl(
      race, race_to_exclude,
      ~ if (is.null(.y)) {FALSE} else {.x %in% .y}
    )
  ) %>% 
  filter(
    !exclude,
    !(nct_id == "NCT01084005" & race == "asian"),
    !(nct_id == "NCT01792518" & race == "black_aa")
  ) %>% 
  select(-c(race_to_exclude, exclude))

# Nest by trial
# Removed NCT00800683 and NCT01084005 due to durn convergence issues
# Removed NCT00819091 because all estimates were NA
# 74 trials remain
b_nest <- b_df_race %>% 
  filter(!(nct_id %in% c("NCT00800683", "NCT01084005", "NCT00819091"))) %>% 
  select(nct_id, id, weeks, event, trttype, age10:n_atc, -category) %>%
  group_by(nct_id) %>% 
  nest() %>% 
  arrange(nct_id)

# Addressing convergence issues ------------------------------------------------

# Trial lookups for omitting comorbs and conmeds
z_rm_comorbs <- c("NCT01064687", "NCT01191268", "NCT01289119", "NCT01648582",
                  "NCT01084005", "NCT01214239")
z_rm_conmeds <- c("NCT00996658", "NCT01194830", "NCT01438814")

# Address convergence issues
c_cnvg <- b_nest %>% 
  mutate(
    
    data_prep = map2(data, nct_id, function(df, id) {
      
      # NA smoke in NCT00996658, only two current smokers
      if (id %in% c("NCT00996658", "NCT00654381", "NCT01084005")) {df$smoke <- NA}
      
      # NA comorbs in set, very few non-zeros, any != 0 completed trial
      if (id %in% z_rm_comorbs) {df$comorbs <- NA}
      
      # NA con meds in set, very few non-zeros, any != 0 vast majority completed
      if (id %in% z_rm_conmeds) {df$n_atc <- NA}
      
      # Remove race for NCT01708902 (two levels, very few in second)
      if (id == "NCT01708902") {df$race <- NA}
      
      # Remove NA variables
      # 8 trials - 1 durn, 1 smoke, 4 comorbs, 3 con meds
      df <- df %>% select(id, weeks, event, where(~ all(!is.na(.))))
      df
      
    }),
    
    data_race = map(data_prep, function(df) {
      
      # Set race to NA if there is one unique (non-NA) level
      if ("race" %in% names(df) && length(unique(df$race)) == 1) {df$race <- NA}
      
      # Factorise race
      if ("race" %in% names(df)) {
        df <- df %>% 
          mutate(
            race = factor(
              race, 
              levels = c("white", setdiff(unique(df$race), "white"))
            )
          )
      }
      
      # Remove NA race
      df <- df %>% select(id, weeks, event, where(~ all(!is.na(.))))
      df
      
    })
  )

# Cox models  ------------------------------------------------------------------

# Fit models
cox <- c_cnvg %>% 
  mutate(
    cox = map(data_race, function(df) {
      vars <- names(df)
      event_index <- which(vars == "event")
      covars <- vars[(event_index + 1):length(vars)]
      form <- as.formula(paste("Surv(weeks, event) ~", paste(covars, collapse = " + ")))
      mod1 <- coxph(form, data = df, control = coxph.control(iter.max = 100))
      mod1
    })
  )

# Get outputs (results, diagnostics)
cox_output <- cox %>% 
  select(nct_id, cox) %>% 
  mutate(
    res = map(cox, ~ if (!is.null(.x)) {tidy(.x)} else {NA}),
    diag = map(cox, ~ if (!is.null(.x)) {glance(.x)} else {NA})
  )

# Extract variance-covariance matrices
cox_vcov <- cox %>% 
  select(nct_id, cox) %>% 
  mutate(
    vcov = map(cox, function(df) {
      vcov(df)
    })
  )

# Extract proportional hazard checks
cox_ph <- cox %>% 
  select(nct_id, cox) %>% 
  mutate(
    ph = map(cox, ~ if (!is.null(.x)) {cox.zph(.x, term = FALSE)} else {NA}),
    ph_df = map(ph, ~ if (inherits(.x, "cox.zph")) {as_tibble(.x$table, rownames = "term")} else {NA})
  )

cox_ph <- cox_ph %>% select(nct_id, ph_df) %>% unnest(c(ph_df))

# Save
saveRDS(cox %>% select(nct_id, data_race, cox), "processed_data/cox_outputs/cox_mdls.rds")
saveRDS(cox_output, "processed_data/cox_outputs/cox_outputs.rds")
saveRDS(cox_vcov, "processed_data/cox_outputs/cox_vcov.rds")
saveRDS(cox_ph, "processed_data/cox_outputs/cox_ph.rds")
