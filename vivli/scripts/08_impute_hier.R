
# Fit hierarchical models
# Predict values for completely missing variables within trials
# Save covariates with full imputation implemented
source("Scripts/00_config.R")
source("Scripts/00_functions.R")

set.seed(123)

library(tidyverse)
library(broom)
library(lme4)

# Import data with imputation for variables that were not completely missing
a_imp_df <- readRDS("Processed_data/imputation_part1.rds")

# Plot distributions
a_imp_df %>% 
  pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "skyblue") +
  facet_wrap(~var, scales = "free")

a_imp_df %>% reframe(across(everything(), ~ sum(is.na(.))))

# Summarise missingness
# 34 trials, different variables
# Save as log
a_log_imp <- a_imp_df %>% 
  group_by(nct_id) %>% 
  reframe(n = n(), across(bmi:n_atc, ~ sum(is.na(.)))) %>% 
  filter(if_any(bmi:n_atc, ~ . == n))

write_csv(a_log_imp, "Scratch_data/attrisk_imp_hier_log.csv")

# Workflow
# Fit hierarchical model for each covariate in above
# Remove trial(s) missing the covariate from specification
# Fit model and get outputs
# Predict covariates for that trial

# Note
# Random intercepts model - random slopes only for age, sex and bmi
# Takes quite long to run (haven't seen one fit)

# BMI --------------------------------------------------------------------------

bmi_df <- a_imp_df %>% 
  anti_join(
    a_log_imp %>% 
      filter(bmi != 0) %>% 
      select(nct_id)
  )

bmi_cmpl <- bmi_df %>% mutate(across(dbp:n_atc, ~ if_else(is.na(.), 0, .)))

bmi_p <- paste0(names(bmi_cmpl[c(4:9, 11:18)]), collapse = " + ")
bmi_lp <- paste0("bmi ~ ", bmi_p, " + (1 + age10 + sex | nct_id)")
bmi_fml <- as.formula(bmi_lp)

bmi_mdl <- lmer(formula = bmi_fml, data = bmi_cmpl)
summary(bmi_mdl)

# New addition
# Replace other NA columns with 0, otherwise can't predict
bmi_new <- a_imp_df %>% 
  inner_join(
    a_log_imp %>% 
      filter(bmi != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, bmi, race), ~ if_else(is.na(.), 0, .)))

bmi_new$bmi_pred <- predict(
  bmi_mdl, 
  newdata = bmi_new, 
  allow.new.levels = TRUE
)

bmi_new <- bmi_new %>% mutate(bmi = bmi_pred) %>% select(-bmi_pred)
added_bmi <- bmi_df %>% full_join(bmi_new)

# Check result
a_imp_df %>% reframe(na = sum(is.na(bmi)))
added_bmi %>% reframe(na = sum(is.na(bmi)))

a_imp_df %>% 
  ggplot(aes(x = bmi)) + 
  geom_density() + 
  geom_density(data = added_bmi, aes(x = bmi), colour = "green") +
  theme_bw()

# Diastolic bp -----------------------------------------------------------------

dbp_df <- added_bmi %>%
  anti_join(
    a_log_imp %>% 
      filter(dbp != 0) %>% 
      select(nct_id)
  )

dbp_cmpl <- dbp_df %>% mutate(across(sbp:n_atc, ~ if_else(is.na(.), 0, .)))

dbp_p <- paste0(names(dbp_cmpl[c(4:10, 12:18)]), collapse = " + ")
dbp_lp <- paste0("dbp ~ ", dbp_p, " + (1 + age10 + sex | nct_id)")
dbp_fml <- as.formula(dbp_lp)

dbp_mdl <- lmer(formula = dbp_fml, data = dbp_cmpl)
summary(dbp_mdl)

dbp_new <- added_bmi %>% 
  inner_join(
    a_log_imp %>% 
      filter(dbp != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, dbp, race), ~ if_else(is.na(.), 0, .)))

dbp_new$dbp_pred <- predict(
  dbp_mdl, 
  newdata = dbp_new, 
  allow.new.levels = TRUE
)

dbp_new <- dbp_new %>% mutate(dbp = dbp_pred) %>% select(-dbp_pred)
added_dbp <- dbp_df %>% full_join(dbp_new)

a_imp_df %>% reframe(na = sum(is.na(dbp)))
added_dbp %>% reframe(na = sum(is.na(dbp)))

a_imp_df %>% 
  ggplot(aes(x = dbp)) + 
  geom_density() + 
  geom_density(data = added_dbp, aes(x = dbp), colour = "green") +
  theme_bw()

# Systolic bp ------------------------------------------------------------------

sbp_df <- added_dbp %>%
  anti_join(
    a_log_imp %>% 
      filter(sbp != 0) %>% 
      select(nct_id)
  )

sbp_cmpl <- sbp_df %>% mutate(across(egfr:n_atc, ~ if_else(is.na(.), 0, .)))

sbp_p <- paste0(names(sbp_cmpl[c(4:11, 13:18)]), collapse = " + ")
sbp_lp <- paste0("sbp ~ ", sbp_p, " + (1 + age10 + sex | nct_id)")
sbp_fml <- as.formula(sbp_lp)

sbp_mdl <- lmer(formula = sbp_fml, data = sbp_cmpl)
summary(sbp_mdl)

sbp_new <- added_dbp %>% 
  inner_join(
    a_log_imp %>% 
      filter(sbp != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, sbp, race), ~ if_else(is.na(.), 0, .)))

sbp_new$sbp_pred <- predict(
  sbp_mdl, 
  newdata = sbp_new, 
  allow.new.levels = TRUE
)

sbp_new <- sbp_new %>% mutate(sbp = sbp_pred) %>% select(-sbp_pred)
added_sbp <- sbp_df %>% full_join(sbp_new)

a_imp_df %>% reframe(na = sum(is.na(sbp)))
added_sbp %>% reframe(na = sum(is.na(sbp)))

a_imp_df %>% 
  ggplot(aes(x = sbp)) + 
  geom_density() + 
  geom_density(data = added_sbp, aes(x = sbp), colour = "green") +
  theme_bw()

# eGFR -------------------------------------------------------------------------

egfr_df <- added_sbp %>%
  anti_join(
    a_log_imp %>% 
      filter(egfr != 0) %>% 
      select(nct_id)
  )

egfr_cmpl <- egfr_df %>% mutate(across(chol_hdl:n_atc, ~ if_else(is.na(.), 0, .)))

egfr_p <- paste0(names(egfr_cmpl[c(4:12, 14:18)]), collapse = " + ")
egfr_lp <- paste0("egfr ~ ", egfr_p, " + (1 + age10 + sex | nct_id)")
egfr_fml <- as.formula(egfr_lp)

egfr_mdl <- lmer(formula = egfr_fml, data = egfr_cmpl)
summary(egfr_mdl)

# New addition
# Replace other NA columns with 0, otherwise can't predict
egfr_new <- added_sbp %>% 
  inner_join(
    a_log_imp %>% 
      filter(egfr != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, egfr, race), ~ if_else(is.na(.), 0, .)))

egfr_new$egfr_pred <- predict(
  egfr_mdl, 
  newdata = egfr_new, 
  allow.new.levels = TRUE
)

egfr_new <- egfr_new %>% mutate(egfr = egfr_pred) %>% select(-egfr_pred)
added_egfr <- egfr_df %>% full_join(egfr_new)

a_imp_df %>% reframe(na = sum(is.na(egfr)))
added_egfr %>% reframe(na = sum(is.na(egfr)))

a_imp_df %>% 
  ggplot(aes(x = egfr)) + 
  geom_density() + 
  geom_density(data = added_egfr, aes(x = egfr), colour = "green") +
  theme_bw()

# HDL cholesterol --------------------------------------------------------------

chol_hdl_df <- added_egfr %>%
  anti_join(
    a_log_imp %>% 
      filter(chol_hdl != 0) %>% 
      select(nct_id)
  )

chol_hdl_cmpl <- chol_hdl_df %>% mutate(across(chol_ldl:n_atc, ~ if_else(is.na(.), 0, .)))

chol_hdl_p <- paste0(names(chol_hdl_cmpl[c(4:13, 15:18)]), collapse = " + ")
chol_hdl_lp <- paste0("chol_hdl ~ ", chol_hdl_p, " + (1 + age10 + sex | nct_id)")
chol_hdl_fml <- as.formula(chol_hdl_lp)

chol_hdl_mdl <- lmer(formula = chol_hdl_fml, data = chol_hdl_cmpl)
summary(chol_hdl_mdl)

# New addition
# Replace other NA columns with 0, otherwise can't predict
chol_hdl_new <- added_egfr %>% 
  inner_join(
    a_log_imp %>% 
      filter(chol_hdl != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, chol_hdl, race), ~ if_else(is.na(.), 0, .)))

chol_hdl_new$chol_hdl_pred <- predict(
  chol_hdl_mdl, 
  newdata = chol_hdl_new, 
  allow.new.levels = TRUE
)

chol_hdl_new <- chol_hdl_new %>% mutate(chol_hdl = chol_hdl_pred) %>% select(-chol_hdl_pred)
added_chol_hdl <- chol_hdl_df %>% full_join(chol_hdl_new)

a_imp_df %>% reframe(na = sum(is.na(chol_hdl)))
added_chol_hdl %>% reframe(na = sum(is.na(chol_hdl)))

a_imp_df %>% 
  ggplot(aes(x = chol_hdl)) + 
  geom_density() + 
  geom_density(data = added_chol_hdl, aes(x = chol_hdl), colour = "green") +
  theme_bw()

# LDL cholesterol --------------------------------------------------------------

chol_ldl_df <- added_chol_hdl %>%
  anti_join(
    a_log_imp %>% 
      filter(chol_ldl != 0) %>% 
      select(nct_id)
  )

chol_ldl_cmpl <- chol_ldl_df %>% mutate(across(chol_total:n_atc, ~ if_else(is.na(.), 0, .)))

chol_ldl_p <- paste0(names(chol_ldl_cmpl[c(4:14, 16:18)]), collapse = " + ")
chol_ldl_lp <- paste0("chol_ldl ~ age10 + sex + (1 + age10 + sex | nct_id)")
chol_ldl_fml <- as.formula(chol_ldl_lp)

chol_ldl_mdl <- lmer(formula = chol_ldl_fml, data = chol_ldl_cmpl)
summary(chol_ldl_mdl)

# New addition
# Replace other NA columns with 0, otherwise can't predict
chol_ldl_new <- added_chol_hdl %>% 
  inner_join(
    a_log_imp %>% 
      filter(chol_ldl != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, chol_ldl, race), ~ if_else(is.na(.), 0, .)))

chol_ldl_new$chol_ldl_pred <- predict(
  chol_ldl_mdl, 
  newdata = chol_ldl_new, 
  allow.new.levels = TRUE
)

chol_ldl_new <- chol_ldl_new %>% mutate(chol_ldl = chol_ldl_pred) %>% select(-chol_ldl_pred)
added_chol_ldl <- chol_ldl_df %>% full_join(chol_ldl_new)

a_imp_df %>% reframe(na = sum(is.na(chol_ldl)))
added_chol_ldl %>% reframe(na = sum(is.na(chol_ldl)))

a_imp_df %>% 
  ggplot(aes(x = chol_ldl)) + 
  geom_density() + 
  geom_density(data = added_chol_ldl, aes(x = chol_ldl), colour = "green") +
  theme_bw()

# Total cholesterol ------------------------------------------------------------

chol_total_df <- added_chol_ldl %>%
  anti_join(
    a_log_imp %>% 
      filter(chol_total != 0) %>% 
      select(nct_id)
  )

chol_total_cmpl <- chol_total_df %>% mutate(across(comorbs:n_atc, ~ if_else(is.na(.), 0, .)))

chol_total_p <- paste0(names(chol_total_cmpl[c(4:15, 17:18)]), collapse = " + ")
chol_total_lp <- paste0("chol_total ~ age10 + sex + chol_hdl + (1 + age10 + sex + chol_hdl | nct_id)")
chol_total_fml <- as.formula(chol_total_lp)

chol_total_mdl <- lmer(formula = chol_total_fml, data = chol_total_cmpl)
summary(chol_total_mdl)

# New addition
# Replace other NA columns with 0, otherwise can't predict
chol_total_new <- added_chol_ldl %>% 
  inner_join(
    a_log_imp %>% 
      filter(chol_total != 0) %>% 
      select(nct_id)
  ) %>% 
  mutate(across(-c(nct_id:trttype, chol_total, race), ~ if_else(is.na(.), 0, .)))

chol_total_new$chol_total_pred <- predict(
  chol_total_mdl, 
  newdata = chol_total_new, 
  allow.new.levels = TRUE
)

chol_total_new <- chol_total_new %>% mutate(chol_total = chol_total_pred) %>% select(-chol_total_pred)
added_chol_total <- chol_total_df %>% full_join(chol_total_new)

a_imp_df %>% reframe(na = sum(is.na(chol_total)))
added_chol_total %>% reframe(na = sum(is.na(chol_total)))

a_imp_df %>% 
  ggplot(aes(x = chol_total)) + 
  geom_density() + 
  geom_density(data = added_chol_total, aes(x = chol_total), colour = "green") +
  theme_bw()

# Comorbidities ----------------------------------------------------------------

comorbs_df <- added_chol_total %>%
  anti_join(
    a_log_imp %>% 
      filter(comorbs != 0) %>% 
      select(nct_id)
  )

comorbs_cmpl <- comorbs_df %>% mutate(n_atc = if_else(is.na(n_atc), 0, n_atc))

comorbs_p <- paste0(names(comorbs_cmpl[c(4:16, 18)]), collapse = " + ")
comorbs_lp <- paste0("comorbs ~ ", comorbs_p, " + (1 + age10 + sex | nct_id)")
comorbs_fml <- as.formula(comorbs_lp)

comorbs_mdl <- lmer(formula = comorbs_fml, data = comorbs_cmpl)
summary(comorbs_mdl)

comorbs_new <- added_chol_total %>% 
  inner_join(
    a_log_imp %>% 
      filter(comorbs != 0) %>% 
      select(nct_id)
  ) %>%
  mutate(across(-c(nct_id:trttype, comorbs, race), ~ if_else(is.na(.), 0, .)))

comorbs_new$comorbs_pred <- predict(
  comorbs_mdl, 
  newdata = comorbs_new, 
  allow.new.levels = TRUE
)

comorbs_new <- comorbs_new %>% mutate(comorbs = comorbs_pred) %>% select(-comorbs_pred)
added_comorbs <- comorbs_df %>% full_join(comorbs_new)

a_imp_df %>% reframe(na = sum(is.na(comorbs)))
added_comorbs %>% reframe(na = sum(is.na(comorbs)))

a_imp_df %>% 
  ggplot(aes(x = comorbs)) + 
  geom_density() + 
  geom_density(data = added_comorbs, aes(x = comorbs), colour = "green") +
  theme_bw()

# Con meds ---------------------------------------------------------------------

n_atc_df <- added_comorbs %>%
  anti_join(
    a_log_imp %>% 
      filter(n_atc != 0) %>% 
      select(nct_id)
  )

n_atc_cmpl <- n_atc_df %>% mutate(n_atc = if_else(is.na(n_atc), 0, n_atc))

n_atc_p <- paste0(names(n_atc_cmpl[c(4:17)]), collapse = " + ")
n_atc_lp <- paste0("n_atc ~ ", n_atc_p, " + (1 + age10 + sex | nct_id)")
n_atc_fml <- as.formula(n_atc_lp)

n_atc_mdl <- lmer(formula = n_atc_fml, data = n_atc_cmpl)
summary(n_atc_mdl)

n_atc_new <- added_comorbs %>% 
  inner_join(
    a_log_imp %>% 
      filter(n_atc != 0) %>% 
      select(nct_id)
  ) %>%
  mutate(across(-c(nct_id:trttype, n_atc, race), ~ if_else(is.na(.), 0, .)))

n_atc_new$n_atc_pred <- predict(
  n_atc_mdl, 
  newdata = n_atc_new, 
  allow.new.levels = TRUE
)

n_atc_new <- n_atc_new %>% mutate(n_atc = n_atc_pred) %>% select(-n_atc_pred)
added_n_atc <- n_atc_df %>% full_join(n_atc_new)

a_imp_df %>% reframe(na = sum(is.na(n_atc)))
added_n_atc %>% reframe(na = sum(is.na(n_atc)))

a_imp_df %>% 
  ggplot(aes(x = n_atc)) + 
  geom_density() + 
  geom_density(data = added_n_atc, aes(x = n_atc), colour = "green") +
  theme_bw()

# Overall checks ---------------------------------------------------------------

a_imp_df %>% reframe(across(age10:n_atc,  ~ sum(is.na(.))))
added_n_atc %>% reframe(across(age10:n_atc,  ~ sum(is.na(.))))

a_imp_df %>% 
  reframe(
    across(
      c(age10, hba1c:n_atc),
      ~ unlist(
        list(
          mean = mean(., na.rm = TRUE), 
          min = min(., na.rm = TRUE), 
          max = max(., na.rm = TRUE)
        )
      )
    )
  )

added_n_atc %>% 
  reframe(
    across(
      c(age10, hba1c:n_atc),
      ~ unlist(
        list(
          mean = mean(., na.rm = TRUE), 
          min = min(., na.rm = TRUE), 
          max = max(., na.rm = TRUE)
        )
      )
    )
  )

# Put comorbidities predicted close to but less than 0 as 0
# Format
d_tidy <- added_n_atc %>% 
  mutate(
    across(c(comorbs:n_atc), ~ if_else(. < 0, 0, .)),
    across(c(sex:smoke, comorbs:n_atc), as.integer)
  )

d_tidy %>% reframe(across(age10:n_atc,  ~ sum(is.na(.))))

plot_dists <- d_tidy %>% 
  pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density() +
  geom_density(
    data = a_imp_df %>% 
      pivot_longer(c(age10, hba1c:n_atc), names_to = "var", values_to = "val"),
    aes(x = val),
    colour = "forestgreen", alpha = 0.2
  ) +
  facet_wrap(~ var, scales = "free")

# Save outputs -----------------------------------------------------------------

# Plots of distribution
ggsave(
  "Outputs/imputation/plot_before_after_hier_imp.png",
  plot_dists,
  width = 8,
  height = 6,
  units = "in"
)

# Hierarchical-imputation dataset
saveRDS(d_tidy, "Processed_data/imputed_variables.rds")
