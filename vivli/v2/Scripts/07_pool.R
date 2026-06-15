
# Pool participant-level effects
# Results not exported - serves to examine effects

source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(metafor)
library(broom)

# Import model outputs
a_imp_out <- readRDS("Processed_data/06_cox_out.rds")

# Trial-level plots
# Truncating at +/-10 logHR to make viewing easier
b_trial <- a_imp_out %>% 
  select(sponsor, nct_id, res) %>% 
  unnest(res) %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
    hr = exp(estimate),
    hr_lci = exp(lci),
    hr_uci = exp(uci),
    across(lci:hr_uci, ~ round(., 8))
  ) %>% 
  filter(!is.na(estimate)) %>% 
  mutate(
    term = gsub("term", "", term),
    term = factor(
      term,
      levels = c(
        "age10", "sex", "raceasian", "raceblack_aa", "raceother_na_unk",
        "smoke", "durn", "bmi", "comorbs", "conmeds", "hba1c", "egfr",
        "sbp", "chol_hdl", "chol_total", "trttype"
      )
    )
  )

b_trial_plot <- b_trial %>% 
  ggplot(aes(x = estimate, y = nct_id)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(xmin = lci, xmax = uci)) +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  facet_wrap(~term, scales = "free_y") +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_classic() +
  theme(axis.text.y = element_blank())

ggsave(
  "Outputs/trial_loghr.pdf",
  b_trial_plot,
  width = 10,
  height = 10,
  units = "in"
)

# Meta-analysis
b_ma_out <- a_imp_out %>% 
  select(sponsor, nct_id, res) %>%
  unnest(res) %>% 
  ungroup()

b_ma_vcov <- a_imp_out %>% 
  select(sponsor, nct_id, vcov) %>% 
  ungroup()

b_ma_id <- b_ma_out$nct_id
b_ma_params <- b_ma_out$term
b_ma_ests <- b_ma_out$estimate
b_ma_vcov_bdiag <- as.matrix(bdiag(b_ma_vcov$vcov))

b_ma <- rma.mv(
  yi = b_ma_ests,
  V = b_ma_vcov_bdiag,
  mods = ~ factor(b_ma_params) - 1,
  random = ~ b_ma_params | b_ma_id,
  # sparse = TRUE,
  verbose = TRUE
)

b_ma_tidy <- tidy(b_ma) %>% 
  mutate(
    term = gsub("factor\\(b_ma_params\\)", "", term),
    term = factor(
      term,
      levels = c(
        "age10", "sex", "raceasian", "raceblack_aa", "raceother_na_unk",
        "smoke", "durn", "bmi", "comorbs", "conmeds", "hba1c", "egfr",
        "sbp", "chol_hdl", "chol_total", "trttype"
      )
    ),
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
    hr = exp(estimate),
    hr_lci = exp(lci),
    hr_uci = exp(uci),
    across(lci:hr_uci, ~ round(., 8))
  )

# Plot pooled log-hazard
b_loghr_plot <- b_ma_tidy %>% 
  arrange(desc(estimate)) %>% 
  ggplot(
    aes(
      x = fct_rev(term), 
      y = estimate, 
      ymin = lci,
      ymax = uci
    )
  ) +
  geom_point() + 
  geom_linerange(size = 0.5) +
  scale_y_continuous(limits = c(-0.6, 0.6), n.breaks = 10) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  theme_classic() +
  labs(x = NULL, y = "Pooled log-hazard ratio")

b_loghr_plot

# Plot pooled HR
b_hr_plot <- b_ma_tidy %>% 
  arrange(desc(estimate)) %>% 
  ggplot(
    aes(
      x = fct_rev(term), 
      y = hr, 
      ymin = hr_lci,
      ymax = hr_uci
    )
  ) +
  geom_point() + 
  geom_linerange(size = 0.5) +
  scale_y_continuous(limits = c(0.5, 1.5), n.breaks = 10) +
  coord_flip() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  theme_classic() +
  labs(x = NULL, y = "Pooled hazard ratio")

b_hr_plot

# Save
ggsave(
  "Outputs/pooled_loghr.pdf",
  b_loghr_plot,
  height = 6,
  width = 6,
  units = "in"
)

ggsave(
  "Outputs/pooled_hr.pdf",
  b_hr_plot,
  height = 6,
  width = 6,
  units = "in"
)
