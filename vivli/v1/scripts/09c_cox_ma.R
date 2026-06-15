
# Cox meta-analysis
# Final check before export
# Cox checks - proportional hazards, non-linearity
source("Scripts/00_config.R")
source("Scripts/00_functions.R")

library(tidyverse)
library(metafor)
library(broom)

# Import results
a_imp_res <- readRDS("processed_data/cox_outputs/cox_outputs.rds") %>% 
  select(nct_id, res) %>% 
  unnest(c(res))

# Import vcov
a_imp_vcov <- readRDS("processed_data/cox_outputs/cox_vcov.rds") %>% 
  select(nct_id, vcov)

# Fit trial-level plots
plot_trial_level <- a_imp_res %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
  ) %>% 
  ggplot(aes(x = nct_id, y = estimate, ymin = lci, ymax = uci)) +
  geom_point() +
  geom_linerange() +
  scale_y_continuous(n.breaks = 15) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~term, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "Pooled log-hazard")

plot_trial_level

ggsave(
  "Outputs/cox_trial_level/res/trial_plots_linear.png",
  plot_trial_level,
  height = 20,
  width = 20,
  units = "in"
)

# Meta-analyse (independent covariates) ----------------------------------------

cox_ma_lst <- split(a_imp_res, a_imp_res$term)

cox_ma <- lapply(cox_ma_lst, function(df) {
  rma.mv(
    yi = df$estimate, 
    V = df$std.error ^ 2, 
    random = ~ 1 | nct_id,
    data = df
  )
})

cox_ma_tidy <- bind_rows(
  lapply(names(cox_ma), function(nm) {
    tidy(cox_ma[[nm]]) %>% 
      mutate(model_term = nm)
  })
)

# Add SDs for labelling
cox_ma_add_sd <- cox_ma_tidy %>% 
  left_join(
    read_csv("Scratch_data/attrisk_var_scaling.csv") %>% 
      rename(
        model_term = var,
        sd = scaling_sd
      )
  ) %>% 
  mutate(
    label = case_when(
      !is.na(sd) ~ paste0(model_term, " [", round(sd, 1), "]"),
      TRUE ~ model_term
    )
  ) %>% 
  select(-sd)

# Plot pooled log-hazard
plot_cox_ma <- cox_ma_add_sd %>% 
  arrange(desc(estimate)) %>% 
  ggplot(
    aes(
      x = label, 
      y = estimate, 
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    )
  ) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 15) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red") +
  theme_bw() +
  labs(x = NULL, y = "Log estimates")

plot_cox_ma

ggsave(
  "Outputs/cox_trial_level/log_hr_ma.png",
  plot_cox_ma,
  height = 6,
  width = 8,
  units = "in"
)

# Meta-analyse (no independence) -----------------------------------------------

cox_ma2_id <- a_imp_res$nct_id
cox_ma2_params <- a_imp_res$term
cox_ma2_ests <- a_imp_res$estimate
cox_ma2_vcov <- as.matrix(bdiag(a_imp_vcov$vcov))

cox_ma2 <- rma.mv(
  yi = cox_ma2_ests,
  V = cox_ma2_vcov,
  mods = ~ factor(cox_ma2_params) - 1,
  random = ~ cox_ma2_params | cox_ma2_id
)

cox_ma2_tidy <- tidy(cox_ma2)

# Plot pooled log-hazard
plot_cox_ma2 <- cox_ma2_tidy %>% 
  mutate(term = gsub("factor\\(cox_ma2_params\\)", "", term)) %>% 
  arrange(desc(estimate)) %>% 
  ggplot(
    aes(
      x = term, 
      y = estimate, 
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    )
  ) +
  geom_point() + 
  geom_linerange() +
  scale_y_continuous(n.breaks = 15) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red") +
  theme_bw() +
  labs(x = NULL, y = "Log estimates")

plot_cox_ma2

ggsave(
  "Outputs/cox_trial_level/log_hr_ma2.png",
  plot_cox_ma2,
  height = 6,
  width = 8,
  units = "in"
)
