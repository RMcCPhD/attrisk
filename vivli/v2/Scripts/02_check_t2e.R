
source("Scripts/00_config.R")
library(tidyverse)
library(survival)
library(broom)

# Import dataset
a_imp_df <- readRDS("Processed_data/01_attr_covariates_n91.rds")

# Inspect time-to-event data ---------------------------------------------------

# Get t2e dataset
b_t2e <- survfit(Surv(weeks, event) ~ nct_id, data = a_imp_df)

b_tidy <- tidy(b_t2e) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

b_plot <- b_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free") +
  theme_classic()

ggsave(
  "Scratch_data/t2e_checks/01_initial.png", 
  b_plot, 
  width = 15, 
  height = 10, 
  units = "in"
)

# b_tidy %>%
#   group_by(nct_id) %>%
#   slice_max(time, n = 1) %>%
#   ungroup() %>%
#   distinct(nct_id, time) %>%
#   write_csv("scratch_data/check_censoring.csv")

# Add protocol-defined end of follow-up for completers
b_fu_fix <- a_imp_df %>% 
  left_join(
    read_csv("Scratch_data/t2e_checks/check_censoring.csv") %>% 
      filter(ok == 0) %>% 
      select(nct_id, time_c)
  ) %>% 
  mutate(
    event = case_when(
      !is.na(time_c)
      & weeks > time_c
      & event == 1 ~ 0,
      TRUE ~ event
    ),
    weeks = case_when(
      !is.na(time_c) & weeks > time_c ~ time_c, 
      TRUE ~ weeks
    )
  ) %>% 
  select(-time_c)

# Reestimate t2e
c_t2e_new <- survfit(Surv(weeks, event) ~ nct_id, data = b_fu_fix)

c_tidy <- tidy(c_t2e_new) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

c_plot <- c_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free") +
  theme_classic()

ggsave(
  "Scratch_data/t2e_checks/02_protocol_fu.png", 
  c_plot, 
  width = 15, 
  height = 10, 
  units = "in"
)

# Final checks
b_fu_fix %>% count(nct_id)
b_fu_fix %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save with fixes to end of follow-up
saveRDS(b_fu_fix, "Processed_data/02_t2e_checked.rds")
