
source("Scripts/00_config.R")
library(tidyverse)
library(survival)
library(broom)

# Import dataset
a_imp_df <- readRDS("Processed_data/n95.rds")

# Get time-to-event data
b_t2e <- survfit(Surv(weeks, event) ~ nct_id, data = a_imp_df)

b_tidy <- tidy(b_t2e) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

b_plot <- b_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free")

ggsave(
  "Scratch_data/plot1.png", 
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
#   write_csv("scratch_data/extra_censoring.csv")

# Adjust censoring where there are extremes (wrong end of follow-up)
c_imp_fix <- read_csv("Scratch_data/extra_censoring.csv")

# 51 trials
c_fix <- a_imp_df %>% 
  left_join(c_imp_fix %>% select(nct_id, time_c)) %>% 
  group_by(nct_id) %>% 
  fill(time_c, .direction = "updown") %>% 
  mutate(event = case_when(weeks > time_c ~ 0, TRUE ~ event)) %>% 
  ungroup()

c_fix_n51 <- c_fix %>% filter(!is.na(time_c))

# Replot cumulative incidence
c_t2e <- survfit(Surv(weeks, event) ~ nct_id, data = c_fix_n51)

c_tidy <- tidy(c_t2e) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

c_plot <- c_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free")

ggsave(
  "Scratch_data/plot2_fix.png", 
  c_plot, 
  width = 15, 
  height = 10, 
  units = "in"
)

d_t2e <- survfit(Surv(weeks, event) ~ nct_id, data = c_fix)

d_tidy <- tidy(d_t2e) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

d_plot <- d_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free")

ggsave(
  "Scratch_data/plot3.png", 
  d_plot, 
  width = 15, 
  height = 10, 
  units = "in"
)

# Correct follow-up period
e_imp_fix <- read_csv("Scratch_data/extra_censoring.csv")

# 16 trials
e_fix <- c_fix %>% 
  left_join(e_imp_fix %>% select(nct_id, time_f)) %>% 
  group_by(nct_id) %>% 
  fill(time_f, .direction = "updown") %>% 
  mutate(weeks = case_when(weeks > time_f ~ time_f, TRUE ~ weeks)) %>% 
  ungroup()

e_t2e <- survfit(Surv(weeks, event) ~ nct_id, data = e_fix)

e_tidy <- tidy(e_t2e) %>% 
  mutate(nct_id = gsub("nct_id=", "", strata)) %>% 
  select(nct_id, time, estimate) %>% 
  mutate(cind = 1 - estimate) %>% 
  arrange(nct_id, time)

e_plot <- e_tidy %>% 
  ggplot(aes(x = time, y = cind)) +
  geom_step(colour = "black", linewidth = 0.5) +
  facet_wrap(~nct_id, scales = "free")

ggsave(
  "Scratch_data/plot4.png", 
  e_plot, 
  width = 15, 
  height = 10, 
  units = "in"
)

# Final checks
e_fix %>% count(nct_id)
e_fix %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save with fixes to end of follow-up
z_df <- e_fix %>% select(-time_c, -time_f)
saveRDS(z_df, "Processed_data/n95_t2e_fixes.rds")
