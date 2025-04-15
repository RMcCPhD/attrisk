
# Add mean follow-up time
# Available for 86 IPD trials, approximated for remaining
source("scripts/00_functions_packages_etc.R")

# 348 trials, 86 have fu ratios (IPD trials)
a_imp_df <- read_csv("processed_data/df_stan_n348.csv")
a_imp_fu <- read_csv("data/fu_ratio_n86.csv")

# Prepare fu ratio
b_fu <- a_imp_fu %>% 
  mutate(collapsed = strsplit(collapsed, " AND ")) %>% 
  unnest(collapsed) %>% 
  separate(
    collapsed,
    into = c("trttype", "mft_ratio"),
    sep = " "
  ) %>% 
  mutate(
    across(trttype:mft_ratio, ~ gsub(",", "", .)),
    across(trttype:mft_ratio, ~ gsub("\\b\\w+=", "", .)),
    trttype = if_else(trttype == "plc", 1, 0)
  ) %>% 
  rename(trial_id = nct_id)

# Join onto dataset for ipd trials
c_join_fu <- a_imp_df %>% 
  left_join(b_fu) %>% 
  select(trial_id:arm_n_attr, duration, mft_ratio, everything())

# Remove duplicates for IPD trials
d_rem_dups <- c_join_fu %>% 
  filter(source == "ipd") %>% 
  group_by(trial_id, trttype) %>% 
  slice_head(n = 1) %>% 
  ungroup()

d_no_dups <- c_join_fu %>% 
  filter(source != "ipd") %>% 
  full_join(d_rem_dups)

# Scale parameters for finding the closest IPD match to get ratio
e_scale <- d_no_dups %>% 
  mutate(
    scale_arm_n = scale(arm_n)[, 1],
    scale_arm_attr = scale(arm_n_attr)[, 1],
    scale_dur = scale(duration)[, 1]
  )

# Split data - placebo arms
e_agg_plc <- e_scale %>% filter(source == "agg", trttype == 1L)
e_ipd_plc <- e_scale %>% filter(source == "ipd", trttype == 1L)

# Find best match between ipd and agg trials per agg trial
# Take mft ratio for that ipd trial as the mft ratio for the agg trial
e_agg_plc_ratio <- e_agg_plc %>% 
  rowwise() %>% 
  mutate(app_mft_ratio = nn_ratio(cur_data(), e_ipd_plc)) %>% 
  ungroup()

# Repeat for active comparator arms
e_agg_ac <- e_scale %>% filter(source == "agg", trttype == 0L)
e_ipd_ac <- e_scale %>% filter(source == "ipd", trttype == 0L)

e_agg_ac_ratio <- e_agg_ac %>% 
  rowwise() %>% 
  mutate(app_mft_ratio = nn_ratio(cur_data(), e_ipd_ac)) %>% 
  ungroup()

# Combine aggregate with approximations
e_agg_bind <- bind_rows(e_agg_plc_ratio, e_agg_ac_ratio) %>% 
  select(-mft_ratio, -contains("scale_")) %>% 
  rename(mft_ratio = app_mft_ratio) %>% 
  select(trial_id:duration, mft_ratio, everything())

# Combine with ipd
e_bind <- bind_rows(e_agg_bind, d_no_dups %>% filter(source == "ipd"))
e_bind %>% count(trial_id, source) %>% count(source)

# Check for NAs
e_bind %>% reframe(across(everything(), ~ sum(is.na(.)))) %>% View()

# Calculate mean person-time
f_pt <- e_bind %>% 
  mutate(
    arm_n = if_else(trial_id == "NCT00313313", 267, arm_n),
    arm_n_attr = if_else(trial_id == "NCT00313313", 82, arm_n_attr),
    mft_ratio = as.numeric(mft_ratio),
    yr = duration/52,
    arm_n_compl = arm_n - arm_n_attr,
    pt = arm_n_compl * yr + arm_n_attr * (mft_ratio * yr),
    pt = pt/arm_n
  ) %>% 
  select(
    trial_id:arm_n, arm_n_compl, arm_n_attr, duration, yr, pt, everything()
  )

# Save
saveRDS(f_pt, "df_fu_n348.csv")
