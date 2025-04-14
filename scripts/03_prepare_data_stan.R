
# Preparing STAN data
# Creating dummy variables to categorise numeric variables like sample size,
# number of sites and number of facilities

source("scripts/00_functions_packages_etc.R")

# Import dataset
a_imp_df <- read_csv("processed_data/df_stan_n348.csv")
glimpse(a_imp_df)

# Inspecting variables to determine cut-offs
summary(a_imp_df$n)
summary(a_imp_df$n_arms)
summary(a_imp_df$n_facilities)
summary(a_imp_df$n_regions)
summary(a_imp_df$duration)
summary(a_imp_df$n_poms)

# Create dummy variables
a_dmy <- a_imp_df %>% 
  mutate(
    
    # Sample size (<250 - 250-<500 - 500-<750 - 750-<1000 - >=1000)
    n_250 = if_else(n < 250, 1, 0),
    n_500 = if_else(n >= 250 & n < 500, 1, 0),
    n_750 = if_else(n >= 500 & n < 750, 1, 0),
    n_1000 = if_else(n >= 750 & n < 1000, 1, 0),
    n_abv = if_else(n >= 1000, 1, 0),
    
    # Number of arms (2 - 3 - <3)
    arm_2 = if_else(n_arms == 2, 1, 0),
    arm_3 = if_else(n_arms == 3, 1, 0),
    arm_abv = if_else(n_arms > 3, 1, 0),
    
    # Number of sites (<50 - 50-<100 - 100-<150 - >=150)
    site_50 = if_else(n_facilities < 50, 1, 0),
    site_100 = if_else(n_facilities >= 50 & n_facilities < 100, 1, 0),
    site_150 = if_else(n_facilities >= 100 & n_facilities < 150, 1, 0),
    site_abv = if_else(n_facilities >= 150, 1, 0),
    
    # Number of regions (1 - 2 - 3 - 4 - 5 - 6)
    reg_1 = if_else(n_regions == 1, 1, 0),
    reg_2 = if_else(n_regions == 2, 1, 0),
    reg_3 = if_else(n_regions == 3, 1, 0),
    reg_4 = if_else(n_regions == 4, 1, 0),
    reg_5 = if_else(n_regions == 5, 1, 0),
    reg_6 = if_else(n_regions == 6, 1, 0),
    
    # Trial duration (<6m - 6m-<1yr - >=1yr)
    dur_6m = if_else(duration < 26, 1, 0),
    dur_1yr = if_else(duration >= 26 & duration < 52, 1, 0),
    dur_abv = if_else(duration >= 52, 1, 0),
    
    # Number of primary outcomes (1 - 2 - 3 - 4 - 5 - 6)
    pom_1 = if_else(n_poms == 1, 1, 0),
    pom_2 = if_else(n_poms == 2, 1, 0),
    pom_3 = if_else(n_poms == 3, 1, 0),
    pom_4 = if_else(n_poms == 4, 1, 0),
    pom_5 = if_else(n_poms == 5, 1, 0),
    pom_6 = if_else(n_poms == 6, 1, 0)
    
  )

# Shorten names
a_shrt <- a_dmy %>% 
  rename(
    p3 = phase_3,
    p4 = phase_4,
    po = phase_notapp,
    b_dlb = double,
    b_ol = open_label,
    b_q = quadruple,
    b_t = triple,
    reg_am = region_of_the_americas,
    reg_eur = european_region,
    reg_pac = western_pacific_region,
    reg_afr = african_region,
    reg_sea = south_east_asia_region,
    reg_med = eastern_mediterranean_region
  ) %>% 
  select(
    trial_id, arm_n, arm_n_attr, duration, trttype, p3, p4, po,
    starts_with("b_"),
    starts_with("reg_"), 
    everything()
  ) %>% 
  select(-c(source:n_regions))

glimpse(a_shrt)

# Save
saveRDS(a_shrt, "processed_data/df_stan_ready.rds")
