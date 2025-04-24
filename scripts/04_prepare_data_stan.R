
# Preparing STAN data
# Creating dummy variables to categorise numeric variables like sample size,
# number of sites and number of facilities

source("scripts/00_functions_packages_etc.R")

# Import dataset
a_imp_df <- readRDS("processed_data/df_fu_n348.rds")
glimpse(a_imp_df)

# Removing n regions as it's redundant (same thing is achieved with
# multiple regions being 1)
b_rm <- a_imp_df %>% 
  select(-n_regions) %>% 
  mutate(
    across(c(trttype:triple, n_poms:eastern_mediterranean_region), as.integer)
  )

# Summarise variables
covs <- names(b_rm[c(3, 4, 7, 11:27)])
for (var in covs) {tidy_sum(b_rm, var)}

# Create dummy variables
c_dmy <- b_rm %>% 
  mutate(

    # Sample size (<250 - >=250-<500, >=500-<750 - >=750-<1000 - >=1000)
    n_250 = if_else(n < 250, 1, 0),
    n_500 = if_else(n >= 250 & n < 500, 1, 0),
    n_750 = if_else(n >= 500 & n < 750, 1, 0),
    n_1000 = if_else(n >= 750 & n < 1000, 1, 0),
    n_abv = if_else(n >= 1000, 1, 0),
    
    # Arm size (<100 - >=100-<200 - >=200-<300 - >=300)
    arm_100 = if_else(n < 100, 1, 0),
    arm_200 = if_else(n >= 100 & n < 200, 1, 0),
    arm_300 = if_else(n >= 200 & n < 300, 1, 0),
    arm_abv = if_else(n >= 300, 1, 0),
    
    # Trial duration (<6m - >=6m-<1yr - >=1yr)
    dur_6m = if_else(duration < 26, 1, 0),
    dur_1yr = if_else(duration >= 26 & duration < 52, 1, 0),
    dur_abv = if_else(duration >= 52, 1, 0),
    
    # Number of arms (2 - 3 - 4 - >4)
    arm_2 = if_else(n_arms == 2, 1, 0),
    arm_3 = if_else(n_arms == 3, 1, 0),
    arm_4 = if_else(n_arms == 4, 1, 0),
    arm_abv = if_else(n_arms > 4, 1, 0),
    
    # Phase (3 - 4)
    p3 = if_else(phase_3 == 1, 1, 0),
    p4 = if_else(phase_4 == 1, 1, 0),
    #po = if_else(phase_notapp == 1, 1, 0),
    
    # Blinding (double - triple - quadruple - open label)
    bdbl = double,
    btrip = triple,
    bquad = quadruple,
    bol = open_label,
    
    # N facilities (<100 - >=100-<200 - >=200)
    f_100 = if_else(n_facilities < 100, 1, 0),
    f_200 = if_else(n_facilities >= 100 & n_facilities < 200, 1, 0),
    f_abv = if_else(n_facilities >= 200, 1, 0),
    
    # N poms (1 - 2 - 3 - >3)
    pom_1 = if_else(n_poms == 1, 1, 0),
    pom_2 = if_else(n_poms == 2, 1, 0),
    pom_3 = if_else(n_poms == 3, 1, 0),
    pom_abv = if_else(n_poms >3, 1, 0),
    
    # Regions (am - eur - wpac - afr - sea - em)
    reg_am = region_of_the_americas,
    reg_eur = european_region,
    reg_wpac = western_pacific_region,
    reg_afr = african_region,
    reg_sea = south_east_asia_region,
    reg_em = eastern_mediterranean_region,
    
    duration = duration/52
    
  ) %>% 
  mutate(across(n_250:reg_em, as.integer))

# Save
saveRDS(c_dmy, "processed_data/dat_covars.rds")

# Check frequencies among dummy variables
covs <- names(c_dmy[c(28:61)])
for (var in covs) {tidy_sum(c_dmy, var)}

covs <- names(c_dmy[c(28:61)])

b_stan_dat <- list(
  I = nrow(c_dmy),
  K = length(covs),
  r = c_dmy$arm_n_attr,
  N = c_dmy$arm_n,
  dt = c_dmy$duration,
  arm = as.integer(c_dmy$trttype),
  X = as.matrix(c_dmy[, covs])
)

glimpse(b_stan_dat)

# Save
saveRDS(b_stan_dat, "processed_data/stan_dat.rds")
