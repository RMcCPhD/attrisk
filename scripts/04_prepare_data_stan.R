
# Preparing STAN data
# Creating dummy variables to categorise numeric variables like sample size,
# number of sites and number of facilities

source("scripts/00_functions_packages_etc.R")

# Import dataset
a_imp_df <- readRDS("processed_data/df_fu_n348.rds")
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
    
    # Sample size (<500 - 500-<750 - >=750)
    n_500 = if_else(n < 500, 1, 0),
    n_750 = if_else(n >= 500 & n < 750, 1, 0),
    n_abv = if_else(n >= 750, 1, 0),
    
    # Number of arms (2 - 3 - >3)
    arm_2 = if_else(n_arms == 2, 1, 0),
    arm_3 = if_else(n_arms == 3, 1, 0),
    arm_abv = if_else(n_arms > 3, 1, 0),
    
    # Number of sites (<100 vs >= 100)
    site_100 = if_else(n_facilities < 100, 1, 0),
    
    # Number of regions (1 or >=1)
    reg_1 = if_else(n_regions == 1, 1, 0),

    # Trial duration (<6m - >=6m)
    dur_6m = if_else(duration < 26, 1, 0),
    
    # Number of primary outcomes (1 - >= 1)
    pom_1 = if_else(n_poms == 1, 1, 0),
    
  )

# Shorten names
# Remove phase 4 or other (either phase 3 or not for model stability)
# Remove region counts and just use WHO regions instead (same thing)
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
    trial_id, arm_n, arm_n_attr, pt, trttype, p3, p4, po,
    starts_with("b_"),
    starts_with("reg_"), 
    everything()
  ) %>% 
  select(-c(source:n_regions)) %>% 
  select(-c(
    p4, po, reg_1
  )) %>% 
  rename(b_dbl = b_dlb)

glimpse(a_shrt)

for (i in seq_along(a_shrt)) {
  nm <- names(a_shrt[i])
  print(nm)
  sum <- table(a_shrt[i])
  print(sum)
}

# Check for collinearity between dummy variables
# Model falls over and suspect this
X_dummies <- a_shrt %>% select(trttype:pom_1)
X_cor <- cor(as.matrix(X_dummies))
X_high_cor <- which(abs(X_cor) > 0.8 & lower.tri(X_cor), arr.ind = TRUE)
colnames(X_dummies)[unique(c(X_high_cor[,1], X_high_cor[,2]))]

covs <- names(a_shrt[-c(1:5)])

a_shrt[, covs] <- lapply(a_shrt[, covs], as.integer)

b_stan_dat <- list(
  I = nrow(a_shrt),
  K = length(covs),
  r = a_shrt$arm_n_attr,
  N = a_shrt$arm_n,
  dt = a_shrt$pt,
  arm = as.integer(a_shrt$trttype),
  X = as.matrix(a_shrt[, covs])
)

glimpse(b_stan_dat)

# Save
saveRDS(b_stan_dat, "processed_data/stan_dat.rds")
