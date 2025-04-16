
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
a_dmy <- b_rm %>% 
  mutate(
    
    # Sample size (<250 - >=250-<500, >=500-<750 - >=750-<1000 - >=1000)
    n_250 = if_else(n < 250, 1, 0),
    n_500 = if_else(n >= 250 & n < 500, 1, 0),
    n_750 = if_else(n >= 500 & n < 750, 1, 0),
    # n_1000 = if_else(n >= 750 & n < 1000, 1, 0),
    # n_abv = if_else(n >= 1000, 1, 0),
    
    # Arm size (<100 - >=100-<200 - >=200-<300 - >=300)
    arm_100 = if_else(arm_n < 100, 1, 0),
    arm_200 = if_else(arm_n >=100 & arm_n < 200, 1, 0),
    arm_300 = if_else(arm_n >= 200 & arm_n < 300, 1, 0),
    # arm_abv = if_else(arm_n >= 300, 1, 0),
    
    # Trial duration (<6m - >=6m-<1yr - >=1yr)
    dur_6m = if_else(duration < 26, 1, 0),
    dur_1yr = if_else(duration >= 26 & duration < 52, 1, 0),
    # dur_abv = if_else(duration >= 52, 1, 0),
    
    # Number of arms (2 - 3 - 4 - >4)
    arm_2 = if_else(n_arms == 2, 1, 0),
    arm_3 = if_else(n_arms == 3, 1, 0),
    # arm_4 = if_else(n_arms == 4, 1, 0),
    # arm_abv = if_else(n_arms > 4, 1, 0),
    
    # Phase (3)
    p3 = if_else(phase_3 == 1, 1, 0),
    # p4 = if_else(phase_4 == 1, 1, 0),
    
    # Blinding (double - triple - quadruple - open label)
    bdbl = double,
    btrip = triple,
    bquad = quadruple,
    # bol = open_label,
    
    # N facilities (<200)
    f_200 = if_else(n_facilities < 200, 1, 0),
    
    # N poms (1 - 2 - 3 - >3)
    pom_1 = if_else(n_poms == 1, 1, 0),
    pom_2 = if_else(n_poms == 2, 1, 0),
    # pom_3 = if_else(n_poms == 3, 1, 0),
    # pom_abv = if_else(n_poms >3, 1, 0),
    
    # Regions (am - eur - wpac - afr - sea - em)
    reg_am = region_of_the_americas,
    reg_eur = european_region,
    reg_wpac = western_pacific_region,
    reg_afr = african_region,
    reg_sea = south_east_asia_region,
    reg_em = eastern_mediterranean_region
    
  )

glimpse(a_dmy)

# Removed n_abv, arm_abv, p4, bol, pom_abv and reg_am
# High collinearity and multicollinearity, removed to serve as references
X_dummies <- a_dmy %>% select(n_250:reg_em)
X_mat <- as.matrix(X_dummies)
X_cor <- cor(X_mat, method = "pearson")
X_high_cor <- which(abs(X_cor) > 0.8 & lower.tri(X_cor), arr.ind = TRUE)
colnames(X_dummies)[unique(c(X_high_cor[,1], X_high_cor[,2]))]

# Check for multicollinearity
library(car)
X_r <- rnorm(nrow(X_dummies))
X_mdl <- lm(X_r ~ ., data = X_dummies)
alias(X_mdl) # Check for aliased covariates
X_vif <- vif(X_mdl)
X_vif[X_vif > 5] # None after commenting out above

kappa(X_mdl)

covs <- names(a_dmy[c(28:58)])

b_stan_dat <- list(
  I = nrow(a_dmy),
  K = length(covs),
  r = a_dmy$arm_n_attr,
  N = a_dmy$arm_n,
  dt = a_dmy$pt,
  arm = as.integer(a_dmy$trttype),
  X = as.matrix(a_dmy[, covs])
)

glimpse(b_stan_dat)

# Save
saveRDS(b_stan_dat, "processed_data/stan_dat.rds")
