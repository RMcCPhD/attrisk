
# Pool distributions for each variable based on best-fitting
# Helper functions written to pool for each trial
# Flexsurv models fitted where not normal to find best-fitting distribution,
# pooled used respective parameters and converted to mean and sd

source("scripts/00_config.R")
library(tidyverse)

# Import imputation iterations
a_imp_lst <- readRDS("processed_data/agg_bl_imp_lst.rds")
a_imp_df <- readRDS("processed_data/agg_bl_design_imputed_iters.rds")

# Check distributions of variable marginals
plot_dist <- a_imp_df %>% 
  pivot_longer(age_mu:former_never) %>% 
  mutate(
    name = case_when(
      grepl("asian|black|other|white", name) ~ paste0("race_", name),
      grepl("curr|form", name) ~ paste0("smoke_", name),
      TRUE ~ name
    )
  ) %>% 
  ggplot(aes(x = value, group = imp)) +
  geom_density(alpha = 0.01) +
  facet_wrap(~name, scales = "free")

plot_dist

ggsave(
  "outputs/conduct_model/agg_bl/dist_before_clp.png",
  height = 10,
  width = 10,
  units = "in"
)

# Find best-fitting distribution per variable ----------------------------------

# Parametric distributions
z_dists <- c("norm", "lnorm", "gamma", "beta", "weibull")

# Function to fit distribution to each variable and get best
find_best_dist <- function(x, families = z_dists) {
  x <- x[is.finite(x)]
  fits <- list()
  
  # Conditionally apply
  
  if ("norm" %in% families) {
    fits$norm  <- fitdistrplus::fitdist(x, "norm")
  }
  
  if ("lnorm" %in% families && all(x > 0)) {
    fits$lnorm <- fitdistrplus::fitdist(x, "lnorm")
  }
  
  if ("gamma" %in% families && all(x > 0)) {
    fits$gamma <- fitdistrplus::fitdist(x, "gamma")
  }
  
  if ("beta" %in% families && all(x >= 0 & x <= 1)) {
    # Avoid boundary issues (i.e. exactly 0 or 1)
    eps <- 1e-6
    x2  <- pmin(pmax(x, eps), 1 - eps)
    fits$beta <- fitdistrplus::fitdist(x2, "beta")
  }
  
  if ("weibull" %in% families && all(x > 0)) {
    fits$weibull <- fitdistrplus::fitdist(x, "weibull")
  }
  
  # Get AICs
  aics <- map_dbl(fits, ~ .x$aic)
  fam <- names(aics)[which.min(aics)]
  list(family = fam, fit = fits[[fam]])
}

# Apply function to variables
z_vars <- names(a_imp_df[5:27])

b_best_dists <- map(z_vars, function(v) {
  best <- find_best_dist(a_imp_df[[v]])
  tibble(
    var = v,
    family = best$family,
    aic = best$aic,
    estimates = list(if (is.null(best$fit)) NA else best$fit$estimate)
  )
})

b_best_df <- bind_rows(b_best_dists) %>% 
  separate(estimates, into = c("est1", "est2"), sep = ", ") %>% 
  mutate(across(est1:est2, ~ gsub("c\\(|\\(|\\)", "", .))) %>% 
  pivot_longer(est1:est2) %>% 
  separate(value, into = c("dist", "value"), sep = " = ") %>% 
  select(-name)







