
# Approximate mean and SD where data given were categorical levels

# NCT03689374 eGFR -------------------------------------------------------------

# Normal, Mild, Moderate, Severe, End-stage
counts <- c(normal = 549, mild = 272, moderate = 52, severe = 1, end = 0)

approx_stats_from_egfr_groups <- function(counts, cap_top = 120, add_within_bin = TRUE) {
  stopifnot(length(counts) == 5)
  N <- sum(counts)
  breaks <- c(0, 15, 30, 60, 90, cap_top) # Define bins
  mids   <- c((90 + cap_top)/2, (60 + 89)/2, (30 + 59)/2, (15 + 29)/2, (0 + 14)/2)
  widths <- c(cap_top - 90,     89 - 60 + 1, 59 - 30 + 1, 29 - 15 + 1, 14 - 0 + 1)
  widths <- c(cap_top - 90, 30, 30, 15, 15)
  mids   <- mids[c(1, 2, 3, 4, 5)]
  widths <- widths[c(1, 2, 3, 4, 5)]
  
  # Weighted mean using midpoints
  mu <- sum(counts * mids) / N
  
  # Weighted SD using midpoints only
  sum_wx2 <- sum(counts * mids^2)
  sd_mid  <- sqrt((sum_wx2 - N * mu^2) / (N - 1))
  
  if (!add_within_bin) {
    tibble(cap_top = cap_top, mean = mu, sd = sd_mid)
  } else {
    # Add within-bin variance (uniform assumption: Var = width^2 / 12)
    within_var_sum <- sum(counts * (widths^2) / 12)
    sum_wx2_adj    <- sum(counts * mids^2) + within_var_sum
    sd_adj         <- sqrt((sum_wx2_adj - N * mu^2) / (N - 1))
    tibble(cap_top = cap_top, mean = mu, sd_midpoint_only = sd_mid, sd_within_added = sd_adj)
  }
}

# Approximate
approx_stats_from_egfr_groups(counts, cap_top = 120, add_within_bin = TRUE)
