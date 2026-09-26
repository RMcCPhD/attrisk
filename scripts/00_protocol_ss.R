
# Exploratory assessment of the number of candidate predictor parameters
# supported by the available conduct-stage IPD

rm(list = ls())
gc()

library(tidyverse)
library(pmsampsize)

# Helper function for survival-outcome calculations ---------------------------

pms_surv <- function(csrsquared,
                     parameters,
                     rate,
                     timepoint,
                     meanfup,
                     shrinkage = 0.90) {
  
  out <- pmsampsize(
    type = "s",
    csrsquared = csrsquared,
    parameters = parameters,
    rate = rate,
    timepoint = timepoint,
    meanfup = meanfup,
    shrinkage = shrinkage
  )
  
  tibble(
    csrsquared = csrsquared,
    parameters = parameters,
    rate = rate,
    timepoint = timepoint,
    meanfup = meanfup,
    shrinkage = shrinkage,
    sample_size = out$sample_size,
    events = out$events,
    epp = out$EPP
  )
  
}

# pmsampsize calculates the minimum sample size for a specified number of
# candidate predictor parameters. This function inverts that calculation by
# finding the largest integer number of parameters whose required sample size
# does not exceed the available sample.
find_max_parameters <- function(available_sample,
                                csrsquared,
                                rate,
                                timepoint,
                                meanfup,
                                shrinkage = 0.90) {
  
  required_n <- function(parameters) {
    pms_surv(
      csrsquared = csrsquared,
      parameters = parameters,
      rate = rate,
      timepoint = timepoint,
      meanfup = meanfup,
      shrinkage = shrinkage
    )
  }
  
  # First identify an upper bound that is not feasible
  lower <- 0L
  upper <- 1L
  upper_res <- required_n(upper)
  
  while (upper_res$sample_size <= available_sample) {
    lower <- upper
    upper <- upper * 2L
    upper_res <- required_n(upper)
  }
  
  # Then use a binary search to locate the exact feasibility boundary
  while ((upper - lower) > 1L) {
    midpoint <- as.integer(floor((lower + upper) / 2))
    midpoint_res <- required_n(midpoint)
    
    if (midpoint_res$sample_size <= available_sample) {
      lower <- midpoint
    } else {
      upper <- midpoint
    }
  }
  
  max_res <- required_n(lower)
  next_res <- required_n(upper)
  
  tibble(
    available_sample = available_sample,
    csrsquared = csrsquared,
    rate = rate,
    timepoint = timepoint,
    meanfup = meanfup,
    shrinkage = shrinkage,
    max_parameters = lower,
    required_sample_at_max = max_res$sample_size,
    expected_events_at_max = max_res$events,
    epp_at_max = max_res$epp,
    expected_events_in_available_sample = available_sample * rate * meanfup,
    expected_epp_in_available_sample =
      (available_sample * rate * meanfup) / lower,
    first_infeasible_parameters = upper,
    required_sample_first_infeasible = next_res$sample_size
  )
  
}

# Conduct-stage model assumptions ---------------------------------------------

# Available IPD
b_ipd_n <- 19023L
b_ipd_events <- 4022L

# Assumptions retained from the original exploratory assessment:
# - one-year prediction horizon and mean follow-up
# - event rate of 0.107 per participant-year, based on the previous attrition
#   rates analysis of trials involving participants with type 2 diabetes
# - anticipated Cox-Snell R-squared of 0.10
# - target shrinkage of 0.90 (the pmsampsize default)
b_ipd_rate <- 0.107
b_ipd_fu <- 1
b_ipd_tmpt <- 1
b_ipd_csr2 <- 0.10
b_ipd_shrinkage <- 0.90

# Begin with the available sample and calculate the largest supportable number 
# of candidate predictor parameters
b_ipd_primary <- find_max_parameters(
  available_sample = b_ipd_n,
  csrsquared = b_ipd_csr2,
  rate = b_ipd_rate,
  timepoint = b_ipd_tmpt,
  meanfup = b_ipd_fu,
  shrinkage = b_ipd_shrinkage
) %>%
  mutate(observed_events = b_ipd_events)

# Retain the two rows surrounding the primary feasibility boundary so that the
# maximum supported value and the first unsupported value can be inspected
b_ipd_boundary <- tibble(
  parameters = c(
    b_ipd_primary$max_parameters,
    b_ipd_primary$first_infeasible_parameters
  )
) %>%
  mutate(
    result = map(
      parameters,
      \(parameters) pms_surv(
        csrsquared = b_ipd_csr2,
        parameters = parameters,
        rate = b_ipd_rate,
        timepoint = b_ipd_tmpt,
        meanfup = b_ipd_fu,
        shrinkage = b_ipd_shrinkage
      )
    )
  ) %>%
  select(-parameters) %>%
  unnest(result) %>%
  mutate(
    available_sample = b_ipd_n,
    feasible = sample_size <= available_sample,
    sample_difference = available_sample - sample_size,
    observed_events = b_ipd_events,
    observed_epp = observed_events / parameters
  )

# Secondary check for the candidate set considered during model development.
# pmsampsize counts candidate predictor parameters rather than raw predictor
# variables; update this value if the 13 characteristics require more than 13
# parameters after accounting for categorical coding or non-linear terms.
b_ipd_planned_parameters <- 13L

b_ipd_planned <- pms_surv(
  csrsquared = b_ipd_csr2,
  parameters = b_ipd_planned_parameters,
  rate = b_ipd_rate,
  timepoint = b_ipd_tmpt,
  meanfup = b_ipd_fu,
  shrinkage = b_ipd_shrinkage
) %>%
  mutate(
    available_sample = b_ipd_n,
    feasible = sample_size <= available_sample,
    sample_difference = available_sample - sample_size,
    observed_events = b_ipd_events,
    observed_epp = observed_events / parameters
  )

# Sensitivity analyses ---------------------------------------------------------

# Explore how the maximum supported number of parameters changes across the
# Cox-Snell R-squared and event-rate scenarios considered previously.
b_ipd_scenarios <- expand_grid(
  csrsquared = seq(0.05, 0.20, by = 0.01),
  rate = c(b_ipd_rate, 0.10, 0.15, 0.20) %>% unique() %>% sort(),
  timepoint = b_ipd_tmpt,
  meanfup = b_ipd_fu
)

b_ipd_sensitivity <- pmap_dfr(
  b_ipd_scenarios,
  \(csrsquared, rate, timepoint, meanfup) {
    find_max_parameters(
      available_sample = b_ipd_n,
      csrsquared = csrsquared,
      rate = rate,
      timepoint = timepoint,
      meanfup = meanfup,
      shrinkage = b_ipd_shrinkage
    )
  }
) %>%
  arrange(rate, csrsquared)

b_ipd_sensitivity_summary <- b_ipd_sensitivity %>%
  summarise(
    min_supported_parameters = min(max_parameters),
    max_supported_parameters = max(max_parameters),
    min_epp_at_boundary = min(epp_at_max),
    max_epp_at_boundary = max(epp_at_max)
  )

# Export results ---------------------------------------------------------------

write_csv(
  b_ipd_sensitivity,
  "outputs/pmsampsize/pmsampsize_ipd_all.csv"
)

write_csv(
  b_ipd_primary,
  "outputs/pmsampsize/pmsampsize_ipd_best.csv"
)

write_csv(
  b_ipd_boundary,
  "outputs/pmsampsize/pmsampsize_ipd_boundary.csv"
)

write_csv(
  b_ipd_planned,
  "outputs/pmsampsize/pmsampsize_ipd_planned.csv"
)