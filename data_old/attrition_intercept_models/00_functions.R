

# Handle weak variables (few events/variability)
# Causes convergence issues leaving them in
filter_weak_vars <- function(data,
                             id_var = "id",
                             time_var = "weeks",
                             status_var = "event",
                             min_events_per_level = 2,
                             min_unique_values = 3,
                             min_sd_threshold = 1e-3,
                             coerce_integer_to_factor = TRUE,
                             max_levels_for_coercion = 5) {
  
  core_vars <- c(id_var, time_var, status_var)
  outcome_vars <- data %>% select(c(all_of(core_vars)))
  covariates <- data %>% select(-all_of(core_vars))
  status <- outcome_vars[[status_var]]
  
  # Logic for which variables to keep
  keep_vars <- map_lgl(covariates, function(x) {
    
    if (all(is.na(x))) return(FALSE)
    
    if (is.numeric(x) || is.integer(x)) {
      unique_vals <- length(unique(na.omit(x)))
      if (unique_vals < min_unique_values) return(FALSE)
      if (sd(x, na.rm = TRUE) < min_sd_threshold) return(FALSE)
      
      if (is.integer(x) && coerce_integer_to_factor && unique_vals <= max_levels_for_coercion) {
        x <- factor(x)
        tab <- table(x, status, useNA = "no")
        if (!"1" %in% colnames(tab)) return(FALSE)
        event_counts <- tab[, "1", drop = FALSE]
        return(all(event_counts >= min_events_per_level))
      }
      
      return(TRUE)
    }
    
    if (is.factor(x) || is.character(x)) {
      tab <- table(x, status, useNA = "no")
      if (!"1" %in% colnames(tab)) return(FALSE)
      event_counts <- tab[, "1", drop = FALSE]
      return(all(event_counts >= min_events_per_level))
    }
    
    return(FALSE)
    
  })
  
  good_covariates <- covariates[which(keep_vars)]
  
  cleaned_data <- bind_cols(outcome_vars, good_covariates)
  
  return(cleaned_data)
}

# Function to iterate over each covariate and fit cox and lr models
# Saves outputs in a list
single_cov_mdls <- function(data, covariates) {
  
  mdl_res <- lapply(covariates, function(covariate) {
    
    # Create formula
    f_cox <- as.formula(paste0("Surv(weeks, event) ~ ", covariate))
    f_lr <- as.formula(paste0("event ~ ", covariate))
    
    # Fit model
    cox <- coxph(f_cox, data = data)
    lr <- glm(f_lr, data = data, family = "binomial")
    
    # Get outputs
    cox_res <- tidy(cox)[1, c("term", "estimate", "std.error")] %>% mutate(model = "cox")
    lr_res <- tidy(lr)[2, c("term", "estimate", "std.error")] %>% mutate(model = "lr")
    res <- bind_rows(cox_res, lr_res)
    
  })
  
  # Get outputs
  final_res <- bind_rows(mdl_res)
  return(final_res)
  
}













