
# Splines ----------------------------------------------------------------------

# Safely select df for splines
safe_df <- function(x, df_min = 1, df_max = 3) {
  x <- as.numeric(x)
  # NA if there are less than or equal to 3 unique values
  # Or if length is less than 130 (1st quantile for trial size)
  if (length(unique(x)) <= 5) return(NA)
  
  for (df in seq(df_max, df_min, by = -1)) {
    valid <- tryCatch({
      suppressWarnings({
        # Fit covariate with splines and get attributes
        b <- ns(x, df = df)
        knots <- attr(b, "knots")
        boundary <- attr(b, "Boundary.knots")
        
        # Check if knots exist and are not at boundaries (min, max)
        knots_ok <- length(knots) > 0 &&
          any(knots > boundary[1] & knots < boundary[2]) &&
          !all(knots == boundary[1]) &&
          !all(knots == boundary[2])
        
        # Check if spline basis has full rank (i.e. independence between transf)
        full_rank <- is.matrix(b) && qr(b)$rank == ncol(b)
        
        knots_ok && full_rank
      })
    },
    error = function(e) FALSE
    )
    if (valid) return(df)
  }
  return(NA)
}
        
# safe_df <- function(x, df_min = 1, df_max = 3) {
#   x <- as.numeric(x)
#   if (length(unique(x)) <= 3) return(NA)
#   
#   for (df in seq(df_max, df_min, by = -1)) {
#     warning_triggered <- FALSE
#     valid <- tryCatch(
#       withCallingHandlers({
#         b <- ns(x, df = df)
#         knots <- attr(b, "knots")
#         boundary <- attr(b, "Boundary.knots")
#         
#         knots_ok <- (df == 1) || (
#           length(knots) > 0 &&
#             any(knots > boundary[1] & knots < boundary[2]) &&
#             !all(knots == boundary[1]) &&
#             !all(knots == boundary[2])
#         )
#         
#         full_rank <- is.matrix(b) && qr(b)$rank == ncol(b)
#         if (knots_ok && full_rank) TRUE else FALSE
#       },
#       warning = function(w) {
#         # Disqualify df if warning triggered
#         warning_triggered <<- TRUE
#         invokeRestart("muffleWarning")
#       }),
#       error = function(e) FALSE
#     )
#     
#     if (valid) return(df)
#   }
#   
#   return(NA)
# }

# Safely select df for splines
# safe_df <- function(x, df_min = 2, df_max = 3) {
#   last_valid_df <- NA
# 
#   for (df in df_min:df_max) {
#     warning_triggered <- TRUE
#     knot_issue <- TRUE
# 
#     tryCatch(
#       withCallingHandlers(
#         {
#           b <- ns(x, df = df)
#           knots <- attr(b, "knots")
#           boundary <- attr(b, "Boundary.knots")
# 
#           if (
#             length(knots) > 0 &&
#             !all(knots == boundary[1]) &&
#             !all(knots == boundary[2])
#           ) {
#             warning_triggered <- FALSE
#             knot_issue <- FALSE
#           }
#         },
#         warning = function(w) {
#           warning_triggered <<- TRUE
#           invokeRestart("muffleWarning")
#         }
#       ),
#       error = function(e) {
#         warning_triggered <<- TRUE
#       }
#     )
# 
#     if (!warning_triggered && !knot_issue) {
#       return(df)
#       break
#     }
#   }
# 
#   return(NA)
# }

# Apply FPs  -------------------------------------------------------------------

# Replace fp terms where estimates were unstable using FP
replace_fp <- function(fml_string, vars_to_fix, nct_id) {
  
  if (is.null(vars_to_fix)) {
    return(fml_string)
  }
  
  for (var in vars_to_fix) {
    if (grepl("df", fml_string)) {
      fml_string <- str_replace_all(
        fml_string,
        fixed(paste0("fp(", var, ", df = 2)")), var
      )
    } else {
      fml_string <- str_replace_all(
        fml_string,
        fixed(paste0("fp(", var, ")")), var
      )
    }
  }
  
  return(fml_string)
  
}

# Contrast computation ---------------------------------------------------------

# Build metadata for power tranformations
metadata_power <- function(res_fix, nct_id = NA) {
  res_fix %>% 
    filter(str_starts(term, "I")) %>% 
    mutate(
      # Take transformation string
      inner = str_match(term, "I\\((.*)\\)")[, 2],
      # Take term
      var = str_extract(term, b_terms),
      # Identify power transformation
      power = str_match_all(inner, "\\^([-]?[0-9.]+)") %>% map(~ as.numeric(.x[,2])),
      # Identify if it has a log-transformation
      has_log = str_detect(inner, "log"),
      # Stack powers where > 1, otherwise return n=1
      full_powers = map2(power, has_log, function(p, log_term) {
        if (length(p) == 1 && log_term) return(rep(p, 2)) else return(p)
      })
    ) %>% 
    group_by(var) %>% 
    summarise(
      power = list(sort(unique(unlist(full_powers)))),
      terms = list(term),
      .groups = "drop"
    ) %>% 
    mutate(
      power = map2(power, terms, function(p, t) {
        if (length(p) == 1 && length(t) == 2) rep(p, 2) else p
      }),
      nct_id = nct_id
    )
}

# Generate basis vector for FP (i.e. what it looks like at a given value
# representing a 1SD increase if scaled, or one-unit increase if unscaled)
fp_basis <- function(x, powers) {
  if (length(powers) == 1) return(c(x^powers[1]))
  if (powers[1] == powers[2]) return(c(x^powers[1], x*powers[1] * log(x)))
  return(c(x^powers[1], x^powers[2]))
}

# Compute contrasts
# Function for computing contrasts
compute_contrast <- function(beta, Sigma, powers, x0, x1) {
  f0 <- fp_basis(x0, powers)
  f1 <- fp_basis(x1, powers)
  f_diff <- f1 - f0
  contrast <- as.numeric(crossprod(f_diff, beta))
  variance <- as.numeric(t(f_diff) %*% Sigma %*% f_diff)
  se <- sqrt(variance)
  cov <- as.numeric(t(f_diff) %*% Sigma)
  tibble(contrast, se, f_diff, variance, cov)
}

# Apply compute_contrast
apply_comcon <- function(res, vmat, terms, powers, nct_id, var) {
  
  # Global x0 and x1 values
  x_row <- b_contrast_points %>% filter(nct_id == !!nct_id & var == !!var)
  
  if (nrow(x_row) != 1) {
    message("No unique value for x0/x1 for ", var, " in trial", nct_id)
    return(tibble(contrast = NA_real_, se = NA_real_, var = var))
  }
  
  x0_val <- x_row$x0
  x1_val <- x_row$x1
  
  # Give transformation names for vcov rows and names (ordered correctly)
  vmat_fixed <- vmat
  rownames(vmat_fixed)[1:nrow(res)] <- res$term
  colnames(vmat_fixed)[1:nrow(res)] <- res$term
  
  # Get betas for FP covariates
  beta <- res %>% 
    filter(term %in% terms) %>% 
    arrange(match(term, terms)) %>% 
    pull(estimate)
  
  # Get 2x2 vcov
  Sigma <- vmat_fixed[terms, terms]
  
  if (is.na(x0_val) || is.na(x1_val) || any(is.na(beta))) {
    return(tibble(contrast = NA_real_, se = NA_real_, var = var))
  }
  
  # Compute contrast for FP covariates
  compute_contrast(beta, Sigma, powers, x0 = x0_val, x1 = x1_val) %>%
    mutate(var = var)
  
}

# Reconstruct vcov with contrast variance-covariance
replace_vcov <- function(vcov, term, var, f_diff) {
  
  # Get names of non-fp variables (keeping)
  other_term <- setdiff(colnames(vcov), term)
  
  # Extract submatrices
  vcov_fp <- vcov[term, term, drop = FALSE]
  vcov_fp_other <- vcov[term, other_term, drop = FALSE]
  
  # Compute contrast variance and covariance
  contrast_var <- as.numeric(t(f_diff) %*% vcov_fp %*% f_diff)
  contrast_cov <- as.numeric(t(f_diff) %*% vcov_fp_other)
  names(contrast_cov) <- other_term
  
  # Drop FP term from matrix
  vcov_reduce <- vcov[other_term, other_term, drop = FALSE]
  
  # Build new matrix including contrast
  n <- length(other_term)
  vcov_new <- matrix(NA_real_, nrow = n + 1, ncol = n + 1)
  rownames(vcov_new) <- colnames(vcov_new) <-  c(other_term, var)
  
  # Fill in values
  vcov_new[other_term, other_term] <- vcov_reduce
  vcov_new[var, var] <- contrast_var
  vcov_new[var, other_term] <- contrast_cov
  vcov_new[other_term, var] <- contrast_cov
  
  return(vcov_new)
}

# Apply replace_vcov for all vars at once per trial
apply_all_contrasts <- function(vcov, contrast_data) {
  vcov_tmp <- vcov
  for (i in seq_len(nrow(contrast_data))) {
    vcov_tmp <- replace_vcov(
      vcov_tmp,
      contrast_data$term[[i]],
      contrast_data$var[[i]],
      contrast_data$f_diff[[i]]
    )
  }
  vcov_tmp
}


# Detect convergence issues and assign NULL for model --------------------------

safe_cox <- function(fml, df, id) {
  warning_msg <- NULL
  
  fit <- withCallingHandlers(
    expr = mfp(fml, family = "cox", data = df),
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("mufflewarning")
    }
  )
  
  if (!is.null(warning_msg)) {
    message("Convergence issue for ", id)
  }
  
  return(fit)
}
  
# Handle weak variables (few events/variability) -------------------------------
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

# Iterate over each covariate and fit cox and lr models ------------------------
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













