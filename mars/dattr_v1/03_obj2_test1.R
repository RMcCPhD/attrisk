
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    class_short = case_match(
      class_short,
      "a_gluc" ~ "agluc",
      "biguanides" ~ "biguanide",
      .default = class_short
    ),
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Import IPD model exports (coefficients, variance-covariance)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Unadjusted treatment (test) --------------------------------------------------

# Setup aggregate data
b_agg_df <- a_imp_df %>% 
  filter(source != "ipd") %>% 
  select(trial_id, class_short, atc_short, ref_class:arm_attr) %>% 
  mutate(
    atc_short = if_else(class_short == "oad", "A10BX", atc_short)
  )

# Setup ipd data
b_ipd_filter <- a_imp_ipd_res %>% 
  filter(
    spec == "unadj",
    modeltype == "glm"
  ) %>% 
  rename(class_short = term) %>% 
  select(-c(std.error:uci)) %>% 
  left_join(
    b_agg_df %>% 
      select(class_short, atc_short) %>% 
      distinct()
  )

# Add placebo rows with NA estimate for reference level
b_ipd_df <- b_ipd_filter %>%
  filter(!grepl("Inter", class_short)) %>% 
  bind_rows(
    b_ipd_filter %>% 
      group_by(nct_id) %>% 
      transmute(
        nct_id,
        modeltype,
        spec,
        class_short = "placebo",
        atc_short = "placebo",
        ref,
        estimate = NA
      )
  ) %>% 
  arrange(nct_id) %>% 
  distinct()

# Prepare vcov, setting NA for intercepts
b_vcov <- a_imp_vcov %>% 
  filter(modeltype == "glm_unadj") %>% 
  mutate(
    cov2 = map(cov_matrix, function(mat) {
      
      rownames(mat) <- gsub("Intercept", "placebo", rownames(mat))
      rownames(mat) <- gsub("\\(|\\)", "", rownames(mat))
      colnames(mat) <- gsub("Intercept", "placebo", colnames(mat))
      colnames(mat) <- gsub("\\(|\\)", "", colnames(mat))
      
      placebo_cell <- which(rownames(mat) == "placebo")
      
      if (length(placebo_cell) > 0) {
        mat[placebo_cell, ] <- NA
        mat[, placebo_cell] <- NA
      }
      
      return(mat)
      
    })
  )

b_vcov_ready <- list(b_vcov$cov2)
b_vcov_ready <- b_vcov_ready[[1]]
names(b_vcov_ready) <- unique(b_vcov$nct_id)

# Keep non-NA cells
b_vcov_prep <- map2(b_vcov_ready, names(b_vcov_ready), function(mat, id) {
  
  est_terms <- b_ipd_df %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(class_short)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

# Setup networks
c_agg <- set_agd_arm(
  data = b_agg_df,
  study = trial_id,
  trt = class_short,
  trt_ref = "placebo",
  r = arm_attr,
  n = arm_n,
  trt_class = atc_short
)
  
plot(c_agg)
  
c_ipd <- set_agd_regression(
  data = b_ipd_df,
  study = nct_id,
  trt = class_short,
  trt_ref = "placebo",
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt,
  trt_class = atc_short
)

plot(c_ipd)

# Combine
c_cmbn <- combine_network(c_agg, c_ipd)
plot(c_cmbn)

# Add integration
d_mdl <- nma(
  c_cmbn,
  trt_effects = "random",
  regression = ~ .trt,
  class_interactions = "common",
  prior_intercept = normal(scale = 10),
  prior_reg = normal(scale = 10), chains = 4, cores = 4
)
