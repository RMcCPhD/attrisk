
source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)
library(pROC)

# Import data (variable names, model outputs)
# a_imp_df <- colnames(readRDS("mars/design_model/latest_stan_design.rds")$X)
a_imp_df <- readRDS("mars/design_model/latest_stan_design.rds")
a_imp_res <- readRDS("mars/design_model/outputs/design_sum_draws.rds")
a_imp_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")
a_imp_rhat <- read_csv("mars/design_model/outputs/sum_rhat_ess.csv")

# Trace plots ------------------------------------------------------------------

# Trace plots for global parameters
trace_global <- traceplot(
  a_imp_mdl, 
  pars = c(
    "lp__", "mean_rate", "sd_rate", "mu_mu1", "tau_mu1", "tau_beta"
  )
)

ggsave(
  "outputs/design_model/ppc/trace_global.png", 
  trace_global, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots for covariate betas
trace_beta <- traceplot(
  a_imp_mdl, 
  pars = "beta"
)

ggsave(
  "outputs/design_model/ppc/trace_beta.png", 
  trace_beta, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots of Gompertz parameters - first 56 of shape and rate
trace_a <- traceplot(a_imp_mdl, pars = paste0("a[", 1:56, "]"))
trace_b <- traceplot(a_imp_mdl, pars = paste0("b[", 1:56, "]"))

ggsave(
  "outputs/design_model/ppc/trace_a.png", 
  trace_a, 
  width = 15, 
  height = 10,
  units = "in"
)

ggsave(
  "outputs/design_model/ppc/trace_b.png", 
  trace_b, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots of predicted attrition probabilities and counts (first 56)
trace_p <- traceplot(a_imp_mdl, pars = paste0("p[", 1:56, "]"))
trace_r <- traceplot(a_imp_mdl, pars = paste0("r_new[", 1:56, "]"))

ggsave(
  "outputs/design_model/ppc/trace_p.png", 
  trace_p, 
  width = 15, 
  height = 10,
  units = "in"
)

ggsave(
  "outputs/design_model/ppc/trace_r.png", 
  trace_r, 
  width = 15, 
  height = 10,
  units = "in"
)

# Rate coefficients ------------------------------------------------------------

betas <- a_imp_res %>% 
  filter(grepl("^beta\\[", variable)) %>% 
  separate(variable, into = c("beta", "order"), sep = "\\[") %>% 
  mutate(order = gsub("\\]", "", order)) %>% 
  cbind(data.frame(var = colnames(a_imp_df$X))) %>% 
  select(var, everything())

betas_plot <- betas %>% 
  mutate(
    var = case_match(
      var,
      "age_lower" ~ "Eligibility: Lower age limit",
      "age_upper" ~ "Eligibility: Upper age limit",
      "blindingdouble" ~ "Double blinding",
      "blindingquadruple" ~ "Quadruple blinding",
      "blindingtriple" ~ "Triple blinding",
      "duration" ~ "Total follow-up time (years)",
      "max_bmi" ~ "Eligibility: Upper BMI limit",
      "max_egfr" ~ "Eligibility: Upper eGFR limit",
      "min_bmi" ~ "Eligibility: Lower BMI limit",
      "min_egfr" ~ "Eligibility: Lower eGFR limit",
      "n_arms" ~ "Number of arms",
      "n_facilities" ~ "Number of facilities/sites",
      "n_log" ~ "Sample size",
      "one_pom" ~ "Single primary outcome",
      "phase" ~ "Phase",
      "reg_afr" ~ "Sites in African regions",
      "reg_am" ~ "Sites in the Americas",
      "reg_em" ~ "Sites in Eastern Medi.",
      "reg_eur" ~ "Sites in Europe",
      "reg_sea" ~ "Sites in SE Asia",
      "reg_wpac" ~ "Sites in Western Pacific",
      "trttype" ~ "Treatment type",
      .default = var
    ),
    var = factor(
      var,
      levels = c(
        "Eligibility: Upper age limit",
        "Eligibility: Lower age limit",
        "Eligibility: Upper BMI limit",
        "Eligibility: Lower BMI limit",
        "Eligibility: Upper eGFR limit",
        "Eligibility: Lower eGFR limit",
        "Total follow-up time (years)",
        "Treatment type",
        "Phase",
        "Sample size",
        "Number of arms",
        "Number of facilities/sites",
        "Single primary outcome",
        "Double blinding",
        "Quadruple blinding",
        "Triple blinding",
        "Sites in African regions",
        "Sites in the Americas",
        "Sites in Eastern Medi.",
        "Sites in Europe",
        "Sites in SE Asia",
        "Sites in Western Pacific"
      )
    )
  ) %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(var))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  scale_x_continuous(n.breaks = 15) +
  theme_classic() +
  labs(x = "Mean effect on Gompertz rate (95% credible interval)", y = NULL)

betas_plot

a_imp_res %>% 
  filter(grepl("^log_b\\[", variable)) %>% 
  separate(variable, into = c("log_rate", "order"), sep = "\\[") %>% 
  mutate(order = gsub("\\]", "", order)) %>% 
  reframe(mean = mean(mean))

ggsave(
  "outputs/design_model/ppc/rate_betas.png",
  betas_plot,
  width = 8,
  height = 8,
  units = "in"
)

# AUC --------------------------------------------------------------------------

# Observed variables (arms, counts)
obs <- tibble(
  iter = 1L:344L,
  r = a_imp_df$r,
  N = a_imp_df$N
)

# Get posterior draws for predicted attrition
pred_draws <- as_draws_df(a_imp_mdl)
pred_df <- pred_draws[, grepl("^p\\[", names(pred_draws)), drop = FALSE]

# Sample posterior draws
set.seed(123)
draw_samples <- sample.int(nrow(pred_df), size = min(2000, nrow(pred_df)), replace = FALSE)

# Build arm-level labels and scores from counts and probabilities
build_lbls_scores <- function(p_vec, r, N) {
  scores <- rep(p_vec, times = N)
  labels <- unlist(Map(function(ri, Ni) c(rep(1L, ri), rep(0L, Ni - ri)), r, N), use.names = FALSE)
  list(labels = labels, scores = scores)
}

# Compute roc and auc per posterior draw across arms
roc_auc_lst <- map(draw_samples, function(did) {
  p <- as.numeric(pred_df[did, ])
  ls <- build_lbls_scores(p, obs$r, obs$N)
  roc_obj <- roc(response = ls$labels, predictor = ls$scores, quiet = TRUE)
  list(roc = roc_obj, auc = as.numeric(auc(roc_obj)))
})

# AUC summary
auc_vals <- vapply(roc_auc_lst, `[[`, numeric(1), "auc")

auc_sum <- tibble(
  mean_auc = mean(auc_vals),
  sd_auc   = sd(auc_vals),
  med_auc  = median(auc_vals),
  min_auc  = min(auc_vals),
  max_auc  = max(auc_vals),
  lo_auc = quantile(auc_vals, 0.025),
  hi_auc = quantile(auc_vals, 0.975)
)

# Mean AUC 0.729 - moderate discrimination
auc_sum

write_csv(auc_sum, "outputs/design_model/valid/auc_sum.csv")

# Plot auc density
auc_dens <- ggplot(tibble(auc = auc_vals), aes(x = auc)) +
  geom_density() +
  theme_classic() +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "AUC", y = "Density")

ggsave(
  "outputs/design_model/valid/auc_dist.png",
  auc_dens,
  width = 8,
  height = 8,
  units = "in"
)

# Aggregate roc across draws
roc_df <- map_dfr(roc_auc_lst, function(x) {
  data.frame(
    fpr = 1 - x$roc$specificities,
    tpr = x$roc$sensitivities
  )
}, .id = "draw")

roc_sum <- roc_df %>%
  mutate(
    fpr_bin = cut(fpr, breaks = seq(0, 1, by = 0.01))
  ) %>% 
  group_by(fpr_bin) %>% 
  reframe(
    fpr = mean(as.numeric(gsub("\\((.+),.*", "\\1", fpr_bin))),
    median_tpr = median(tpr),
    q2.5 = quantile(tpr, 0.025),
    q97.5 = quantile(tpr, 0.975),
  )

# Plot roc - narrow confidence intervals indicative of accuracy
roc_plot <- roc_sum %>% 
  ggplot(aes(x = fpr, y = median_tpr, ymin = q2.5, ymax = q97.5)) +
  geom_line() +
  geom_ribbon(alpha = 0.3, fill = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_classic() +
  labs(
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)"
  )

ggsave(
  "outputs/design_model/valid/roc_plot.png",
  roc_plot,
  width = 8,
  height = 8,
  units = "in"
)

# Brier score ------------------------------------------------------------------

# Join observed and predicted attrition probabilities
# Add treatment type
join_p <- obs %>% 
  mutate(obs_p = r / N) %>% 
  left_join(
    pred_df %>% 
      pivot_longer(everything()) %>% 
      mutate(name = as.integer(str_extract(name, "\\d+"))) %>% 
      rename(iter = name, p = value)
  ) %>% 
  left_join(
    tibble(iter = 1L:344L, trttype = a_imp_df$X[,1])
  )

# Compute Brier score
brier <- join_p %>%
  group_by(iter) %>% 
  reframe(brier = mean((obs_p - p)^2))

# Repeat by arm
brier_arm <- join_p %>%
  group_by(iter, trttype) %>% 
  reframe(brier = mean((obs_p - p)^2))

# Summarise
brier_sum <- brier %>% 
  reframe(
    mean_brier = mean(brier),
    sd_brier = sd(brier),
    med_brier = median(brier),
    min_brier = min(brier),
    max_brier = max(brier),
    low_brier = quantile(brier, 0.025),
    hi_brier = quantile(brier, 0.975)
  )

brier_sum_arm <- brier_arm %>% 
  group_by(trttype) %>% 
  reframe(
    mean_brier = mean(brier),
    sd_brier = sd(brier),
    med_brier = median(brier),
    min_brier = min(brier),
    max_brier = max(brier),
    low_brier = quantile(brier, 0.025),
    hi_brier = quantile(brier, 0.975)
  )

write_csv(brier_sum, "outputs/design_model/valid/brier_sum.csv")
write_csv(brier_sum_arm, "outputs/design_model/valid/brier_sum_arm.csv")

# Plot density
brier_plot <- brier %>% 
  ggplot(aes(x = brier)) +
  geom_density() +
  theme_classic() +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "Brier score", y = "Density")

brier_arm_plot <- brier_arm %>% 
  mutate(trttype = if_else(trttype == 0L, "Placebo", "Active")) %>% 
  ggplot(aes(x = brier, fill = factor(trttype, levels = c("Placebo", "Active")))) +
  geom_density(colour = "black", alpha = 0.3) +
  theme_classic() +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "Brier score", y = "Density", fill = "Treatment type")

ggsave(
  "outputs/design_model/valid/brier_dist.png",
  brier_plot,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  "outputs/design_model/valid/brier_arm_dist.png",
  brier_arm_plot,
  width = 8,
  height = 8,
  units = "in"
)
