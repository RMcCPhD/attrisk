
rm(list = ls())
gc()

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Config
options(scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Import stan data and flattened arrays
a_imp_dat <- readRDS("latest_stan_citl.rds")
a_imp_mu1 <- readRDS("citl_mu1.rds")
a_imp_beta <- readRDS("citl_beta.rds")
a_imp_error <- readRDS("citl_error.rds")

# Get dims
total_draws <- nrow(a_imp_mu1)
I <- ncol(a_imp_mu1)
K <- ncol(a_imp_beta)

# Get subset of posterior samples
D_target <- 120L
n_chains <- 4L
iters <- total_draws / n_chains
per_chain <- max(1L, floor(D_target / n_chains))
chain_idx <- unique(as.integer(round(seq(1, iters, length.out = per_chain))))
ids <- as.integer(unlist(lapply(0:(n_chains - 1L), function(c) c * iters + chain_idx)))
ids <- ids[ids >= 1 & ids <= total_draws]
ids <- unique(ids)
message(sprintf("Selected %d stage-1 draws (%d per chain).", length(ids), per_chain))

# Compile stan model
b_compile <- stan_model("mdl_conduct_citl.stan")

# Build data with the pre-CITL draws per trial
build_stage2 <- function(d) {
  within(a_imp_dat, {
    mu1_draw   <- a_imp_mu1[d, ]    # length I
    beta_draw  <- a_imp_beta[d, ]   # length K
    error_draw <- a_imp_error[d, ]  # length I
  })
}

# Iterate model per trial
run_batch <- function(draws_ids, out_dir = "citl_fits") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (d in draws_ids) {
    dat <- build_stage2(d)
    
    fit <- sampling(
      b_compile,
      data = dat,
      warmup = 250,
      iter = 500,
      chains = 4,
      seed = 123,
      control = list(adapt_delta = 0.95, max_treedepth = 12)
    )
    
    saveRDS(fit, file.path("citl_fits", sprintf("citl_fit_d%05d.rds", d)))
  }
}

run_batch(ids)