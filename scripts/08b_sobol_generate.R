
# Generate Sobol sequences representing synthetic IPD participant sets
# Sequences generated for each AgD trial

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import prepared Sobol materials
imp_dist <- readRDS("processed_data/sobol_dist.rds")
imp_cor <- readRDS("processed_data/sobol_cor_mat.rds")
imp_best <- read_csv("scratch_data/sobol_covariate_dist.csv")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import pooled log-hazard estimates from IPD models
imp_loghr <- readRDS("outputs/conduct/ma/out/ma_tidy.rds")

# Import covariate scaling
imp_scale <- read_csv("vivli/v2/Summaries/covariate_mu_d.csv")

# Helper: Centre variables on mean and scale by standard deviation
scale_centre_cov <- function(dat, imp_scale) {
  
  scale_lkp <- imp_scale %>% filter(name %in% names(dat))
  mu <- setNames(scale_lkp$mu, scale_lkp$name)
  d <- setNames(scale_lkp$d, scale_lkp$name)
  
  dat %>% 
    mutate(
      across(
        all_of(scale_lkp$name),
        ~ (.x - mu[cur_column()]) / d[cur_column()]
      )
    )
}


# Prepare distribution parameters ----------------------------------------------

# Check counts
imp_best

# Helper functions
get_lnorm_meanlog <- function(mu, d) {
  log(mu^2 / sqrt(d^2 + mu^2))
}

get_lnorm_sdlog <- function(mu, d) {
  sqrt(log(1 + d^2 / mu^2))
}

get_gamma_shape <- function(mu, d) {
  (mu / d)^2
}

get_gamma_rate <- function(mu, d) {
  mu / d^2
}

get_weibull_shape <- function(mu, d) {
  
  cv <- d / mu
  
  map_dbl(cv, ~ {
    uniroot(
      function(k) {
        sqrt(gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1) - .x
      },
      interval = c(0.1, 50)
    )$root
  })
}

get_weibull_scale <- function(mu, shape) {
  mu / gamma(1 + 1 / shape)
}

# Add parameters
prep_dist <- imp_dist %>% 
  mutate(
    across(ends_with("_d"), ~ pmax(.x, .Machine$double.eps)),
    
    # log-normal: bmi, chol_hdl, sbp, tot_hdl_ratio
    bmi_meanlog = get_lnorm_meanlog(bmi_mu, bmi_d),
    bmi_sdlog   = get_lnorm_sdlog(bmi_mu, bmi_d),
    
    chol_hdl_meanlog = get_lnorm_meanlog(chol_hdl_mu, chol_hdl_d),
    chol_hdl_sdlog   = get_lnorm_sdlog(chol_hdl_mu, chol_hdl_d),
    
    sbp_meanlog = get_lnorm_meanlog(sbp_mu, sbp_d),
    sbp_sdlog   = get_lnorm_sdlog(sbp_mu, sbp_d),
    
    # gamma: chol_total
    chol_total_shape = get_gamma_shape(chol_total_mu, chol_total_d),
    chol_total_rate  = get_gamma_rate(chol_total_mu, chol_total_d),
    
    # Weibull: durn, egfr
    durn_shape = get_weibull_shape(durn_mu, durn_d),
    durn_scale = get_weibull_scale(durn_mu, durn_shape),
    
    egfr_shape = get_weibull_shape(egfr_mu, egfr_d),
    egfr_scale = get_weibull_scale(egfr_mu, egfr_shape)
  )

# Generate Sobol sequences -----------------------------------------------------

# Nest AgD trials, generated 64 Sobol sequences per trial
# Join sobol data (distribution summaries, correlation matrices)
sobol_agd <- imp_agd %>% 
  group_by(trial_id) %>% 
  nest(.key = "agd") %>% 
  cbind(prep_dist %>% nest(.key = "dist")) %>% 
  cbind(imp_cor %>% select(cor = out))

# Generate 64 Sobol sequences per AgD trial
# Unnest Sobol points
sobol_gen <- sobol_agd %>% 
  mutate(
    sobol = map2(dist, cor, ~ {
      add_integration(
        .x,
        age           = distr(qnorm,    mean = .x$age_mu, sd = .x$age_d),
        female_p      = distr(qbern,    prob = .x$female_p),
        bmi           = distr(qlnorm,   meanlog = .x$bmi_meanlog, sdlog = .x$bmi_sdlog),
        durn          = distr(qweibull, shape = .x$durn_shape, scale = .x$durn_scale),
        sbp           = distr(qlnorm,   meanlog = .x$sbp_meanlog, sdlog = .x$sbp_sdlog),
        chol_hdl      = distr(qlnorm,   meanlog = .x$chol_hdl_meanlog, sdlog = .x$chol_hdl_sdlog),
        chol_total    = distr(qgamma,   shape = .x$chol_total_shape, rate = .x$chol_total_rate),
        hba1c         = distr(qnorm,    mean = .x$hba1c_mu, sd = .x$hba1c_d),
        egfr          = distr(qweibull, shape = .x$egfr_shape, scale = .x$egfr_scale),
        cursmoke      = distr(qbern,    prob = .x$cursmoke),
        asian         = distr(qbern,    prob = .x$asian),
        other         = distr(qbern,    prob = .x$other),
        black_aa      = distr(qbern,    prob = .x$black_aa),
        comorbs       = distr(qpois,    lambda = .x$comorbs_mu),
        conmeds       = distr(qpois,    lambda = .x$conmeds_mu),
        cor = .y,
        n_int = 64
      )
    }),
    sobol_unnest = map(sobol, ~ {
      .x %>% 
        select(starts_with(".int")) %>% 
        rename_with(., ~ gsub(".int_", "", .)) %>% 
        unnest(cols = everything())
    })
  )

# Save outputs
saveRDS(sobol_gen, "processed_data/sobol_out.rds")

# Transform sobol points - centre on mean, scale by standard deviation
# Summaries for transformation exported from Vivli, across all trials
sobol_unnest <- sobol_gen %>% 
  mutate(
    sobol_rdy = map(sobol_unnest, ~ {
      .x %>% 
        rename(
          age10 = age,
          sex = female_p
        ) %>%
        mutate(age10 = age10 / 10) %>%
        scale_centre_cov(imp_scale)
    })
  )

# Inspect distributions
sobol_plot <- sobol_unnest %>% 
  select(trial_id, sobol_rdy) %>% 
  unnest(sobol_rdy) %>% 
  ungroup() %>% 
  pivot_longer(-trial_id) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free") +
  theme_classic()

ggsave(
  "outputs/conduct/sobol/sobol_distributions.png",
  width = 8,
  height = 6,
  units = "in"
)

# Prepare CITL inputs ----------------------------------------------------------

citl_ma <- imp_loghr %>% 
  select(term, loghr) %>% 
  pivot_wider(names_from = "term", values_from = "loghr") %>% 
  rename_with(., ~ paste0(., "_loghr"))

# Rename pooled betas
citl_beta <- c(
  age10       = citl_ma$age10_loghr,
  bmi         = citl_ma$bmi_loghr,
  chol_hdl    = citl_ma$chol_hdl_loghr,
  chol_total  = citl_ma$chol_total_loghr,
  cursmoke    = citl_ma$smoke_loghr,
  durn        = citl_ma$durn_loghr,
  egfr        = citl_ma$egfr_loghr,
  hba1c       = citl_ma$hba1c_loghr,
  asian       = citl_ma$raceasian_loghr,
  black_aa    = citl_ma$raceblack_aa_loghr,
  other       = citl_ma$raceother_na_unk_loghr,
  sbp         = citl_ma$sbp_loghr,
  sex         = citl_ma$sex_loghr,
  comorbs     = citl_ma$comorbs_loghr,
  conmeds     = citl_ma$conmeds_loghr
)

# Parameters
I <- nrow(imp_agd)    # Number of trial arms
J <- 64L              # Number of Sobol points
Pnm <- names(citl_ma) # Predictor names
pK <- length(Pnm)     # Number of predictors

# Construct a matrix of sobol sequences per trial (IxJ) x pK
Z_blocks <- vector("list", I)

for (i in seq_len(I)) {
  trial <- imp_agd$trial_id[i]
  idx <- match(trial, sobol_unnest$trial_id)
  M <- sobol_unnest$sobol_rdy[[idx]]
  Z_blocks[[i]] <- as.matrix(M)
}

Z <- do.call(rbind, Z_blocks)

# List of outputs
sobol_lst <- list(
  I = I,
  J = J,
  Pnm = Pnm,
  pK = pK,
  Z = Z,
  loghr = citl_beta
)

str(sobol_lst)

# Save
saveRDS(sobol_lst, "processed_data/sobol_lst.rds")
