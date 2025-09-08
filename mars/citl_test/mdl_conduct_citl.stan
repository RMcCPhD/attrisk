// Clamp p to prevent exact 0 or 1s, allows STAN to sample in stable regions
functions {
  real clamp_p(real x) {
    return fmin(fmax(x, 1e-12), 1 - 1e-12);
  }
}

data {
  // Non-CITL inputs
  int<lower=1> I;                       // number of data points (arms in trials)
  int<lower=1> K;                       // number of covariates
  array[I] int<lower=0> r;              // observed attrition counts
  array[I] int<lower=1> N;              // number of participants
  vector[I] dt;                         // total follow-up time in years
  vector[I] lp;                         // covariate linear predictor for shape
  real beta_rate;                       // rate beta from linear model of shape
  real<lower=0> resid_sd;               // residual SD from linear model of shape
  matrix[I, K] X;                       // covariate matrix
  
  // CITL inputs
  int<lower=1> J;                       // Participant per arm
  int<lower=1> pK;                      // Number of sobol covariates
  matrix[I * J, pK] Z;                  // Stacked matrix of Sobol sequences
  vector[pK] pB;                        // Pooled Cox log-HRs
}

transformed data {
  // Participant-level linear predictor across arms
  vector[I * J] hr_all = exp(Z * pB);   // Participant-level HR
  
  // Helper - starting row for each arm
  array[I] int arm_start;
  for (i in 1:I) arm_start[i] = (i - 1) * J + 1;
}

parameters {
  // Hierarchical structure for rate
  real mu_mu1;                          // global log-mean for intercept (mu1)
  real<lower=0> tau_mu1;                // sd for log(mu1)
  vector[I] z_mu1;                      // standard normal latent mu1 variables
  
  // Global-local shrinkage for design features
  real<lower=0> tau_beta;               // global shrinkage factor
  vector<lower=0>[K] lambda_beta;       // local shrinkage factor
  vector[K] z_beta;                     // latent normal betas
  
  // Error for shape
  vector[I] error;                      // N(0, resid_sd)
  
  // Arm-specific CITL intercept shift (log scale)
  real gamma_citl;                      // Global CITL factor
  vector[I] gamma_i_raw;                // Arm-specific calibration
  real<lower=0> gamma_sigma;            // Shrinkage scale
}

transformed parameters {
  
  // Centre gamma
  vector[I] gamma_i = gamma_i_raw - mean(gamma_i_raw);
  
  // Arm-specific log-rate intercepts and betas
  vector[I] mu1 = mu_mu1 + tau_mu1 * z_mu1;
  vector[K] beta = tau_beta * (lambda_beta .* z_beta);
  
  // Calibrated rate
  vector[I] log_b_cal = mu1 + X * beta + gamma_citl + gamma_i;
  vector[I] b_cal = exp(log_b_cal);      
  
  // Shape updated with calibrated log_rate
  vector[I] a_cal = log_b_cal * beta_rate + lp + error;
  
  // Compute average arm risk from Sobol sequences
  vector[I] p;
  for (i in 1:I) {
    int s = arm_start[i]; // Starting row index
    vector[J] hr_i = segment(hr_all, s, J);  // Slice HR rows
    
    // Safely compute Gompertz constant
    real a_i = a_cal[i];
    real b_i = b_cal[i];
    real t_i = dt[i];
    real C_i;
    
    if (fabs(a_i) < 1e-4) {
      // Limit case: Reduce to exp if a ≈ 0
      C_i = b_i * t_i;
    } else {
      C_i = (b_i / a_i) * (exp(a_i * t_i) - 1);
    }
    
    real p_i = 1 - mean(exp(-hr_i * C_i));
    p[i] = clamp_p(p_i);
  }
}

model {
  // b has a log-normal prior via exp(mu1 + X * beta)
  mu_mu1 ~ normal(0, 1);
  tau_mu1 ~ normal(0, 1);               // truncated normal via <lower=0>
  z_mu1 ~ normal(0, 1);
  
  // Hierarchical global-local gaussian prior for beta shrinkage
  tau_beta ~ normal(0, 1);
  lambda_beta ~ normal(0, 1);
  z_beta ~ normal(0, 1);
  
  // Normal prior for error, scale is residual sd from linear model of shape
  error ~ normal(0, resid_sd);
  
  // Normal priors for CITL gamma
  gamma_citl ~ normal(0, 1);
  gamma_i_raw ~ normal(0, gamma_sigma);
  gamma_sigma ~ normal(0, 1);
  
  // Likelihood using average risk from sobol sequences
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;                   // predicted attrition counts
  vector[I] log_lik;                    // log-likelihood
  
  // Pre and post-CITL intercepts
  vector[I] log_b_pre  = mu1 + X * beta;
  vector[I] log_b_post = log_b_cal;

  // Expected versus observed attrition counts and rates
  vector[I] exp_count;                  // N[i] * p[i]
  vector[I] exp_rate;                   // p[i]
  vector[I] obs_count;                  // r[i]
  vector[I] obs_rate;                   // r[i] / N[i]

  // Monitoring CITL
  real expected_total = 0;
  int  observed_total = 0;
  int  total_N = 0;

  for (i in 1:I) {
    // Posterior predictive and log-likelihood
    r_new[i]   = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);

    // per-arm diagnostics
    exp_count[i] = N[i] * p[i];
    exp_rate[i]  = p[i];
    obs_count[i] = r[i];
    obs_rate[i]  = r[i] / N[i];

    // accumulate totals
    expected_total += exp_count[i];
    observed_total += r[i];
    total_N        += N[i];
  }

  // Overall CITL summaries
  real avg_expected = expected_total / total_N;
  real avg_observed = observed_total / total_N;

  // (Optional) expose calibrated parameters for downstream use
  vector[I] b_cal_out = b_cal;
  vector[I] a_cal_out = a_cal;
}

