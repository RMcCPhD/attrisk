functions {
  // Safely use Gompertz cumulative distribution
  // If a is close to 0, collapse to exponential
  // Otherwise use Gompertz
  vector gompertz_prob(vector b, vector a, vector dt) {
    int I = rows(b);
    vector[I] p;
    for (i in 1:I) {
      if (fabs(a[i]) < 1e-4) {
        // Limit case: Exponential model when a ≈ 0
        p[i] = 1 - exp(-b[i] * dt[i]);
      } else {
        // General Gompertz model
        p[i] = 1 - exp(-b[i] / a[i] * (exp(a[i] * dt[i]) - 1));
      }
    }
    return p;
  }
  
  // Safely obtain effective rate for a given target p, shape and time
  // Use to obtain log-diff between effective and actual log-rate
  // Delta produced, added to log-rate to calibrate intercept to participant
  // risk profile per arm
  real get_eff_rate(real p_tgt, real a, real dt) {
    real p = fmin(1 - 1e-12, fmax(1e-12, p_tgt)); // Clamped to stabilise
    // Inverted Gompertz CDF
    if (fabs(a) < 1e-4) {
      // Reduce to inverted exp when a ≈ 0
      return -log1m(p) / fmax(1e-12, dt);
    } else {
      // Inverted Gompertz CDF
      real num = a * (-log1m(p));
      real den = fmax(1e-12, exp(a * dt) - 1);
      return num / den;
    }
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
  int<lower=1> pK;                      // Participant-level covariates
  matrix[I * J, pK] Z;                  // Stacked Sobol sequences for covariates
  vector[pK] pB;                        // Pooled Cox log-HRs
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
}

transformed parameters {
  // (Pre-CITL) Arm-specific log-rate intercepts and betas
  vector[I] mu1 = mu_mu1 + tau_mu1 * z_mu1;
  vector[K] beta = tau_beta * (lambda_beta .* z_beta);
  
  // (Pre-CITL) Arm-specific baseline log-rate and rate
  vector[I] log_b_base = mu1 + X * beta;
  vector[I] b_base = exp(log_b_base);      
  
  // Arm-specific shape (conditioned on log-rate, linear predictor and error)
  // parameters from shape model fit externally
  vector[I] a_base = log_b_base * beta_rate + lp + error;
  
  // Participant-level linear predictor across arms
  vector[I * J] eta_all = Z * pB;       // Stacked across arms
  vector[I] avg_p;                      // Average participant risk per arm
  
  // Compute arm-specific average risk
  for (i in 1:I) {
    int s = (i - 1) * J + 1;                              // Starting row index
    vector[J] eta_i = segment(eta_all, s, J);             // Slice participant row
    vector[J] hr_i  = exp(eta_i);                         // Get participant HR
    vector[J] b_i   = rep_vector(b_base[i], J) .* hr_i;   // Participant rate
    vector[J] a_i   = rep_vector(a_base[i], J);           // Participant shape
    vector[J] dt_i  = rep_vector(dt[i], J);               // Duration
    vector[J] p_i   = gompertz_prob(b_i, a_i, dt_i);      // Compute risk
    avg_p[i]        = mean(p_i);                          // Average arm risk
  }
  
  // Invert to effective baseline rate (rate without participant risk profiles)
  vector[I] log_b_eff;
  for (i in 1:I) {
    real b_eff_i = fmax(get_eff_rate(avg_p[i], a_base[i], dt[i]), 1e-12);
    log_b_eff[i] = log(b_eff_i);
  }
  
  // Apply delta to intercept, obtain updated log-rate
  vector[I] delta     = log_b_eff - log_b_base; // arm-specific shift on log-rate
  vector[I] mu1_cal   = mu1 + delta;            // updated intercept
  vector[I] log_b_cal = mu1_cal + X * beta;     // construct log-rate
  vector[I] b_cal     = exp(log_b_cal);         // new rate

  // 5) Calibrated probabilities used in likelihood
  vector[I] p = gompertz_prob(b_cal, a_base, dt);
  
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
  
  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;                   // predicted attrition counts
  vector[I] log_lik;                    // log-likelihood
  
  // Calibrated intercept components
  vector[I] log_b_base_out = log_b_base;
  vector[I] log_b_cal_out  = log_b_cal;
  vector[I] delta_out      = delta;
  vector[I] avg_risk_out   = avg_p;
  
  // Rate marginals (pre-CITL and calibrated)
  real mean_rate = mean(b_base);
  real sd_rate = sd(b_base);
  real mean_rate_cal = mean(b_cal);
  real sd_rate_cal = sd(b_cal);

  // Posterior predictive and log-likelihood
  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

}
