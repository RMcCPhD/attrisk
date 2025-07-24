// Safely estimate Gompertz cumulative distribution
// If a is close to 0, collapse to exponential
// Otherwise use Gompertz
functions {
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
}

data {
  int<lower=0> I;                       // number of data points (arms in trials)
  int<lower=0> K;                       // number of covariates
  array[I] int<lower=0> r;              // observed attrition counts
  array[I] int<lower=0> N;              // number of participants
  vector[I] dt;                         // total follow-up time in years
  vector[I] lp;                         // covariate linear predictor for shape
  real beta_rate;                       // rate beta from linear model of shape
  real resid_sd;                        // residual SD from linear model of shape
  matrix[I, K] X;                       // covariate matrix
}

parameters {
  real mu_mu1;                          // global log-mean for intercept (mu1)
  real<lower=0> tau_mu1;                // sd for log(mu1)
  vector[I] z_mu1;                      // standard normal latent mu1 variables
  
  real<lower=0> tau_beta;               // global shrinkage factor
  vector<lower=0>[K] lambda_beta;       // local shrinkage factor
  vector[K] z_beta;                     // latent normal betas
}

transformed parameters {
  // log-rate intercept per arm|trial (row)
  vector[I] mu1 = mu_mu1 + tau_mu1 * z_mu1;
  
  // beta with shrinkage
  vector[K] beta = tau_beta * lambda_beta .* z_beta;
  
  // b for each row (exp of mu1)
  vector[I] b = exp(mu1 + X * beta);      
  
  // a for each row (conditioned on b and covariates)
  vector[I] a = b * beta_rate + lp + resid_sd;       
  
  // p for each row
  vector[I] p = gompertz_prob(b, a, dt);                          
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
  
  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;                   // predicted attrition counts
  vector[I] log_lik;                    // log-likelihood

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

  real mean_rate = mean(b);             // mean rate
  real sd_rate = sd(b);                 // sd rate
}
