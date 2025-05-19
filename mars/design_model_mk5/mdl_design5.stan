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
  matrix[I, K] X;                       // covariate matrix
}

parameters {
  real mu_mu1;                          // global log-mean for intercept (mu1)
  real<lower=0> tau_mu1;                // sd for log(mu1)
  vector[I] z_mu1;                      // standard normal latent mu1 variables
  
  real mu_a;                            // global log-mean for a
  real<lower=0> tau_a;                  // sd for log(a)
  vector[I] z_a;                        // standard normal latent a variables
  
  real<lower=0> tau_beta;               // global shrinkage factor
  
  // local shrinkage factors (a and b)
  vector<lower=0>[K] lambda_beta_a;
  vector<lower=0>[K] lambda_beta_b;
  
  // latent betas (a and b)
  vector[K] z_beta_a;
  vector[K] z_beta_b;                     
}

transformed parameters {
  // log-rate intercept for each arm/trial
  vector[I] mu1 = mu_mu1 + tau_mu1 * z_mu1;                        
  vector[I] log_a_raw = mu_a + tau_a * z_a; // log-a for each arm/trial 
  vector[I] a;                          // a for each arm/trial
  vector[I] b;                          // b for each arm/trial 
  vector[I] p;                          // p for each arm/trial 
  
  // Shrinkage applied to betas
  vector[K] beta_a = z_beta_a .* (tau_beta * lambda_beta_a);
  vector[K] beta_b = z_beta_b .* (tau_beta * lambda_beta_b);
  
  for (i in 1:I) { // Iteratively compute a and b
    a[i] = log_a_raw[i] + dot_product(X[i], beta_a);
    b[i] = exp(mu1[i] + dot_product(X[i], beta_b));
  }
  
  p = gompertz_prob(b, a, dt);          // safely compute p (see function)
}

model {
  // Priors
  // b has a log-normal prior via exp(mu1 + X * beta)
  mu_mu1 ~ normal(0, 1);
  tau_mu1 ~ normal(0, 1);               // Truncated automatically by <lower=0>
  z_mu1 ~ normal(0, 1);
  
  mu_a ~ normal(0, 1);
  tau_a ~ normal(0, 1);                 // Truncated automatically by <lower=0>
  z_a ~ normal(0, 1);
  
  // Shrinkage
  tau_beta ~ normal(0, 1);
  lambda_beta_a ~ normal(0, 1);
  lambda_beta_b ~ normal(0, 1);
  z_beta_a ~ normal(0, 1);
  z_beta_b ~ normal(0, 1);

  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

  real mean_rate = mean(b);
  real sd_rate = sd(b);
}
