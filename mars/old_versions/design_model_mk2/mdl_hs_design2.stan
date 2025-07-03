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
  vector[I] arm;                        // treatment arm (0 = plc, 1 = ac)
  matrix[I, K] X;                       // covariate matrix
}

parameters {
  real mu_mu1;                          // global log-mean for intercept (mu1)
  real<lower=0> tau_mu1;                // sd for log(mu1)
  vector[I] z_mu1;                      // standard normal latent mu1 variables
  
  real mu_a;                            // global log-mean for a
  real<lower=0> tau_a;                  // sd for log(a)
  vector[I] z_a;                        // standard normal latent a variables
  
  real beta_arm;                        // arm coefficient
  
  // horseshoe prior (prev. gaussian global-local)
  real<lower=0> tau_beta;               // global shrinkage factor
  vector<lower=0>[K] lambda_beta;       // local shrinkage factor
  vector[K] z_beta;                     // latent betas
}

transformed parameters {
  vector[I] mu1;                        // intercept for each arm per trial (row)
  vector[I] a;                          // a for each row
  vector[I] b;                          // b for each row
  vector[I] p;                          // p for each row
  vector[K] beta;
  
  mu1 = mu_mu1 + tau_mu1 * z_mu1;       // compute mu1
  a = mu_a + tau_a * z_a;               // compute a on real scale
  beta = z_beta .* (tau_beta * lambda_beta); // shrinkage applied to beta
  
  for (i in 1:I) { // Iteratively compute b
    b[i] = exp(mu1[i] + arm[i] * beta_arm + dot_product(X[i], beta));
  }
  
  p = gompertz_prob(b, a, dt);          // safely compute p (see function)
}

model {
  // Priors
  // b has a log-normal prior via exp(mu1 + arm * beta_arm + X * beta)
  mu_mu1 ~ normal(0, 1);
  tau_mu1 ~ normal(0, 1);               
  z_mu1 ~ normal(0, 1);
  
  mu_a ~ normal(0, 1);
  tau_a ~ normal(0, 1);                 
  z_a ~ normal(0, 1);
  
  beta_arm ~ normal(0, 1);
  tau_beta ~ cauchy(0, 1);              // half-Cauchy for hyperparameters
  lambda_beta ~ cauchy(0, 1);
  z_beta ~ normal(0, 1);                // latent normal betas

  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;
  vector[I] a_out = a;
  vector[I] b_out = b;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

  real mean_rate = mean(b);
  real sd_rate = sd(b);
}
