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
  array[I] int<lower=0> r;              // observed attrition counts
  array[I] int<lower=0> N;              // number of participants
  vector[I] dt;                         // total follow-up time in years
}

parameters {
  real mu_mu1;                          // global centre for log(mu1)
  real<lower=0> tau_mu1;                // global scale for log(mu1)
  vector[I] z_mu1;                      // standard normal latent mu1
  
  real mu_a;                            // fixed shape (a)
}

transformed parameters {
  vector[I] mu1;                        // intercept for each arm per trial (row)
  vector[I] a;                          // a for each row
  vector[I] b;                          // b for each row
  vector[I] p;                          // p for each row
  
  mu1 = mu_mu1 + tau_mu1 * z_mu1;       // compute mu1 from hyperparameters
  a = mu_a + tau_a * z_a;               // shared shape
  
  for (i in 1:I) { 
    b[i] = exp(mu1[i]);                 // Iteratively compute b
  }
  
  p = gompertz_prob(b, a, dt);          // safely compute p (see function)
}

model {
  // Priors for log-rate
  mu_mu1 ~ normal(-1.42, 0.54);         // Intercept and MAD from rate ~ shape
  tau_mu1 ~ normal(0.54, 1.08);         // MAD and MAD*2 from rate ~ shape
  z_mu1 ~ normal(0, 1);                 // Standard latent normal
  
  // Priors for log-rate
  mu_a ~ normal(-3.85, 1.31);           // Fixed shape
  
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