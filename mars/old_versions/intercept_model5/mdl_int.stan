# Conditional use of Gompertz cumulative distribution function
functions {
  vector gompertz_prob(vector b, vector a, vector dt) {
    int I = rows(b);
    vector[I] p;
    for (i in 1:I) {
      if (fabs(a[i]) < 1e-4) {
        // Reduce to exponential when absolute difference in a is <0.0001
        p[i] = 1 - exp(-b[i] * dt[i]);
      } else {
        // Otherwise use Gompertz CDF
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
  // Hyperparameters for hierarchical priors
  // Log-normal rate prior inferred from hyperparameters
  real mu_mu1;                          // centre for log-rate
  real<lower=0> tau_mu1;                // scale for log-rate
  vector[I] z_mu1;                      // standard normal latent log-rate
  
  // Normal shape prior inferred from hyperparameters
  real mu_a;                            // centre for a
  real<lower=0> tau_a;                  // scale for a
  vector[I] z_a;                        // standard normal latent a
}

transformed parameters {
  vector[I] mu1;                        // log-rate for each arm per trial
  vector[I] a;                          // a for each row
  vector[I] b;                          // b for each row
  vector[I] p;                          // p for each row
  
  mu1 = mu_mu1 + tau_mu1 * z_mu1;       // compute log-rate from hyperparameters
  a = mu_a + tau_a * z_a;               // compute a from hyperparameters
  
  for (i in 1:I) { 
    b[i] = exp(mu1[i]);                 // Iteratively compute b (arm i in 1:trial)
  }
  
  p = gompertz_prob(b, a, dt);          // safely compute p (see function)
}

model {
  // Priors for log-rate
  mu_mu1 ~ normal(-1.42, 0.54);         // Centre = Intercept from shape ~ rate
  // Scale = Median absolute deviance (MAD) from rate ~ shape
  
  tau_mu1 ~ normal(0.54, 1.08);         // Centre = MAD, Scale = MAD*2 for flexibility
  z_mu1 ~ normal(0, 1);                 // Standard latent normal
  
  // Priors for log-rate
  mu_a ~ normal(-3.85, 1.31);           // Centre = Intercept from rate ~ shape
  tau_a ~ normal(1.31, 2.62);           // MAD and MAD*2 from shape ~ rate
  z_a ~ normal(0, 1);                   // Standard latent normal
  
  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  // Predicted attrition counts
  array[I] int r_new;
  // Log-likelihood
  vector[I] log_lik;
  
  // Estimate attrition counts and log-likelihood for each arm in trials
  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }
  
  // Compute mean and sd rate
  real mean_rate = mean(b);
  real sd_rate = sd(b);
} 