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
  cholesky_factor_corr[2] L;            // fixed Cholesky correlation factor
}

parameters {
  vector[2] mu;                         // mean of log-rate and shape
  vector<lower=0>[2] tau;               // sd of log-rate and shape
  matrix[I, 2] z;                       // latent uncorrelated standard normal
}

transformed parameters {
  vector[I] mu1;                        // log-rate for each arm per trial
  vector[I] a;                          // a for each row
  vector[I] b;                          // b for each row
  vector[I] p;                          // p for each row
  
  matrix[I, 2] theta_raw = z * L';
  matrix[I, 2] theta = rep_matrix(mu', I) + theta_raw .* rep_matrix(tau', I);
  
  mu1 = theta[, 1];
  a = theta[, 2];
  b = exp(mu1);
  p = gompertz_prob(b, a, dt);          // safely compute p (see function)
}

model {
  // Priors
  mu[1] ~ normal(-1.05, 0.73);          // mean and sd log-rate
  mu[2] ~ normal(-1.79, 2.26);          // mean and sd shape
  
  tau[1] ~ normal(0.73, 1.46);          // sd and sd*2 log-rate
  tau[2] ~ normal(2.26, 4.52);          // sd and sd*2 shape
  
  to_vector(z) ~ normal(0, 1);          // independent standard normal
  
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
  
  // Compute mean and sd rate and shape
  real mean_rate = mean(b);
  real sd_rate = sd(b);
  real mean_shape = mean(a);
  real sd_shape = sd(a);
} 
