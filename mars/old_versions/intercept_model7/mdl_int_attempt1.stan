// Conditional use of Gompertz cumulative distribution function
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
  vector[I] z_raw;                      // raw latent parameter
  real<lower=0> s;                      // learned scale of latent variation
  
  real mu_mu1;
  real mu_shape;
}

transformed parameters {
  vector[I] z = z_raw;                  // standardised latent parameters
  vector[I] mu1;                        // log-rate for each arm
  vector[I] a;                          // shape for each arm
  vector[I] b;                          // rate for each arm
  vector[I] p;                          // outcome probability for each arm
  
  for (i in 1:I) {
    mu1[i] = mu_mu1 + s * 0.55 * z[i];  
    a[i]   = mu_shape + s * (-1.63) * z[i];   
    b[i]   = exp(mu1[i]);             
  }

  p = gompertz_prob(b, a, dt);          // CDF-based attrition probability
}

model {
  // Empirical priors on intercepts
  mu_mu1 ~ normal(-1.05, 0.08);
  mu_shape ~ normal(-1.79, 0.26);
  
  // Standard prior on z
  z_raw ~ normal(0, 1);  
  
  // Weakly informative, positive-only prior on learned scale
  s ~ normal(0, 0.5);
  
  r ~ binomial(N, p);                  
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
  }

  // Posterior summaries
  real mean_rate = mean(b);
  real sd_rate = sd(b);
  real mean_shape = mean(a);
  real sd_shape = sd(a);
}

