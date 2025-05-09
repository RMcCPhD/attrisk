functions {
  vector gompertz_prob(vector b, vector a, vector dt) {
    int I = rows(b);
    vector[I] p;
    for (i in 1:I) {
      if (fabs(a[i]) < 1e-5) {
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
  int<lower=0> I;                       // number of data points
  array[I] int<lower=0> r;              // observed counts
  array[I] int<lower=0> N;              // number of participants
  vector[I] dt;                         // total follow-up time
}

parameters {
  vector[I] mu1;                        // intercepts per trial
  real mu_a;                            // global log-mean for a
  real<lower=0> tau_a;                  // sd for log(a)
  vector[I] z_a;                        // standard normal latent a variables
}

transformed parameters {
  vector[I] a;
  vector[I] b;
  vector[I] p;

  a = exp(mu_a + tau_a * z_a);          // allow a to vary around zero
  b = exp(mu1);                         // No covariate contribution
  p = gompertz_prob(b, a, dt);          // Safely compute p with limit handling
}

model {
  // Priors
  mu1 ~ normal(0, 1);
  mu_a ~ normal(0, 1);
  tau_a ~ normal(0, 1);                 // Truncated automatically by <lower=0>
  z_a ~ normal(0, 1);

  // Likelihood
  r ~ binomial(N, p);
}

generated quantities {
  array[I] int r_new;
  vector[I] p_new = p;
  vector[I] a_out = a;
  vector[I] b_out = b;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }

  real rate = mean(b);
}
