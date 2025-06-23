functions {
  vector gompertz_prob(vector b, vector a, vector dt) {
    int I = rows(b);
    vector[I] p;
    for (i in 1:I) {
      if (fabs(a[i]) < 1e-4) {
        p[i] = 1 - exp(-b[i] * dt[i]);
      } else {
        p[i] = 1 - exp(-b[i] / a[i] * (exp(a[i] * dt[i]) - 1));
      }
    }
    return p;
  }
}

data {
  int<lower=0> I;                    // number of arms
  array[I] int<lower=0> r;           // observed attrition
  array[I] int<lower=0> N;           // total enrolled
  vector[I] dt;                      // follow-up time
}

parameters {
  vector[I] z_raw;                   // latent z (noncentered)
  real<lower=0> s;                   // scale of z
  real mu_mu1;                       // mean log-rate

  real beta_0;                       // shape intercept
  real beta_1;                       // shape slope
}

transformed parameters {
  vector[I] z;                       // scaled latent
  vector[I] mu1;                     // log-rate
  vector[I] a;                       // shape
  vector[I] b;                       // rate
  vector[I] p;                       // attrition probability
  
  z = s * z_raw;                     // scaled latent

  for (i in 1:I) {
    mu1[i] = mu_mu1 + z[i];                 // log-rate
    a[i]   = beta_0 + beta_1 * mu1[i];      // shape from learned regression
    b[i]   = exp(mu1[i]);                   // rate
  }

  p = gompertz_prob(b, a, dt);             // compute probabilities
}

model {
  // Empirical priors
  mu_mu1 ~ normal(-1.05, 0.08);            // log-rate intercept
  s ~ normal(0, 0.5);                      // scale of latent z
  z_raw ~ normal(0, 1);                    // latent trial deviation

  beta_0 ~ normal(-3.85, 0.5);             // shape intercept
  beta_1 ~ normal(-0.21, 0.1);             // shape slope on log-rate

  r ~ binomial(N, p);                      // binomial likelihood
}

generated quantities {
  array[I] int r_new;
  vector[I] log_lik;
  vector[I] z_raw_out;
  vector[I] z_out;
  vector[I] s_out;

  for (i in 1:I) {
    r_new[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(r[i] | N[i], p[i]);
    z_raw_out[i] = z_raw[i];
    z_out[i] = z[i];
    s_out[i] = s;
  }

  // Posterior summaries
  real mean_rate = mean(b);
  real sd_rate = sd(b);
  real mean_shape = mean(a);
  real sd_shape = sd(a);
}
