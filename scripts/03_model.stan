
data {
  int<lower=0> I; // Number of trials
  int<lower=0> r[I]; // Observed attrition events per trial
  int<lower=0> N[I]; // Number of participants per trial
  vector[I] dt; // Time in weeks per trial

  // Covariates
  vector[I] arm; // Treatment arm (0 = placebo, 1 = active)
  vector[I] var1; // Phase 3 (0 = No, 1 = Yes, same for remaining variables)
  vector[I] var2; // Phase 4 
  vector[I] var3; // Other phase (i.e. not applicable given) 
  vector[I] var4; // Double blinding 
  vector[I] var5; // Open label 
  vector[I] var6; // Quadruple blinding 
  vector[I] var7; // Triple blinding 
  vector[I] var8; // 1 region 
  vector[I] var9; // 2 region 
  vector[I] var10; // 3 region 
  vector[I] var11; // 4 region 
  vector[I] var12; // 5 region 
  vector[I] var13; // 6 region 
  vector[I] var14; // In region of the americas 
  vector[I] var15; // In European region 
  vector[I] var16; // In Western-pacific region 
  vector[I] var17; // In African region 
  vector[I] var18; // In South-East Asian region 
  vector[I] var19; // In Eastern-Mediterranean region 
  vector[I] var20; // <250 sites 
  vector[I] var21; // 250 - <500 sites 
  vector[I] var22; // 500 - <750 sites 
  vector[I] var23; // 750 - <1000 sites 
  vector[I] var24; // >= 1000 sites 
  vector[I] var25; // Two arms 
  vector[I] var26; // Three arms 
  vector[I] var27; // >= 3 arms 
  vector[I] var28; // Duration < 6 months 
  vector[I] var29; // Duration 6 months - < 1 year 
  vector[I] var30; // Duration >= 1 year 
  vector[I] var31; // 1 POMs 
  vector[I] var32; // 2 POMs 
  vector[I] var33; // 3 POMs 
  vector[I] var34; // 4 POMs 
  vector[I] var35; // 5 POMs 
  vector[I] var36; // 6 POMs 
  vector[I] var37; // < 50 sites
  vector[I] var38; // 50 - <100 sites
  vector[I] var39; // 100 - <150 sites
  vector[I] var40; // >= 150 sites
}

parameters {
  real mu1; // Intercept
  real a_raw; // Raw shape parameter

  real beta_arm; // Treatment effect
  vector[40] beta; // Covariate betas
  real<lower=0> tau; // Shrinkage factor
}

transformed parameters {
  // Safely transform a if near 0
  real a;
  real epsilon = 1e-4;
  a = (a_raw < 0) ? a_raw - epsilon : (a_raw > 0 ? a_raw + epsilon : 0);

  vector[I] b; // Rate parameter
  vector[I] p; // Probability (cumulative incidence of attrition)

  for (i in 1:I) {
    b[i] = exp(
      mu1 +
      beta_arm * arm[i] +
      beta[1] * var1[i] +  beta[2] * var2[i] +  beta[3] * var3[i] +  beta[4] * var4[i] +
      beta[5] * var5[i] +  beta[6] * var6[i] +  beta[7] * var7[i] +  beta[8] * var8[i] +
      beta[9] * var9[i] +  beta[10] * var10[i] + beta[11] * var11[i] + beta[12] * var12[i] +
      beta[13] * var13[i] + beta[14] * var14[i] + beta[15] * var15[i] + beta[16] * var16[i] +
      beta[17] * var17[i] + beta[18] * var18[i] + beta[19] * var19[i] + beta[20] * var20[i] +
      beta[21] * var21[i] + beta[22] * var22[i] + beta[23] * var23[i] + beta[24] * var24[i] +
      beta[25] * var25[i] + beta[26] * var26[i] + beta[27] * var27[i] + beta[28] * var28[i] +
      beta[29] * var29[i] + beta[30] * var30[i] + beta[31] * var31[i] + beta[32] * var32[i] +
      beta[33] * var33[i] + beta[34] * var34[i] + beta[35] * var35[i] + beta[36] * var36[i] +
      beta[37] * var37[i] + beta[38] * var38[i] + beta[39] * var39[i] + beta[40] * var40[i]
    );
    p[i] = 1 - exp(-b[i] * (expm1(a * dt[i]) / a));
  }
}

model {
  mu1 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  beta_arm ~ normal(0, 1);
  beta ~ normal(0, tau); // Shared shrinkage scale
  tau ~ cauchy(0, 0.5); // Half-cauchy

  for (i in 1:I) {
    r[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities {
  vector[I] r_new;
  vector[I] p_new;

  for (i in 1:I) {
    p_new[i] = 1 - exp(-b[i] * (expm1(a * dt[i]) / a));
    r_new[i] = binomial_rng(N[i], p_new[i]);
  }

  real rate = exp(mu1); // Baseline hazard
}


