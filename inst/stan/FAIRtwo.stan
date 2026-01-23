// FAIR(2) model 

data {
  int<lower=0> N;          // Number of observations
  vector[N] Y;             // Dependent variable
  vector[N] Z;             // Narrative shock
  int<lower=1> H;          // Horizon
  
  // Priors
  real prior_a1; real prior_b1; real prior_c1;
  real prior_a2; real prior_b2; real prior_c2;
  real<lower=0> sigma_scale;
  real<lower=0> sdpar;
}

transformed data {
  matrix[N - H, H + 1] Z_lags;
  vector[N - H] Y_trim;
  
  for (t in 1:(N - H)) {
    Y_trim[t] = Y[t + H]; 
    for (h in 0:H) {
      Z_lags[t, h + 1] = Z[t + H - h];
    }
  }
}

parameters {
  real a1;
  real b1;
  real<lower=0> c1;
  real a2;
  real<lower=b1> b2; // Constraint: Peak 2 must be after Peak 1
  real<lower=0> c2;
  real<lower=0> sigma;
}

transformed parameters {
  vector[H + 1] psi;
  
  for (h in 0:H) {
    psi[h + 1] = a1 * exp(-square((h - b1) / c1) / 2.0) + 
                 a2 * exp(-square((h - b2) / c2) / 2.0);
  }
}

model {
  // --- Priors ---
  a1 ~ normal(prior_a1, sigma_scale);
  b1 ~ normal(prior_b1, sdpar);
  c1 ~ normal(prior_c1, sdpar);
  
  a2 ~ normal(prior_a2, sigma_scale);
  b2 ~ normal(prior_b2, sdpar);
  c2 ~ normal(prior_c2, sdpar);
  
  sigma ~ exponential(1 / sigma_scale);
  
  // --- Likelihood (Vectorized) ---
  Y_trim ~ normal(Z_lags * psi, sigma);

}
