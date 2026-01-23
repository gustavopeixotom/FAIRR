// FAIR(1) model 

data {
  int<lower=0> N;          // Number of observations
  vector[N] Y;             // Dependent variable
  vector[N] Z;             // shock
  int<lower=1> H;          // Horizon
  
  // Priors
  real prior_a;
  real prior_b;
  real prior_c;
  real<lower=0> sigma_scale;
  real<lower=0> sdpar;
}

transformed data {
  // Pre-compute the Lag Matrix for fast matrix multiplication
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
  real a;
  real b;
  real<lower=0> c;
  real<lower=0> sigma;
}

transformed parameters {
  vector[H + 1] psi; 
  
  for (h in 0:H) {
    psi[h + 1] = a * exp(-square((h - b) / c) / 2.0);
  }
}

model {
  // VAR priors
  a ~ normal(prior_a, sigma_scale);
  b ~ normal(prior_b, sdpar);
  c ~ normal(prior_c, sdpar);
  
  sigma ~ exponential(1 / sigma_scale);
  
  // Likelihood (Vectorized)
  Y_trim ~ normal(Z_lags * psi, sigma);

}
