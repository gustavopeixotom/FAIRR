// Exponential decay basis function for FAIR - 1.1.0

data {
  int<lower=0> N;          // Number of observations
  vector[N] Y;             // Dependent variable
  vector[N] Z;             // shock
  int<lower=1> H;          // Horizon
  
  // Priors
  real prior_a;            // Scale prior
  real prior_l;            // Decay rate prior
  real<lower=0> sigma_scale; // Standard deviation prior for scale
  real<lower=0> sdpar;    // Standard deviation prior for decay rate
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

// Parameters 
parameters {
  real a; // Scale
  real<lower=0> l; // Decay MINIMUM 0
  real<lower=0> sigma; // Error
}

// Transformed parameters
transformed parameters {
  vector[H + 1] psi; 
  
  for (h in 0:H) {
    psi[h + 1] = a * exp(-l * h); // Exponential decay function
  }
}

// Model
model {
  // Priors
  a ~ normal(prior_a, sigma_scale);
  l ~ normal(prior_l, sdpar);
  
  sigma ~ exponential(1 / sigma_scale);
  
  // Likelihood (Vectorized)
  Y_trim ~ normal(Z_lags * psi, sigma);

}
