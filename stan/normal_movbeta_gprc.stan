
data {
  int<lower=0> N; // number of data points, length(y)
  int<lower=0> K; // number of monthly lags
  vector[N] y;    // response
  matrix[N,K] X;  // matrix of climate covariates
}

transformed data {
  real month[K]; // vector of months, to create distance matrix
  for(k in 1:K) month[k] = k;
}

parameters {
  real alpha;
  real mu_beta;
  real<lower=0> sigma;
  real<lower=0> eta; // maximum covariance for betas
  real<lower=0> rho; // degree of temporal autocorrelation for betas
  vector[K] z;       // unit normal prior for non-centered term
}

transformed parameters {
  vector[N] yhat;
  vector[K] beta;
  matrix[K,K] Sigma; // covariance matrix
  matrix[K,K] L;     // cholesky of covariance matrix
  
  // covariance
  Sigma = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, K));
  L = cholesky_decompose(Sigma);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // linear predictor
  yhat = alpha + X * beta;
}

model {
  z ~ normal(0, 1);
  
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  sigma ~ normal(0, 3);

  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
