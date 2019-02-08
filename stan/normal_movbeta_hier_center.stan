
data {
  int<lower=0> N; // number of data points, length(y)
  int<lower=0> K; // number of monthly lags
  vector[N] y;    // response
  matrix[N,K] X;  // matrix of climate covariates
}

parameters {
  real alpha;
  real<lower=0> sigma;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[K] beta;
}

transformed parameters {
  vector[N] yhat;
  
  // linear predictor
  yhat = alpha + X * beta;
}

model {
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
