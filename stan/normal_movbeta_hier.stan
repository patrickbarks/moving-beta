
data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N,K] X;
}

parameters {
  real alpha;
  real<lower=0> sigma;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[K] z;   // unit normal prior for non-centered term
}

transformed parameters {
  vector[N] yhat;
  vector[K] beta;
  
  // non-centered parameterization
  beta = mu_beta + sigma_beta * z;
  
  // linear predictor
  yhat = alpha + X * beta;
}

model {
  z ~ normal(0, 1);
  
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
