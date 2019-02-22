
data {
  int N;         // number of data points, length(y)
  int K;         // number of monthly lags
  vector[N] y;   // response
  matrix[N,K] X; // matrix of climate covariates
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  simplex[K] theta;
}

transformed parameters {
  vector[K] beta_wt;
  vector[N] yhat;

  beta_wt = theta * beta;
  
  yhat = alpha + X * beta_wt;
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
  theta ~ dirichlet(rep_vector(1.0, K));
 
  // model
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
