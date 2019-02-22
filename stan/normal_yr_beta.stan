
data {
  int N;         // number of data points, length(y)
  int J;         // number of years
  vector[N] y;   // response
  matrix[N,J] X; // matrix of climate covariates
}

parameters {
  real alpha;
  vector[J] beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] yhat;
  yhat = alpha + X * beta;
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
