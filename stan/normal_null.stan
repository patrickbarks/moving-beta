
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] yhat;
  for(i in 1:N) yhat[i] = alpha;
}

model {
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 3);
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
