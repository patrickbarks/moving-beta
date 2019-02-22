functions {
  real dnorm(real x, real mu, real sigma) {
    real x1 = (1 / sqrt(2 *pi()*pow(sigma, 2)));
    real x2 = exp(-((x-mu)^2) / (2*pow(sigma, 2)));
    return(x1 * x2);
  }
  
  vector dnorm_vec(int n, real mu, real sigma) {
    vector[n] x;
    for (i in 1:n) { x[i] = dnorm(i, mu, sigma); }
    return(x / sum(x));
  }
}

data {
  int N;         // number of data points, length(y)
  int K;         // number of monthly lags
  vector[N] y;   // response
  matrix[N,K] X; // matrix of climate covariates
}

parameters {
  real<lower=0, upper=K> sens_mu;
  real<lower=0.5> sens_sd;
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  simplex[K] theta;
  vector[K] beta_wt;
  vector[N] yhat;
  
  theta = dnorm_vec(K, sens_mu, sens_sd);
  beta_wt = theta * beta;
  
  yhat = alpha + X * beta_wt; 
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  sens_mu ~ uniform(0, K);
  sens_sd ~ normal(0.5, K);
 
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
