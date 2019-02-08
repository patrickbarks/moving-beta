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
  int N;        // number of data points, length(y)
  int M;        // number of months-within-years
  int J;        // number of years
  vector[N] y;  // response
  matrix[N,M*J] X; // matrix of climate covariates
}

transformed data {
  int K = M*J; // total number months
  
  int m1[M];   // month-year indices
  int m2[M];
  int m3[M];
  
  for (i in 1:M) {
    m1[i] = i;
    m2[i] = i + M;
    m3[i] = i + 2*M;
  }
}

parameters {
  real<lower=0, upper=M> sens_mu;
  real<lower=0.5> sens_sd;
  simplex[J] theta_y;
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] yhat;
  simplex[M] theta_m;
  vector[K] beta_wt;
  
  theta_m = dnorm_vec(M, sens_mu, sens_sd);
  
  beta_wt[m1] = theta_m * theta_y[1] * beta;
  beta_wt[m2] = theta_m * theta_y[2] * beta;
  beta_wt[m3] = theta_m * theta_y[3] * beta;
  
  yhat = alpha + X * beta_wt; 
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  sens_mu ~ uniform(0, M);
  sens_sd ~ normal(0.5, M);
 
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
