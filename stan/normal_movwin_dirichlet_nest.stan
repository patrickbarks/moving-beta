
data {
  int N;        // number of data points, length(y)
  int M;        // number of months-within-years
  int J;        // number of years
  vector[N] y;  // response
  matrix[N,M*J] X; // matrix of climate covariates
}

transformed data {
  int K = M*J; // total number of monthly lags
  
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
  simplex[J] theta_y;
  simplex[M] theta_m;
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] yhat;
  vector[K] beta_wt;

  beta_wt[m1] = theta_m * theta_y[1] * beta;
  beta_wt[m2] = theta_m * theta_y[2] * beta;
  beta_wt[m3] = theta_m * theta_y[3] * beta;
  
  yhat = alpha + X * beta_wt;
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
 
  // model
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
