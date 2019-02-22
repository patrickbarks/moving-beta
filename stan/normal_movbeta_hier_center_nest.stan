
data {
  int N;        // number of data points, length(y)
  int M;        // number of months-within-years
  int J;        // number of years
  vector[N] y;  // response
  matrix[N,M*J] X; // matrix of climate covariates
}

transformed data {
  int m1[M];
  int m2[M];
  int m3[M];
  
  // month-year indices
  for (i in 1:M) {
    m1[i] = i;
    m2[i] = i + M;
    m3[i] = i + 2*M;
  }
}

parameters {
  real alpha;
  real<lower=0> sigma;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[M] beta;
  simplex[J] theta_j;
}

transformed parameters {
  vector[N] yhat;
  vector[M*J] beta_wt;
  
  // apply yearly weights
  beta_wt[m1] = theta_j[1] * beta;
  beta_wt[m2] = theta_j[2] * beta;
  beta_wt[m3] = theta_j[3] * beta;
  
  // linear predictor
  yhat = alpha + X * beta_wt;
}

model {
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  mu_beta ~ normal(0, 3);
  sigma_beta ~ normal(0, 3);
  
  beta ~ normal(mu_beta, sigma_beta);
  theta_j ~ dirichlet(rep_vector(1.0, J));
  
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
