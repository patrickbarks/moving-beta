
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
  real month[M];
  
  // month-year indices
  for (i in 1:M) {
    m1[i] = i;
    m2[i] = i + M;
    m3[i] = i + 2*M;
  }
  
  // vector of months, to create distance matrix
  for(m in 1:M) month[m] = m;
}

parameters {
  real alpha;
  real mu_beta;
  real<lower=0> sigma;
  real<lower=0> eta; // maximum covariance for betas
  real<lower=0> rho; // degree of temporal autocorrelation for betas
  vector[M] z;       // unit normal prior for non-centered term
  simplex[J] theta_j;
}

transformed parameters {
  vector[N] yhat;
  vector[M] beta;
  vector[M*J] beta_wt;
  matrix[M,M] Sigma; // covariance matrix
  matrix[M,M] L;     // cholesky of covariance matrix
  
  // covariance
  Sigma = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.001, M));
  L = cholesky_decompose(Sigma);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z;
  
  // apply yearly weights
  beta_wt[m1] = theta_j[1] * beta;
  beta_wt[m2] = theta_j[2] * beta;
  beta_wt[m3] = theta_j[3] * beta;
  
  // linear predictor
  yhat = alpha + X * beta_wt;
}

model {
  z ~ normal(0, 1);
  
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  alpha ~ normal(0, 5);
  mu_beta ~ normal(0, 5);
  sigma ~ normal(0, 3);
  
  theta_j ~ dirichlet(rep_vector(1.0, J));
  
  y ~ normal(yhat, sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | yhat[n], sigma);
  }
}
