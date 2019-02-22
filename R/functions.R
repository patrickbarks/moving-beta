

ctrl1 <- list(adapt_delta = 0.95, stepsize = 0.05)
ctrl2 <- list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 11)
ctrl3 <- list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 12)


stan_diagnostics <- function(fit) {
  ll <- extract_log_lik(fit)
  stan_summary <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(stan_summary$Rhat > 1.1))
  n_eff_low <- length(which(stan_summary$n_eff / nrow(ll) < 0.1))
  mcse_high <- length(which(stan_summary$se_mean / stan_summary$sd > 0.1))
  n_diverg <- rstan::get_num_divergent(fit)
  return(tibble(rhat_high, n_eff_low, mcse_high, n_diverg))
}


stanfn <- function(object, control = ctrl1, data = dat) {
  fit <- sampling(object = object, data = data, warmup = 1500,
                  iter = 3500, thin = 2, chains = 2, control = control,
                  pars = c("z", "Sigma", "L"), include = FALSE)
  
  # if signs of poor convergence, re-fit with ctrl2
  if (any(stan_diagnostics(fit) > 0)) {
    fit <- sampling(object = object, data = data, warmup = 1500,
                    iter = 3500, thin = 2, chains = 2, control = ctrl2,
                    pars = c("z", "Sigma", "L"), include = FALSE)
    
    # if signs of poor convergence, re-fit with ctrl3
    if (any(stan_diagnostics(fit) > 0)) {
      fit <- sampling(object = object, data = data, warmup = 1500,
                      iter = 3500, thin = 2, chains = 2, control = ctrl3,
                      pars = c("z", "Sigma", "L"), include = FALSE)
    }
  }
  
  return(fit)
}


summarize_yhat <- function(fit, label) {
  yhat <- rstan::extract(fit, 'yhat')$yhat
  
  tibble(n = seq_len(ncol(yhat)),
         y = y,
         model = label,
         yhat_med = apply(yhat, 2, function(x) quantile(x, 0.50)),
         yhat_low90 = apply(yhat, 2, function(x) quantile(x, 0.05)),
         yhat_upp90 = apply(yhat, 2, function(x) quantile(x, 0.95)))
}


summarize_beta <- function(fit, label, wt = FALSE) {
  if (wt == TRUE) {
    beta <- rstan::extract(fit, 'beta_wt')$beta_wt
  } else {
    beta <- rstan::extract(fit, 'beta')$beta
  }
  
  tibble(lag = seq_len(ncol(beta)),
         model = label,
         beta_med = apply(beta, 2, function(x) quantile(x, 0.50)),
         beta_low90 = apply(beta, 2, function(x) quantile(x, 0.05)),
         beta_upp90 = apply(beta, 2, function(x) quantile(x, 0.95)))
}


summarize_fit <- function(fit, label) {
  ll <- extract_log_lik(fit)
  lppd <- sum(log(colMeans(exp(ll))))
  
  loo_mat <- suppressWarnings(loo(fit)$estimates)
  elpd_loo <- loo_mat[1,1]
  elpd_loo_se <- loo_mat[1,2]
  
  waic_mat <- suppressWarnings(waic(ll)$estimates)
  elpd_waic <- waic_mat[1,1]
  elpd_waic_se <- waic_mat[1,2]
  
  cbind(tibble(model = label),
        stan_diagnostics(fit),
        tibble(elpd_loo, elpd_loo_se, elpd_waic, elpd_waic_se))
}


summarize_xval <- function(fit, label) {
  ll_test <- extract_log_lik(fit, "log_lik_test")
  lppd_test <- sum(log(colMeans(exp(ll_test))))
  
  yhat_test <- rstan::extract(fit, "yhat_test")$yhat_test
  yhat_test_median <- apply(yhat_test, 2, median)
  
  bind_cols(tibble(model = label),
            stan_diagnostics(fit),
            tibble(lppd_test, yhat_test = yhat_test_median))
}

