
## preliminaries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
library(rstan)
library(loo)
library(shinystan)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# load compadre and data from Adler
comp <- cdb_fetch("~/COMPADRE_v.X.X.X.RData")
ppt <- read_csv("data/monthly_ppt.csv")
tmp <- read_csv("data/monthly_temp.csv")


## de-trend and tidy climate data
ppt_tidy <- ppt %>%
  filter(YEAR >= 1930) %>%
  mutate_at(vars(JAN:DEC), function(x) as.numeric(scale(x))) %>% 
  setNames(c("year", 1:12)) %>% 
  gather(month, ppt, -year) %>% 
  mutate(month = as.integer(month)) %>% 
  mutate(ppt = ifelse(is.na(ppt), 0, ppt)) %>% 
  arrange(year, month) %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

tmp_tidy <- tmp %>%
  filter(YEAR >= 1930) %>%
  mutate_at(vars(JAN:DEC), function(x) as.numeric(scale(x))) %>% 
  setNames(c("year", 1:12)) %>% 
  gather(month, tmp, -year) %>% 
  mutate(month = as.integer(month)) %>% 
  mutate(tmp = ifelse(is.na(tmp), 0, tmp)) %>% 
  arrange(year, month) %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

clim_tidy <- left_join(ppt_tidy, tmp_tidy)


## select focal species
spp <- "Cirsium_undulatum"
# Sphaeralcea_coccinea, Psoralea_tenuiflora, Paronychia_jamesii,
# Echinacea_angustifolia, Solidago_mollis, Cirsium_undulatum,
# Thelesperma_megapotamicum, Ratibida_columnifera, Hedyotis_nigricans,
# Lesquerella_ovalifolia


## filter compadre to individual annual mpms from species of interest
comp_spp <- comp %>% 
  filter(MatrixComposite == "Individual",
         MatrixTreatment == "Unmanipulated",
         SpeciesAuthor %in% spp) %>% 
  cdb_flag(c("check_ergodic", "check_zero_U")) %>%
  filter(check_zero_U == FALSE) %>% # one matrix all 0s
  cdb_unnest() %>% 
  mutate(lambda = map_dbl(matA(.), popbio::lambda)) %>% 
  mutate(log_lambda = log(lambda)) %>% 
  mutate(year = as.Date(paste0(MatrixStartYear, "-01-01")))


## prepare data for stan
response <- "log_lambda"
year <- comp_spp$MatrixEndYear
y <- comp_spp[[response]]
N <- length(y)
K <- 36              # total number of monthly lags 
month_start <- "06"  # month 1 for Adler data


## assemble matrix of climate data
X <- matrix(0, N, K)

for(i in seq_along(y)) {
  yr_focal <- year[i]
  date_start <- as.Date(paste(yr_focal-2, month_start, "01", sep = "-"))
  dates_focal <- seq(date_start, by = "month", length.out = K)
  X[i,] <- filter(clim_tidy, date %in% dates_focal)$tmp  # $tmp or $ppt
}


## arrange data for stan
dat_stan <- list(
  N = N, K = K, X = X, y = y,
  x = rowMeans(X[,1:12]), # first-yr climate average for ctrl_yr model
  M = 12, J = 3 # n_years and n_years for nested models
)


## convenience functions for fitting stan models
stanfn <- function(file, control = list()) {
  stan(file = paste0("stan/", file), data = dat_stan, warmup = 1500,
       iter = 3500, thin = 2, chains = 2, control = control,
       pars = c("z", "Sigma", "L"), include = FALSE)
}

ctrl1 <- list(adapt_delta = 0.95, stepsize = 0.05)
ctrl2 <- list(adapt_delta = 0.99, stepsize = 0.01)
ctrl3 <- list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 12)


## fit stan models
fit_control <- stanfn("normal_null.stan")
fit_cyear_1 <- stanfn("normal_yr.stan")
fit_mb_hier <- stanfn("normal_movbeta_hier.stan", control = ctrl2)
# fit_mb_hier <- stanfn("normal_movbeta_hier_center.stan", control = ctrl2)
fit_mb_gpr1 <- stanfn("normal_movbeta_gpr1.stan")
fit_mb_gpr2 <- stanfn("normal_movbeta_gpr2.stan", control = ctrl2)
fit_mw_gaus <- stanfn("normal_movwin_gaussian_nest.stan")
fit_mw_diri <- stanfn("normal_movwin_dirichlet_nest.stan")


## model type labels
mod_lev <- c("control (null)",
             "control (average in year t)",
             "moving-window (gaussian, nested)",
             "moving-window (dirichlet, nested)",
             "moving-beta (hierarchical)",
             "moving-beta (gaus proc 1)",
             "moving-beta (gaus proc 2)")


## extract posterior samples
beta_mb_hier <- rstan::extract(fit_mb_hier, 'beta')$beta
beta_mb_gpr1 <- rstan::extract(fit_mb_gpr1, 'beta')$beta
beta_mb_gpr2 <- rstan::extract(fit_mb_gpr2, 'beta')$beta
beta_mw_gaus <- rstan::extract(fit_mw_gaus, 'beta_wt')$beta_wt
beta_mw_diri <- rstan::extract(fit_mw_diri, 'beta_wt')$beta_wt

yhat_control <- rstan::extract(fit_control, 'yhat')$yhat
yhat_cyear_1 <- rstan::extract(fit_cyear_1, 'yhat')$yhat
yhat_mb_hier <- rstan::extract(fit_mb_hier, 'yhat')$yhat
yhat_mb_gpr1 <- rstan::extract(fit_mb_gpr1, 'yhat')$yhat
yhat_mb_gpr2 <- rstan::extract(fit_mb_gpr2, 'yhat')$yhat
yhat_mw_gaus <- rstan::extract(fit_mw_gaus, 'yhat')$yhat
yhat_mw_diri <- rstan::extract(fit_mw_diri, 'yhat')$yhat

ll_control <- extract_log_lik(fit_control)
ll_cyear_1 <- extract_log_lik(fit_cyear_1)
ll_mb_hier <- extract_log_lik(fit_mb_hier)
ll_mb_gpr1 <- extract_log_lik(fit_mb_gpr1)
ll_mb_gpr2 <- extract_log_lik(fit_mb_gpr2)
ll_mw_gaus <- extract_log_lik(fit_mw_gaus)
ll_mw_diri <- extract_log_lik(fit_mw_diri)


## plot lagged betas by model type
extract_beta <- function(post_samples, model) {
  tibble(
    month = seq_len(ncol(post_samples)),
    model = model,
    beta_med = apply(post_samples, 2, function(x) quantile(x, 0.500)),
    beta_low90 = apply(post_samples, 2, function(x) quantile(x, 0.05)),
    beta_upp90 = apply(post_samples, 2, function(x) quantile(x, 0.95))
  )
}

df_betas <- rbind(
  extract_beta(beta_mb_hier, "moving-beta (hierarchical)"),
  extract_beta(beta_mb_gpr1, "moving-beta (gaus proc 1)"),
  extract_beta(beta_mb_gpr2, "moving-beta (gaus proc 2)"),
  extract_beta(beta_mw_gaus, "moving-window (gaussian, nested)"),
  extract_beta(beta_mw_diri, "moving-window (dirichlet, nested)")
) %>% mutate(model = factor(model, levels = mod_lev))


ggplot(df_betas) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_point(aes(month, beta_med)) +
  geom_linerange(aes(month, ymin = beta_low90, ymax = beta_upp90), size = 0.7) +
  facet_wrap(~ model) +
  labs(x = "Months before present", y = "beta") +
  theme(panel.grid = element_blank())


## plot observed vs. predicted values
extract_yhat <- function(post_samples, model) {
  tibble(
    n = ncol(post_samples),
    y = y,
    model = model,
    yhat_med = apply(post_samples, 2, function(x) quantile(x, 0.500)),
    yhat_low90 = apply(post_samples, 2, function(x) quantile(x, 0.05)),
    yhat_upp90 = apply(post_samples, 2, function(x) quantile(x, 0.95))
  )
}

df_yhat <- rbind(
  extract_yhat(yhat_control, "control (null)"),
  extract_yhat(yhat_cyear_1, "control (average in year t)"),
  extract_yhat(yhat_mb_hier, "moving-beta (hierarchical)"),
  extract_yhat(yhat_mb_gpr1, "moving-beta (gaus proc 1)"),
  extract_yhat(yhat_mb_gpr2, "moving-beta (gaus proc 2)"),
  extract_yhat(yhat_mw_gaus, "moving-window (gaussian, nested)"),
  extract_yhat(yhat_mw_diri, "moving-window (dirichlet, nested)")
) %>% mutate(model = factor(model, levels = mod_lev))

ggplot(df_yhat) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_point(aes(y, yhat_med)) +
  geom_smooth(aes(y, yhat_med), method = "lm", se = FALSE) +
  geom_linerange(aes(x = y, ymin = yhat_low90, ymax = yhat_upp90)) +
  facet_wrap(~ model) +
  labs(x = "Observed", y = "Predicted")



### within-model measures of fit

# lppd
sum(log(colMeans(exp(ll_control))))
sum(log(colMeans(exp(ll_cyear_1))))
sum(log(colMeans(exp(ll_mb_hier))))
sum(log(colMeans(exp(ll_mb_gpr1))))
sum(log(colMeans(exp(ll_mb_gpr2))))
sum(log(colMeans(exp(ll_mw_gaus))))
sum(log(colMeans(exp(ll_mw_diri))))

# loo estimates
loo(fit_control)$estimates
loo(fit_cyear_1)$estimates
loo(fit_mb_hier)$estimates
loo(fit_mb_gpr1)$estimates
loo(fit_mb_gpr2)$estimates
loo(fit_mw_gaus)$estimates
loo(fit_mw_diri)$estimates

# waic 
waic(extract_log_lik(fit_control))
waic(extract_log_lik(fit_cyear_1))
waic(extract_log_lik(fit_mb_hier))
waic(extract_log_lik(fit_mb_gpr1))
waic(extract_log_lik(fit_mb_gpr2))
waic(extract_log_lik(fit_mw_gaus))
waic(extract_log_lik(fit_mw_diri))

