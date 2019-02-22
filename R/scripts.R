
## preliminaries
library(tidyverse)
library(Rcompadre)
library(popbio)
library(rstan)
library(loo)
library(shinystan)
source("R/functions.R")

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# compile stan models
mod_null <- stan_model("stan/normal_null.stan")
mod_year <- stan_model("stan/normal_yr.stan")
mod_mw_gaus <- stan_model("stan/normal_movwin_gaussian.stan")
mod_mw_diri <- stan_model("stan/normal_movwin_dirichlet.stan")
mod_mwn_gaus <- stan_model("stan/normal_movwin_gaussian_nest.stan")
mod_mwn_diri <- stan_model("stan/normal_movwin_dirichlet_nest.stan")
mod_mb_hier <- stan_model("stan/normal_movbeta_hier.stan")
# mod_mb_hier <- stan_model("stan/normal_movbeta_hier_center.stan") # centered
mod_mb_gprc <- stan_model("stan/normal_movbeta_gprc.stan")
mod_mbn_hier <- stan_model("stan/normal_movbeta_hier_nest.stan")
# mod_mbn_hier <- stan_model("stan/normal_movbeta_hier_center_nest.stan") # centered
mod_mbn_gprc <- stan_model("stan/normal_movbeta_gprc_nest.stan")



# load compadre and climate data from Adler
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

clim_tidy <- left_join(ppt_tidy, tmp_tidy, by = c("year", "month", "date"))


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
M <- 12
J <- 3
K <- 36              # total number of monthly lags 
month_start <- "06"  # month 0 for Adler data


## assemble matrix of climate data
X <- matrix(0, N, K)

for(i in seq_along(y)) {
  yr_focal <- year[i]
  date_start <- as.Date(paste(yr_focal-2, month_start, "01", sep = "-"))
  dates_focal <- seq(date_start, by = "month", length.out = K)
  X[i,] <- filter(clim_tidy, date %in% dates_focal)$ppt  # $tmp or $ppt
}



## arrange data for stan
dat <- list(
  N = N, K = K, X = X, y = y,
  M = M, J = J # n_months and n_years for nested models
)

dat_yr1 <- list(N = N, K = K, y = y, x = rowMeans(X[,1:12]))
dat_yr2 <- list(N = N, K = K, y = y, x = rowMeans(X[,13:24]))
dat_yr3 <- list(N = N, K = K, y = y, x = rowMeans(X[,25:36]))


## fit stan models
fit_control <- stanfn(mod_null)
fit_cyear_1 <- stanfn(mod_year, data = dat_yr1)
fit_cyear_2 <- stanfn(mod_year, data = dat_yr2)
fit_cyear_3 <- stanfn(mod_year, data = dat_yr3)
fit_mw_gaus <- stanfn(mod_mw_gaus, control = ctrl1)
fit_mw_diri <- stanfn(mod_mw_diri)
fit_mwn_gaus <- stanfn(mod_mwn_gaus)
fit_mwn_diri <- stanfn(mod_mwn_diri)
fit_mb_hier <- stanfn(mod_mb_hier, control = ctrl2)
fit_mb_gpr2 <- stanfn(mod_mb_gprc, control = ctrl2)

# new nested moving-beta models
fit_mbn_hier <- stanfn(mod_mbn_hier, control = ctrl2)
fit_mbn_gpr2 <- stanfn(mod_mbn_gprc, control = ctrl2)


## model type labels
mod_lev <- c("control (null)",
             "year t",
             "year t-1",
             "year t-2",
             "moving-window (gaussian)",
             "moving-window (dirichlet)",
             "nested moving-window (gaussian)",
             "nested moving-window (dirichlet)",
             "moving-beta (hierarchical)",
             "moving-beta (gaus proc 2)",
             "nested moving-beta (hierarchical)",
             "nested moving-beta (gaus proc 2)")


## plot lagged betas by model type
df_betas <- rbind(
  summarize_beta(fit_mw_gaus, "moving-window (gaussian)", wt = TRUE),
  summarize_beta(fit_mw_diri, "moving-window (dirichlet)", wt = TRUE),
  summarize_beta(fit_mwn_gaus, "nested moving-window (gaussian)", wt = TRUE),
  summarize_beta(fit_mwn_diri, "nested moving-window (dirichlet)", wt = TRUE),
  summarize_beta(fit_mb_hier, "moving-beta (hierarchical)"),
  summarize_beta(fit_mb_gpr2, "moving-beta (gaus proc 2)"),
  summarize_beta(fit_mbn_hier, "nested moving-beta (hierarchical)", wt = TRUE),
  summarize_beta(fit_mbn_gpr2, "nested moving-beta (gaus proc 2)", wt = TRUE)
) %>% mutate(model = factor(model, levels = mod_lev))


ggplot(df_betas) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_point(aes(lag, beta_med)) +
  geom_linerange(aes(lag, ymin = beta_low90, ymax = beta_upp90), size = 0.7) +
  facet_wrap(~ model) +
  labs(x = "Months before present", y = "beta") +
  theme(panel.grid = element_blank())


## plot observed vs. predicted values
df_yhat <- rbind(
  summarize_yhat(fit_control, "control (null)"),
  summarize_yhat(fit_cyear_1, "year t"),
  summarize_yhat(fit_cyear_2, "year t-1"),
  summarize_yhat(fit_cyear_3, "year t-2"),
  summarize_yhat(fit_mw_gaus, "moving-window (gaussian)"),
  summarize_yhat(fit_mw_diri, "moving-window (dirichlet)"),
  summarize_yhat(fit_mwn_gaus, "nested moving-window (gaussian)"),
  summarize_yhat(fit_mwn_diri, "nested moving-window (dirichlet)"),
  summarize_yhat(fit_mb_hier, "moving-beta (hierarchical)"),
  summarize_yhat(fit_mb_gpr2, "moving-beta (gaus proc 2)"),
  summarize_yhat(fit_mbn_hier, "nested moving-beta (hierarchical)"),
  summarize_yhat(fit_mbn_gpr2, "nested moving-beta (gaus proc 2)")
) %>% mutate(model = factor(model, levels = mod_lev))

ggplot(df_yhat) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_point(aes(y, yhat_med)) +
  geom_smooth(aes(y, yhat_med), method = "lm", se = FALSE) +
  geom_linerange(aes(x = y, ymin = yhat_low90, ymax = yhat_upp90)) +
  facet_wrap(~ model) +
  labs(x = "Observed", y = "Predicted")



### measures of fit
rbind(
  summarize_fit(fit_control, "control (null)"),
  summarize_fit(fit_cyear_1, "year t"),
  summarize_fit(fit_cyear_2, "year t-1"),
  summarize_fit(fit_cyear_3, "year t-2"),
  summarize_fit(fit_mw_gaus, "moving-window (gaussian)"),
  summarize_fit(fit_mw_diri, "moving-window (dirichlet)"),
  summarize_fit(fit_mwn_gaus, "nested moving-window (gaussian)"),
  summarize_fit(fit_mwn_diri, "nested moving-window (dirichlet)"),
  summarize_fit(fit_mb_hier, "moving-beta (hierarchical)"),
  summarize_fit(fit_mb_gpr2, "moving-beta (gaus proc 2)"),
  summarize_fit(fit_mbn_hier, "nested moving-beta (hierarchical)"),
  summarize_fit(fit_mbn_gpr2, "nested moving-beta (gaus proc 2)")
)




# ## check equivalence in weights * x * beta
# N <- 4
# K <- 3
# X <- matrix(rnorm(N*K), nrow = N, ncol = K)
# beta <- 0.3
# 
# weights <- runif(K)
# weights <- weights / sum(weights)
# 
# x_ant <- numeric(N)
# x_ant[1] <- sum(X[1,] * weights)
# x_ant[2] <- sum(X[2,] * weights)
# x_ant[3] <- sum(X[3,] * weights)
# x_ant[4] <- sum(X[4,] * weights)
# 
# sum(x_ant %*% beta)
# 
# beta_weighted <- weights * beta
# sum(X %*% beta_weighted)

