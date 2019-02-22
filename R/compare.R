
## preliminaries
library(tidyverse)
library(Rcompadre)
library(popbio)
library(rstan)
library(loo)
library(shinystan)
library(gridExtra)
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
mod_mb_gprc <- stan_model("stan/normal_movbeta_gprc.stan")
mod_mbn_hier <- stan_model("stan/normal_movbeta_hier_nest.stan")
mod_mbn_gprc <- stan_model("stan/normal_movbeta_gprc_nest.stan")


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

clim_tidy <- left_join(ppt_tidy, tmp_tidy, by = c("year", "month", "date"))


## select focal species
spp <- c("Sphaeralcea_coccinea", "Psoralea_tenuiflora", "Paronychia_jamesii",
         "Echinacea_angustifolia", "Solidago_mollis", "Cirsium_undulatum",
         "Thelesperma_megapotamicum", "Ratibida_columnifera",
         "Hedyotis_nigricans", "Lesquerella_ovalifolia")

## select climate variable
clim_var <- "tmp"

## filter to Adler species, individual unmanipulated matrices
comp_full <- comp %>% 
  filter(MatrixComposite == "Individual",
         MatrixTreatment == "Unmanipulated",
         SpeciesAuthor %in% spp) %>% 
  cdb_flag(c("check_ergodic", "check_zero_U")) %>%
  filter(check_zero_U == FALSE) %>% # one matrix all 0s
  cdb_unnest() %>% 
  mutate(lambda = map_dbl(matA(.), popbio::lambda)) %>% 
  mutate(log_lambda = log(lambda)) %>% 
  mutate(year = as.Date(paste0(MatrixStartYear, "-01-01")))


## constant params
M <- 12
J <- 3
K <- 36              # total number of monthly lags 
month_start <- "06"  # month 1 for Adler data


for (i in 1:length(spp)) {
  
  spp_focal <- spp[i]
  
  ## filter compadre to individual annual mpms from species of interest
  comp_spp <- filter(comp_full, SpeciesAuthor == spp_focal)
  
  ## prepare data for stan
  year <- comp_spp$MatrixEndYear
  y <- comp_spp$log_lambda
  N <- length(y)
  
  ## assemble matrix of climate data
  X <- matrix(0, N, K)
  
  for(i in seq_along(y)) {
    yr_focal <- year[i]
    date_start <- as.Date(paste(yr_focal-2, month_start, "01", sep = "-"))
    dates_focal <- seq(date_start, by = "month", length.out = K)
    X[i,] <- filter(clim_tidy, date %in% dates_focal)[[clim_var]]
  }
  
  ## arrange data for stan
  dat <- list(
    N = N, K = K, X = X, y = y,
    M = M, J = J # n_months and n_years for nested models
  )
  
  dat_yr1 <- list(N = N, y = y, x = rowMeans(X[,1:12]))
  dat_yr2 <- list(N = N, y = y, x = rowMeans(X[,13:24]))
  dat_yr3 <- list(N = N, y = y, x = rowMeans(X[,25:36]))
  
  ## fit stan models
  mod <- tibble(
    fit = vector(mode = "list", length = 12),
    model = c("null", "year t", "year t-1", "year t-2",
              "move-win (gaussian)", "move-win (dirichlet)",
              "move-win-nest (gaussian)", "move-win-nest (dirichlet)",
              "move-beta (hier)", "move-beta (gprc)", "move-beta-nest (hier)",
              "move-beta-nest (gprc)"),
    wt = ifelse(grepl("nest|win", model), TRUE, FALSE)
  )
  
  mod$fit[[1]] <- stanfn(mod_null)
  mod$fit[[2]] <- stanfn(mod_year, data = dat_yr1)
  mod$fit[[3]] <- stanfn(mod_year, data = dat_yr2)
  mod$fit[[4]] <- stanfn(mod_year, data = dat_yr3)
  mod$fit[[5]] <- stanfn(mod_mw_gaus)
  mod$fit[[6]] <- stanfn(mod_mw_diri)
  mod$fit[[7]] <- stanfn(mod_mwn_gaus)
  mod$fit[[8]] <- stanfn(mod_mwn_diri)
  mod$fit[[9]] <- stanfn(mod_mb_hier, control = ctrl1)
  mod$fit[[10]] <- stanfn(mod_mb_gprc, control = ctrl2)
  mod$fit[[11]] <- stanfn(mod_mbn_hier, control = ctrl2)
  mod$fit[[12]] <- stanfn(mod_mbn_gprc, control = ctrl2)
  
  # write stanfit model objects to file
  file_mod <- paste0("stanfit/fit_", clim_var, "_", spp_focal, ".RData")
  save(mod, file = file_mod)
  
  ## plot observed vs. predicted
  df_yhat <- mod %>% 
    rowwise() %>% 
    do(summarize_yhat(.$fit, .$model)) %>% 
    ungroup() %>% 
    mutate(model = factor(model, levels = mod$model))
  
  p_yhat <- ggplot(df_yhat) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    geom_point(aes(y, yhat_med), alpha = 0.8) +
    geom_linerange(aes(x = y, ymin = yhat_low90, ymax = yhat_upp90), alpha = 0.8) +
    geom_smooth(aes(y, yhat_med), method = "lm", se = FALSE) +
    facet_wrap(~ model) +
    labs(x = "Observed", y = "Predicted") +
    theme(panel.grid = element_blank())
  
  file_yhat <- paste0("img/yhat_", clim_var, "_", spp_focal, ".png")
  ggsave(file_yhat, p_yhat, width = 7.5, height = 6, units = "in", dpi = 150)
  
  ## plot lagged betas by model type
  df_betas <- mod %>% 
    filter(grepl("move", model)) %>% 
    rowwise() %>% 
    do(summarize_beta(.$fit, .$model, .$wt)) %>% 
    ungroup() %>% 
    mutate(model = factor(model, levels = mod$model))
  
  p_betas <- ggplot(df_betas) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_point(aes(lag, beta_med)) +
    geom_linerange(aes(lag, ymin = beta_low90, ymax = beta_upp90), size = 0.7) +
    facet_wrap(~ model) +
    labs(x = "Months before present", y = "beta") +
    theme(panel.grid = element_blank())
  
  file_betas <- paste0("img/betas_", clim_var, "_", spp_focal, ".png")
  ggsave(file_betas, p_betas, width = 7, height = 6, units = "in", dpi = 150)
  
  ## summarize model fit and diagnostics
  df_fit <- mod %>% 
    rowwise() %>% 
    do(summarize_fit(.$fit, .$model)) %>% 
    ungroup() %>% 
    mutate(model = factor(model, levels = mod$model))
  
  file_fit <- paste0("analysis/fit_", clim_var, "_", spp_focal, ".csv")
  write.csv(df_fit, file = file_fit, row.names = FALSE)
}





## analysis
df_mod <- tibble(
  model = c("null", "year t", "year t-1", "year t-2",
            "move-win (gaussian)", "move-win (dirichlet)",
            "move-win-nest (gaussian)", "move-win-nest (dirichlet)",
            "move-beta (hier)", "move-beta (gprc)",
            "move-beta-nest (hier)", "move-beta-nest (gprc)"),
  mod_lab = c("null", "year t", "year t-1", "year t-2",
              "mw (gaus)", "mw (dirich)",
              "mw-n (gaus)", "mw-n (dirich)",
              "m-beta (hier)", "m-beta (gprc)",
              "m-beta-nest (hier)", "m-beta-nest (gprc)"))

readfn <- function(file) {
  df <- read_csv(file)
  df$species <- paste(strsplit(file, "_|\\.")[[1]][3:4], collapse = "_")
  return(df)
}


# get file names
fit_files <- paste0("analysis/", list.files("analysis"))
fit_files_ppt <- fit_files[grepl("ppt", fit_files)]
fit_files_tmp <- fit_files[grepl("tmp", fit_files)]

# read and bind files
df_ppt <- bind_rows(lapply(fit_files_ppt, readfn))
df_tmp <- bind_rows(lapply(fit_files_tmp, readfn))

# check convergence issues
filter(df_ppt, rhat_high > 0 | n_eff_low > 0)
filter(df_tmp, rhat_high > 0 | n_eff_low > 0)

# scale waic relative to null, and arrange for plotting
arrange_waic <- function(df) {
  df %>% 
    group_by(species) %>% 
    mutate(elpd_loo = elpd_loo - elpd_loo[which(model == "null")],
           elpd_waic = elpd_waic - elpd_waic[which(model == "null")]) %>% 
    ungroup() %>% 
    mutate(model = factor(model,
                          levels = df_mod$model,
                          labels = df_mod$mod_lab)) %>% 
    filter(model != "null")
}

waic_ppt <- arrange_waic(df_ppt)
waic_tmp <- arrange_waic(df_tmp)

# which models have waic >1 SE above null
sig_ppt <- filter(waic_ppt, elpd_waic - elpd_waic_se > 0)
sig_tmp <- filter(waic_tmp, elpd_waic - elpd_waic_se > 0)


# plot heatmap of WAIC values
p_ppt <- ggplot(waic_ppt, aes(model, species)) +
  geom_tile(aes(fill = elpd_waic)) +
  geom_point(data = sig_ppt, col = "yellow", shape = 4, size = 1.2, stroke = 1) +
  scale_fill_gradient2(name = expression(paste(Delta, elpd[waic])), guide = F, limits = c(-3, 13)) +
  labs(x = NULL, y = NULL) +
  ggtitle(expression(paste("log-", lambda, " ~ precipitation"))) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(vjust = -0.5),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

p_tmp <- ggplot(waic_tmp, aes(model, species)) +
  geom_tile(aes(fill = elpd_waic)) +
  geom_point(data = sig_tmp, col = "yellow", shape = 4, size = 1.2, stroke = 1) +
  scale_fill_gradient2(name = expression(paste(Delta, elpd[waic])), limits = c(-3, 13)) +
  labs(x = NULL, y = NULL) +
  ggtitle(expression(paste("log-", lambda, " ~ temperature"))) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -0.5),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

g <- cbind(ggplotGrob(p_ppt), ggplotGrob(p_tmp), size = "last")

dev.off()
quartz(height = 5, width = 8, dpi = 150)
grid.arrange(g)

# ggsave("img/elpd_waic.png", g, height = 5, width = 8, units = "in", dpi = 200)

