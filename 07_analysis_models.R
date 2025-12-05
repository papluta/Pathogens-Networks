library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
source('01_plotting_functions.R')
source('05_main_data_file.R')


### MODELS
unloadNamespace("bridgesampling") # causes hurdle models to error out
library(brms)
library(tidybayes)
library(tidyr)
library(ggplot2)
library(rstan)
library(bayesplot)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(99)

posterior_correlation <- function(fit = fit, cor = 0.5) {
  post <- posterior_samples(fit) %>% 
    select(starts_with("b_"))
  
  cor_mat <- cor(post)
  
  cor_df <- cor_mat %>%
    as.data.frame() %>%
    mutate(var1 = row.names(.)) %>%
    pivot_longer(
      cols = -var1,
      names_to = "var2",
      values_to = "correlation"
    )
  
  cor_df_clean <- cor_df %>%
    filter(var1 < var2) 
  
  cor_df_clean %>% filter(correlation > cor |
                            correlation < -1*cor)
  
}


#### SITE-LEVEL MODELS ####

# hb.dens.m = brm(log(hb.dens) ~ Density + Year, prior = prior(normal(0,5), class = 'b'), family = 'Gaussian', data.site)
#  
# summary(hb.dens.m)
# pp_check(hb.dens.m, ndraws = 100)
# p_direction(hb.dens.m)
# conditional_effects(hb.dens.m)

con.m <- brm(Connectance ~ sum_nodes + flower_dens.s + total_bee_dens.s + Year + (1|Site), family = Beta(), # INTERACTION NOT SIGNIFICANT
            prior = c(prior(normal(0,2), class = 'b'),
                      prior(exponential(1), class = "sd")), 
            control = list(adapt_delta = 0.95), 
            iter = 5000, warmup = 2000, 
            data = data.site)

summary(con.m)
pp_check(con.m, ndraws = 100)
p_direction(con.m)
pairs(con.m)
posterior_correlation(con.m)


# Beta distribution does not give a good pp_check!
# morisita.hb.m.rich <- brm(Morisita.hb.ns ~ total_bee_dens.s + FL_rich + Year + (1|Species) + (1|Site), family = Beta(),
#                           prior = c(prior(normal(0,2), class = 'b'),
#                                     prior(exponential(1), class = 'sd')), data = data.morisita.hb)
# morisita.bl.m.rich <- brm(Morisita.bl.ns ~ total_bee_dens.s + FL_rich + Year + (1|Species) + (1|Site), family = Beta(), 
#                           prior = c(prior(normal(0,2), class = 'b'),
#                                     prior(exponential(1), class = 'sd')), data = data.morisita.bl)
# 

morisita.hb.m.s <- brm(Morisita.hb.s ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = "Gaussian",
                          prior = c(prior(normal(0,5), class = 'b')), data = data.morisita.hb, iter = 5000, warmup = 2000, 
)
morisita.bl.m.s <- brm(Morisita.bl.s ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = "Gaussian", 
                          prior = c(prior(normal(0,5), class = 'b')), data = data.morisita.bl, iter = 5000, warmup = 2000, 
)


summary(morisita.hb.s)
summary(morisita.bl.s)
pp_check(morisita.hb.m.s, ndraws = 100)
pp_check(morisita.bl.m.s, ndraws = 100)
pairs(morisita.hb.s)
pairs(morisita.bl.s)
p_direction(morisita.hb.m.s)
p_direction(morisita.bl.m.s)
posterior_correlation(morisita.hb.m.s)
posterior_correlation(morisita.bl.m.s)

site_models_500 <- list(con.m = con.m,
                        morisita.bl.m = morisita.bl.m.s, morisita.hb.m = morisita.hb.m.s)

save(site_models_500, file = paste0('Data/Results/',date,'_comm_models.RData'))

#### HURDLE MODELS ####

model_priors <- c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd', dpar = 'mu'))

dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species,
         hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species), 
         family = hurdle_lognormal(), 
         data = data.bb, 
         prior = model_priors, 
         sample_prior = TRUE, 
        iter = 5000, warmup = 2000)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species),
         hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species)), 
         family = hurdle_lognormal(), 
         data = data.wb, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000, 
         control = list(adapt_delta = 0.95))

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species,
         hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species), 
         family = hurdle_lognormal(), 
         data = data.bb, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species),
         hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species)), 
         family = hurdle_lognormal(), 
         data = data.wb, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000, 
         control = list(adapt_delta = 0.95))

abpv.b <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species,
         hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + Species), 
         family = hurdle_lognormal(), 
         data = data.bb.nolp, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000)
abpv.w <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species),
         hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site) + (1|Species)), 
         family = hurdle_lognormal(), 
         data = data.wb, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000, 
         control = list(adapt_delta = 0.95))
abpv.h <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site),
         hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + (1|Site)), 
         family = hurdle_lognormal(), 
         data = data.hb, 
         prior = model_priors, 
         sample_prior = TRUE, 
         iter = 5000, warmup = 2000)

WB_models_combined_raw <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
               bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
               abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h)

save(WB_models_combined_raw, file = paste0('Data/Results/', date, '_models_combined_raw.RData'))

######### SENSITIVITY ANALYSIS ##########
#### HURDLE MODELS with network size ####

dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + Species,
                 hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + Species), 
              family = hurdle_lognormal(), 
              data = data.bb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s  + Year + (1|Site) + (1|Species),
                 hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + (1|Species)), 
              family = hurdle_lognormal(), 
              data = data.wb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000, control = list(adapt_delta = 0.95))

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + Species,
                 hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + Species), 
              family = hurdle_lognormal(), 
              data = data.bb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + (1|Species),
                 hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + (1|Species)), 
              family = hurdle_lognormal(), 
              data = data.wb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000, control = list(adapt_delta = 0.95))

abpv.b <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + Species,
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + sum_nodes +  total_bee_dens.s  + Year + (1|Site) + Species), 
              family = hurdle_lognormal(), 
              data = data.bb.nolp, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000)
abpv.w <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s  + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + (1|Species),
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site) + (1|Species)), 
              family = hurdle_lognormal(), 
              data = data.wb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000, control = list(adapt_delta = 0.95))
abpv.h <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site),
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + sum_nodes +  total_bee_dens.s + Year + (1|Site)), 
              family = hurdle_lognormal(), 
              data = data.hb, 
              prior = model_priors, 
              sample_prior = TRUE, iter = 5000, warmup = 2000)

WB_models_combined_raw_size <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                               bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                               abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h)

save(WB_models_combined_raw_size, file = paste0('Data/Results/',date,'_models_combined_raw_size.RData'))

#### SENSITIVITY ANALYSIS ####
### GP

dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species,
                 hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species),
                 hu ~ dwvb.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE, control = list(adapt_delta = 0.95))

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species,
                 hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species), family = hurdle_lognormal(), data = data.bb, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f.s * Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species),
                 hu ~ bqcv.f.s + Morisita.hb.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE, control = list(adapt_delta = 0.95))

abpv.b <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species,
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + Species), family = hurdle_lognormal(), data = data.bb.nolp, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s  + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species),
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE, control = list(adapt_delta = 0.95))
abpv.h <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F),
                 hu ~ abpv.bl.f.s + Morisita.bl.s + Connectance.s + total_bee_dens.s + Year + gp(x_km, y_km, scale = F)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(student_t(3, 2, 0.5), class = "lscale", coef = "gpx_kmy_km"), 
                  prior(exponential(1), class = "sdgp", coef = "gpx_kmy_km")), sample_prior = TRUE)

WB_models_combined_raw_gp <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                           bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                           abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h)

save(WB_models_combined_raw_gp, file = paste0('Data/Results/',date,'_models_combined_raw_gp.RData'))



posterior_correlation(WB_models_combined_raw[["abpv.h"]])

### PAIRS PLOTS
mcmc_pairs(as.matrix(WB_models_combined_raw[['bqcv.w']]),pars = vars(starts_with("b_") & !contains('hu')))
mcmc_pairs(as.matrix(WB_models_combined_raw[['bqcv.b']]),pars = vars(starts_with("b_hu_")))

pp_check(WB_models_combined_raw[["dwvb"]])
