library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
source('00_plotting_functions.R')
source('05_main_data_file.R')

### MODELS

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

#### SITE-LEVEL MODELS ####

ab.m <- brm(log(sum.ab) ~ FL.sum + Year + (1|Site), family = 'Gaussian',  # INTERACTION NOT SIGNIFICANT
            prior = prior(normal(0,5), class = 'b'), data = data.site)
ab.m <- brm(log(sum.ab) ~ FL_rich + Year + (1|Site), family = 'Gaussian',  # INTERACTION NOT SIGNIFICANT
            prior = prior(normal(0,5), class = 'b'), data = data.site)

summary(ab.m)
pp_check(ab.m, ndraws = 100)
p_direction(ab.m)
conditional_effects(ab.m)

con.m <- brm(Connectance ~ FL_rich * log(sum.ab.d) + Year + (1|Site), family = Beta(), # INTERACTION NOT SIGNIFICANT
            prior = prior(normal(0,5), class = 'b'), data = data.site)
summary(con.m)
pp_check(con.m, ndraws = 100)
p_direction(con.m)
pairs(con.m)


morisita.hb.m <- brm(Morisita.z ~ sum.ab * FL_rich + Year + (1|Species) + (1|Site), family = 'Gaussian',
                  prior = prior(normal(0,5), class = 'b'), data = data.morisita.hb)
morisita.bl.m <- brm(Morisita.bl.z ~ sum.ab * FL_rich + Year + (1|Species) + (1|Site), family = 'Gaussian', # INTERACTION NOT SIGNIFICANT
                  prior = prior(normal(0,5), class = 'b'), data = data.morisita.bl)
summary(morisita.hb.m)
pp_check(morisita.hb.m, ndraws = 100)
p_direction(morisita.hb.m)
pairs(morisita.bl.m)
conditional_effects(morisita.bl.m)

site_models_500 <- list(ab.m = ab.m, con.m = con.m,
                        morisita.bl.m = morisita.bl.m, morisita.hb.m = morisita.hb.m)


#save(site_models_500, file = '250130_comm_models.RData')


#### HURDLE MODELS ####

dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.z.s + Year + (1 |Site) + Species,
                 hu ~ dwvb.f.s + Morisita.z.s + Year + (1 |Site) + Species), family = hurdle_lognormal(), data = data.bb.forhb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f.s + Morisita.z.s + Year + (1 |Site) + (1|Species),
                 hu ~ dwvb.f.s + Morisita.z.s + Year + (1 |Site)+ (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE, control = list(adapt_delta = 0.95))

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f.s + Morisita.z.s + Year + (1 |Site) + Species,
                 hu ~ bqcv.f.s + Morisita.z.s + Year + (1 |Site)+ Species), family = hurdle_lognormal(), data = data.bb.nolp, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f.s + Morisita.z.s + Year + (1 |Site) + (1|Species),
                 hu ~ bqcv.f.s + Morisita.z.s + Year + (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE, control = list(adapt_delta = 0.95))

abpv.b <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.z.s + Year + (1 |Site)+ Species,
                 hu ~ abpv.bl.f.s + Morisita.bl.z.s + Year + (1 |Site)+ Species), family = hurdle_lognormal(), data = data.bb.nolp, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.z.s  + Year + (1 |Site) + (1|Species),
                 hu ~ abpv.bl.f.s + Morisita.bl.z.s + Year + (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE, control = list(adapt_delta = 0.95))
abpv.h <- brm(bf(ABPV.abs ~ abpv.bl.f.s + Morisita.bl.z.s + Year + (1 |Site),
                 hu ~ abpv.bl.f.s + Morisita.bl.z.s + Year + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)


WB_models_Morisita <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                           bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                           abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h)


###### connectance

dwvb.b <- brm(bf(DWVB.abs ~ sum.ab.fl + Connectance + Year + (1 |Site) + Species,
                 hu ~ sum.ab.fl + Connectance + Year + (1 |Site) + Species), family = hurdle_lognormal(), data = data.bb.forhb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site) + (1|Species),
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site)+ (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
dwvb.h <- brm(bf(DWVB.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site),
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)


bqcv.b <- brm(bf(BQCV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site) + Species,
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site)+ Species), family = hurdle_lognormal(), data = data.bb.forhb, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site) + (1|Species),
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior =
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
bqcv.h <- brm(bf(BQCV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site)), family = lognormal(), data = data.hb, prior =
                c(prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)

abpv.b <- brm(bf(ABPV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site) + Species,
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site) + Species), family = hurdle_lognormal(), data = data.bb.nolp, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site) + (1|Species),
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)
abpv.h <- brm(bf(ABPV.abs ~ sum.ab.fl +  Connectance + Year + (1 |Site),
                 hu ~ sum.ab.fl +  Connectance + Year + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,2), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b'),
                  prior(exponential(1), class = 'sd', dpar = 'hu'),
                  prior(exponential(1), class = 'sd')), sample_prior = TRUE)

WB_models_Connectance <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, dwvb.h = dwvb.h, 
                           bqcv.b = bqcv.b, bqcv.w = bqcv.w, bqcv.h = bqcv.h,
                           abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h, abpv.bl = abpv.bl)


#save(WB_models_Connectance, WB_models_Morisita,WB_models_Connectance_sp, file = '250304_models_mor_con_sep.RData')


### PAIRS PLOTS
mcmc_pairs(as.matrix(WB_models_Morisita[['sbv.b']]),pars = vars(starts_with("b_") & !contains('hu')))
mcmc_pairs(as.matrix(WB_models_Morisita[['sbv.b']]),pars = vars(starts_with("b_hu_")))

### MODEL COMPARISON
loo.list = lapply(WB_models_combined2_x, waic)
loo.list.con = lapply(WB_models_Connectance_sp2, waic)

est.list = lapply(loo.list, function(x) x$estimates[3,]) %>% bind_rows(.id = 'mod') %>% rename(Est.mor = Estimate, SE.mor = SE)
est.list.con = lapply(loo.list.con, function(x) x$estimates[3,]) %>% bind_rows() %>% rename(Est.int = Estimate, SE.int = SE)
waic.sum = cbind(est.list, est.list.con) %>% mutate(mor_con = Est.mor - Est.int)