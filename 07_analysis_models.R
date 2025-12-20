source('01_custom_functions.R')
source('05_main_data_file.R') # reruns also previous scripts (01 - 05)


### MODELS
unloadNamespace("brms") # unloading brms to unload bridgesampling
unloadNamespace("bridgesampling") # unloading because it causes hurdle models to error out
library(brms)
library(tidybayes)
library(rstan)
library(bayesplot)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(99)


#### SITE-LEVEL MODELS ####

# for Fig S1
hb.dens.m = brm(log(hb_dens) ~ Density + Year, prior = prior(normal(0,5), class = 'b'), family = 'Gaussian', data.site)

summary(hb.dens.m)
pp_check(hb.dens.m, ndraws = 100)
p_direction(hb.dens.m)
conditional_effects(hb.dens.m)

## Connectance

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

## Resource overlap

# Beta distribution does not give a good pp_check!
# morisita.hb.m.rich <- brm(Morisita.hb.ns ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = Beta(),
#                           prior = c(prior(normal(0,2), class = 'b'),
#                                     prior(exponential(1), class = 'sd')), data = data.morisita.hb)
# morisita.bl.m.rich <- brm(Morisita.bl.ns ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = Beta(), 
#                           prior = c(prior(normal(0,2), class = 'b'),
#                                     prior(exponential(1), class = 'sd')), data = data.morisita.bl)
# 

morisita.hb.m.s <- brm(Morisita.hb.s ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = "Gaussian",
                          prior = c(prior(normal(0,5), class = 'b')), data = data.morisita.hb, iter = 5000, warmup = 2000, 
)
morisita.bl.m.s <- brm(Morisita.bl.s ~ total_bee_dens.s + flower_dens.s + sum_nodes + Year + (1|Species) + (1|Site), family = "Gaussian", 
                          prior = c(prior(normal(0,5), class = 'b')), data = data.morisita.bl, iter = 5000, warmup = 2000, 
)

summary(morisita.hb.m.s)
summary(morisita.bl.m.s)
pp_check(morisita.hb.m.s, ndraws = 100)
pp_check(morisita.bl.m.s, ndraws = 100)
pairs(morisita.hb.m.s)
pairs(morisita.bl.m.s)
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

######################### Sensitivity ##########################
## removing datapoints without recorded interaction at a site ##
################################################################

data.both <- data.pathogen.both %>% 
  left_join(dens, by = 'Site') %>% 
  mutate(Density = ifelse(Year == 2021, 0, Density),
         Year = as.numeric(Year)) %>%
  left_join(traits, by = "Species") %>%
  left_join(flower.both %>% select(Site, Year, flower_dens), by = c('Year', 'Site')) %>%
  left_join(hb.exposure, by = c('Site', 'Year')) %>%
  left_join(bl.exposure, by = c('Site', 'Year')) %>%
  left_join(abundance.site, by = c('Site', 'Year')) %>%
  mutate(across(DWVB.abs:ABPV.abs, function(x) x/BUFFER)) %>%
  mutate(total_bee_dens = log(total_bee_abundance/flower_dens/10000 + 1)) %>%
  left_join(network.metrics.both, by = c('Site', 'Year')) %>%
  left_join(morisita.raw %>%
              rename(Morisita.hb = `Apis mellifera`,
                     Morisita.bl = `Bombus lapidarius`) %>%
              mutate(Species = sub(" agg.", "", Species)), by = join_by("Site", "Year", "Species")) %>%
  left_join(coord2, by = 'Site') %>%
  mutate(Year = as.factor(Year), Density = as.factor(Density)) %>%
  #removing bees with no interaction recorded in transects
  filter(!(is.na(Morisita.hb) & is.na(Morisita.bl)))



####  SPLITTING DATA INTO BEE GROUPS, SCALE AND FILTER THE OUTLIER
data.hb <- data.both %>% filter(Species == 'Apis mellifera') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.bb <- data.both %>% filter(Group == 'bb') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.bb.nolp <- data.both %>% filter(Group == 'bb') %>% filter(Species != 'Bombus lapidarius') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.wb <- data.both %>% filter(Group == 'wb') %>% ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

##

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

WB_models_combined_raw_noimputedint <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                                            bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                                            abpv.b = abpv.b, abpv.w = abpv.w, abpv.h = abpv.h)

save(WB_models_combined_raw_noimputedint, file = paste0('Data/Results/', date, '_models_combined_raw_no_imputed_interactions.RData'))


#### REVIEWER's COMMENT ####
### GP INSTEAD OF SITE AS RANDOM INTERCEPT ####

## per Reviewer 2 comment

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



######################
#### MODEL CHECKS ####
######################

## change the model name to explore the checks

# posterior correlation
posterior_correlation(WB_models_combined_raw[["abpv.h"]])

### PAIRS PLOTS
mcmc_pairs(as.matrix(WB_models_combined_raw[['bqcv.w']]),pars = vars(starts_with("b_") & !contains('hu')))
mcmc_pairs(as.matrix(WB_models_combined_raw[['bqcv.b']]),pars = vars(starts_with("b_hu_")))

# posterior predictive check
pp_check(WB_models_combined_raw[["dwvb"]])
pp_hurdle(WB_models_combined_raw[["dwvb"]])

# 
