source('03_data_processing.R')
source('05_plotting_functions.R')

library(brms)
library(tidybayes)
library(tidyr)
library(ggplot2)
library(rstan)
library(bayesplot)
library(bayestestR)
library(readr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### options: 
# analyse full data set and use bee species/genus as a random factor (this will obscure individual responses of different bee species)
# analyse separately hb, bb, wb (for wb random factor with genus)
# force of infection (abundance * prevalence)
# dimensionality of networks: 10.1111/ELE.14383
# sem: landscape -> networks -> prevalence


data.hb <- data %>% filter(Group == 'hb') %>% 
  mutate(across(Ann.fl:fl.cv, ~ scale(.)[,1]))

data.bb <- data %>% filter(Group == 'bb') %>% 
  mutate(across(Ann.fl:fl.cv, ~ scale(.)[,1]))

data.wb <- data %>% filter(Group == 'wb') %>% 
  mutate(across(Ann.fl:fl.cv, ~ scale(.)[,1])) %>%
  mutate(Genus = gsub(' .*', '', Species)) %>%
  # filtering out erroneous barcodes
  filter(Genus != '') %>% filter(Genus != ('Sipha')) %>% filter(Genus != ('Trypoxylon')) %>% filter(Genus != ('Orisarma')) %>% 
  filter(Genus != ('Lindenius')) %>% filter(Genus != ('Oedogonium')) %>% filter(Genus != ('Orasema')) %>% filter(Genus != ('Megalocoleus')) %>% filter(Genus != ('Bombus'))

data.wb %>% group_by(Genus) %>% summarise(n = n()) # 236 Andrena, 55 Lasioglossum
print(data.wb %>% group_by(Site) %>% summarise(n = n()), n = Inf) # Nor 174 only 6 samples, Wm1316 7


### AES ###
# separate models for each virus and bee group
# default (improper) priors

dwvb.h <- brm(bf(DWVB.abs ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.b <- brm(bf(DWVB.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.h <- brm(bf(BQCV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.b <- brm(bf(BQCV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.h <- brm(bf(ABPV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.b <- brm(bf(ABPV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.b <- brm(bf(SBV.abs ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ (Org.farm + SNH + Ann.fl) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm + SNH + Ann.fl) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

AES_dens_models <- list(dwvb.h, dwvb.b, dwvb.w, 
                        bqcv.h, bqcv.b, bqcv.w, 
                        abpv.h, abpv.b, abpv.w, 
                        sbv.h, sbv.b, sbv.w)

names(AES_dens_models) <- c('dwvb.h', 'dwvb.b', 'dwvb.w', 
                            'bqcv.h', 'bqcv.b', 'bqcv.w', 
                            'abpv.h', 'abpv.b', 'abpv.w', 
                            'sbv.h', 'sbv.b', 'sbv.w')

#save(AES_dens_models, file = 'AES_dens_models240501.RData')

## plotting and diagnosing

load('AES_dens_models240501.RData')

cs.all <- lapply(AES_dens_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

pp.all <- lapply(AES_dens_models, pp_hurdle) %>% bind_rows(.id = 'id')

dens_plot_hu(AES_dens_models[[3]], effects = 'Ann.fl:Density') # Nor578 is an outlier for org farming
dens_plot_mu(AES_dens_models[[1]], effects = 'SNH:Density')



  