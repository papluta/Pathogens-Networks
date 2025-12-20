### CALCULATING SPECIES R0 BASED ON FENTON ET AL. 2015###
library(glue)
library(brms)
library(tidybayes)
library(rstan)
library(bridgesampling)
set.seed(99)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('05_main_data_file.R') # reruns also previous scripts (01 - 05)

# re-load data if available
#load('Data/251126_true_prev_R0.RData')


subset <- c('Goe1392', 'Goe1425', 'Goe235', 'Goe288', 'Goe47', 'Goe595',
            'Gos1', 'Gos2', 'Nor1', 'Nor1070', 'Nor1145', 'Nor264', 'Nor508', 'Nor918', 'WM1249', 'WM630')


data.both2 = data.both %>% 
  mutate(across(DWVB.abs:ABPV.abs, function(x) x * BUFFER)) # restore the viral load per bee, not per uL of buffer

data2021 <- data.both2 %>% filter(Year == 2021)
data2022 <- data.both2 %>% filter(Year == 2022)

abundance.2021 <- abundance.both %>% filter(Year == 2021) %>% select(-Year)
abundance.2022 <- abundance.both %>% filter(Year == 2022) %>% select(-Year)


###### ESTIMATING ADJUSTED PREVALENCE

### 2021

obs.prev.2021 = data2021 %>% select(Sample, Site, Year, Species, dwvb, bqcv, abpv) %>%
  mutate(Group = case_match(Species, 
                            'Apis mellifera' ~ 'hb', 
                            'Bombus lapidarius' ~ 'bl', 
                            'Bombus terrestris' ~ 'bt', 
                            'Bombus pascuorum' ~ 'bp', 
                            .default = 'wb')) %>%
  pivot_longer(dwvb:abpv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Group, Virus) %>%
  summarise(prev = mean(Presence), prev.sd = sd(Presence)) %>% 
  mutate(prev.inv = logit_scaled(prev)) %>% 
  # assingning a value to 1
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% # adding a 0.99 prevalence approximation to swap the Inf
  mutate(prev.inv = round(prev.inv, 2)) %>% 
  mutate(Group = factor(Group, levels = c('hb', 'bl', 'bp', 'bt', 'wb'))) %>%
  arrange(Group)

op1 = as.data.frame(obs.prev.2021)

op1[op1$Virus == 'dwvb',]

# cannot assign prior values dynamically because prior function evaluates them as a character, 
# adding them manually from the tibble above

prior.dwvb <- c(prior(normal(-0.47,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-1.22,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-3.31,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.54,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-3.07,3), class = 'b'))

op1[op1$Virus == 'bqcv',]

prior.bqcv <- c(prior(normal(3.50,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(2.09,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-0.11,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(0.48,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-2.01,3), class = 'b'))

op1[op1$Virus == 'abpv',]

prior.abpv <- c(prior(normal(-2.48,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-0.39,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-1.97,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.35,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-2.64,3), class = 'b'))


prev.pred.dwvb1 <- brm(dwvb ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.dwvb, sample_prior = T, data2021)
prev.pred.abpv1 <- brm(abpv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.abpv, sample_prior = T, data2021)
prev.pred.bqcv1 <- brm(bqcv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.bqcv, sample_prior = T, data2021)

true.prev.2021 <- data.frame(Species = data2021$Species, Site = data2021$Site, dwvb.true = fitted(prev.pred.dwvb1)[,1],
                             abpv.true = fitted(prev.pred.abpv1)[,1],
                             bqcv.true = fitted(prev.pred.bqcv1)[,1]) %>% 
  group_by(Species, Site) %>% mutate(n = n())


### 2022
obs.prev.2022 = data2022 %>% 
  select(Sample, Site, Year, Species, dwvb, bqcv, abpv) %>%
  mutate(Group = case_match(Species, 
                            'Apis mellifera' ~ 'hb', 
                            'Bombus lapidarius' ~ 'bl', 
                            'Bombus pascuorum' ~ 'bp', 
                            .default = 'wb')) %>%
  pivot_longer(dwvb:abpv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Group, Virus) %>%
  summarise(prev = mean(Presence), prev.sd = sd(Presence)) %>% 
  mutate(prev.inv = logit_scaled(prev)) %>% 
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% # adding a 0.99 prevalence approximation to swap the Inf
  mutate(prev.inv = round(prev.inv, 2)) %>% 
  mutate(Group = factor(Group, levels = c('hb', 'bl', 'bp', 'bt', 'wb'))) %>%
  arrange(Group)

op2 = as.data.frame(obs.prev.2022)

op2[op2$Virus == 'dwvb',]

prior.dwvb <- c(prior(normal(-0.45,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-1.04,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-2.85,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.20,3), class = 'b'))

op2[op2$Virus == 'bqcv',]

prior.bqcv <- c(prior(normal(3.50,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(3.12,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-0.19,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.24,3), class = 'b'))

op2[op2$Virus == 'abpv',]

prior.abpv <- c(prior(normal(-1.01,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(0.70,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-1.33,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.38,3), class = 'b'))

prev.pred.dwvb2 <- brm(dwvb ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.dwvb, sample_prior = T, data2022)
prev.pred.abpv2 <- brm(abpv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.abpv, sample_prior = T, data2022)
prev.pred.bqcv2 <- brm(bqcv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.bqcv, sample_prior = T, data2022)

true.prev.2022 <- data.frame(Species = data2022$Species, Site = data2022$Site, 
                             dwvb.true = fitted(prev.pred.dwvb2)[,1],
                        abpv.true = fitted(prev.pred.abpv2)[,1],
                        bqcv.true = fitted(prev.pred.bqcv2)[,1]) %>% 
  group_by(Species, Site) %>% 
  mutate(n = n())

save(true.prev.2021, true.prev.2022, file = paste0('Data/', date, '_true_prev_R0.RData'))


### DWV-B

### 2021 (subset)
  
R0.dwvb.noAM.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "dwvb", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Apis mellifera") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Apis mellifera")

R0.dwvb.noLP.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "dwvb", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Lasioglossum pauxillum")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Lasioglossum pauxillum")
  
## 2022

R0.dwvb.noAM.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "dwvb", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Apis mellifera") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Apis mellifera")

R0.dwvb.noLP.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "dwvb", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Lasioglossum pauxillum")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Lasioglossum pauxillum")

### BQCV

## 2021

R0.bqcv.noAM.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "bqcv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Apis mellifera") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Apis mellifera")

R0.bqcv.noBL.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "bqcv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Bombus lapidarius")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Bombus lapidarius")

## 2022

R0.bqcv.noAM.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "bqcv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Apis mellifera") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Apis mellifera")

R0.bqcv.noBL.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "bqcv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Bombus lapidarius")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Bombus lapidarius")

### ABPV

## 2021

R0.abpv.noAM.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "abpv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Andrena minutula") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Andrena minutula")

R0.abpv.noBL.2021 <- r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "abpv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Bombus lapidarius")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Bombus lapidarius")

## 2022

R0.abpv.noAM.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "abpv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Andrena minutula") %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Andrena minutula")

R0.abpv.noBL.2022 <- r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "abpv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Bombus lapidarius")  %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Bombus lapidarius")


#######

r0.results.long <- list(dwvb = bind_rows(R0.dwvb.noAM.2021, R0.dwvb.noAM.2022, R0.dwvb.noLP.2021, R0.dwvb.noLP.2022), 
                     bqcv = bind_rows(R0.bqcv.noAM.2021, R0.bqcv.noAM.2022, R0.bqcv.noBL.2021, R0.bqcv.noBL.2022), 
                     abpv = bind_rows(R0.abpv.noAM.2021, R0.abpv.noAM.2022, R0.abpv.noBL.2021, R0.abpv.noBL.2022)
                     )

r0.dwvb <- r0.results.long$dwvb %>% 
  #removing duplicated values for R0, filtering by either main host will bear the same result
  filter(main.host == "Apis mellifera") %>% 
  group_by(Species) %>% 
  summarise(r0_mean = mean(r0), r0_sd = sd(r0), n_networks = n())

r0.sim.prev.dwvb <- r0.results.long$dwvb %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  rename(sim_prev_noAM = `Apis mellifera`, sim_prev_noBL = `Lasioglossum pauxillum`) %>%
  group_by(Species) %>% 
  summarise(obs_prev_mean = mean(prev), obs_prev_sd = sd(prev),
            sim_prev_noAM_mean = mean(sim_prev_noAM), sim_prev_noAM_sd = sd(sim_prev_noAM),
            sim_prev_noBL_mean = mean(sim_prev_noBL), sim_prev_noBL_sd = sd(sim_prev_noBL))%>%
  left_join(r0.dwvb, by = "Species")


r0.bqcv <- r0.results.long$bqcv %>% 
  #removing duplicated values for R0, filtering by either main host will bear the same result
  filter(main.host == "Apis mellifera") %>% group_by(Species) %>% 
  summarise(r0_mean = mean(r0), r0_sd = sd(r0), n_networks = n())

r0.sim.prev.bqcv <- r0.results.long$bqcv %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  rename(sim_prev_noAM = `Apis mellifera`, sim_prev_noBL = `Bombus lapidarius`) %>%
  group_by(Species) %>% 
  summarise(obs_prev_mean = mean(prev), obs_prev_sd = sd(prev),
            sim_prev_noAM_mean = mean(sim_prev_noAM), sim_prev_noAM_sd = sd(sim_prev_noAM),
            sim_prev_noBL_mean = mean(sim_prev_noBL), sim_prev_noBL_sd = sd(sim_prev_noBL))%>%
  left_join(r0.bqcv, by = "Species")


r0.abpv <- r0.results.long$abpv %>% 
  #removing duplicated values for R0, filtering by either main host will give the same result
  filter(main.host == "Bombus lapidarius") %>% group_by(Species) %>% 
  summarise(r0_mean = mean(r0), r0_sd = sd(r0), n_networks = n())

r0.sim.prev.abpv <- r0.results.long$abpv %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  rename(sim_prev_noAM = `Andrena minutula`, sim_prev_noBL = `Bombus lapidarius`) %>%
  group_by(Species) %>% 
  summarise(obs_prev_mean = mean(prev), obs_prev_sd = sd(prev),
            sim_prev_noAM_mean = mean(sim_prev_noAM), sim_prev_noAM_sd = sd(sim_prev_noAM),
            sim_prev_noBL_mean = mean(sim_prev_noBL), sim_prev_noBL_sd = sd(sim_prev_noBL)) %>%
  left_join(r0.abpv, by = "Species")


save(r0.results.long, r0.sim.prev.dwvb, r0.sim.prev.bqcv, r0.sim.prev.abpv, file = "Data/Results/251126_R0_results.RData")

### comparing simulated prevalence

#### DWV-B
dwvb.data <- r0.results.long$dwvb %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  pivot_longer(cols = c(prev, `Apis mellifera`, `Lasioglossum pauxillum`), names_to = "Type", values_to = "Prev")  %>%
  distinct(Site, Year, Type, .keep_all = T) %>%
  mutate(Type = case_when(
    Type == "prev" ~ "observed",
    Type == "Apis mellifera" ~ "no AM",
    Type == "Lasioglossum pauxillum" ~ "no LP",
    TRUE ~ NA
  ) %>% factor(., levels = c("observed", "no AM", "no LP")),
  Year = as.factor(Year),
  Site = as.factor(Site)) %>%
  filter(!is.na(Prev))

## Honeybee

dwvb.data.hb <- dwvb.data %>% filter(Type != "no LP", Species != "Apis mellifera")

mod.dwvb.hb = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = dwvb.data.hb)
mod.dwvb.hb.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = dwvb.data.hb)

bridge_full <- bridge_sampler(mod.dwvb.hb)
bridge_null <- bridge_sampler(mod.dwvb.hb.null)

bf.dwvb.hb <- bayes_factor(bridge_full, bridge_null)

## L. pauxillum

dwvb.data.lp <- dwvb.data %>% filter(Type != "no AM", Species != "Lasioglossum pauxillum")

mod.dwvb.lp = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = dwvb.data.lp)
mod.dwvb.lp.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = dwvb.data.lp)

bridge_full <- bridge_sampler(mod.dwvb.lp)
bridge_null <- bridge_sampler(mod.dwvb.lp.null)

bf.dwvb.lp <- bayes_factor(bridge_full, bridge_null)

#### BQCV

bqcv.data <- r0.results.long$bqcv %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  pivot_longer(cols = c(prev, `Apis mellifera`, `Bombus lapidarius`), names_to = "Type", values_to = "Prev")  %>%
  distinct(Site, Year, Type, .keep_all = T) %>%
  mutate(Type = case_when(
    Type == "prev" ~ "observed",
    Type == "Apis mellifera" ~ "no AM",
    Type == "Bombus lapidarius" ~ "no BL",
    TRUE ~ NA
  ) %>% factor(., levels = c("observed", "no AM", "no BL")),
  Year = as.factor(Year),
  Site = as.factor(Site)) %>%
  filter(!is.na(Prev))

## Honeybee

bqcv.data.hb <- bqcv.data %>% filter(Type != "no BL", Species != "Apis mellifera")

mod.bqcv.hb = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = bqcv.data.hb)
mod.bqcv.hb.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = bqcv.data.hb)

bridge_full <- bridge_sampler(mod.bqcv.hb)
bridge_null <- bridge_sampler(mod.bqcv.hb.null)

bf.bqcv.hb <- bayes_factor(bridge_full, bridge_null)

## B. lapidarius

bqcv.data.bl <- bqcv.data %>% filter(Type != "no AM", Species != "Bombus lapidarius")

mod.bqcv.bl = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = bqcv.data.bl)
mod.bqcv.bl.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = bqcv.data.bl)

bridge_full <- bridge_sampler(mod.bqcv.bl)
bridge_null <- bridge_sampler(mod.bqcv.bl.null)

bf.bqcv.bl <- bayes_factor(bridge_full, bridge_null)

#### ABPV

abpv.data <- r0.results.long$abpv %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  pivot_longer(cols = c(prev, `Andrena minutula`, `Bombus lapidarius`), names_to = "Type", values_to = "Prev")  %>%
  distinct(Site, Year, Type, .keep_all = T) %>%
  mutate(Type = case_when(
    Type == "prev" ~ "observed",
    Type == "Andrena minutula" ~ "no AM",
    Type == "Bombus lapidarius" ~ "no BL",
    TRUE ~ NA
  ) %>% factor(., levels = c("observed", "no AM", "no BL")),
  Year = as.factor(Year),
  Site = as.factor(Site)) %>%
  filter(!is.na(Prev))

## Honeybee

abpv.data.am <- abpv.data %>% filter(Type != "no BL", Species != "Andrena minutula")

mod.abpv.am = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = abpv.data.am)
mod.abpv.am.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = abpv.data.am)

bridge_full <- bridge_sampler(mod.abpv.am)
bridge_null <- bridge_sampler(mod.abpv.am.null)

bf.abpv.am <- bayes_factor(bridge_full, bridge_null)

## B. lapidarius

abpv.data.bl <- abpv.data %>% filter(Type != "no AM", Species != "Bombus lapidarius")

mod.abpv.bl = brm(Prev ~ Type + Year + (1|Species) + (1|Site), family = Beta(), 
                  sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                              prior(normal(0,1), class = 'Intercept'),
                                              prior(exponential(1), class = "sd")), 
                  save_pars = save_pars(all = TRUE), data = abpv.data.bl)
mod.abpv.bl.null = brm(Prev ~ Year + (1|Species) + (1|Site), family = Beta(), 
                       sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                   prior(normal(0,1), class = 'Intercept'),
                                                   prior(exponential(1), class = "sd")), 
                       save_pars = save_pars(all = TRUE), data = abpv.data.bl)

bridge_full <- bridge_sampler(mod.abpv.bl)
bridge_null <- bridge_sampler(mod.abpv.bl.null)

bf.abpv.bl <- bayes_factor(bridge_full, bridge_null)

bf.list <- list(bf.dwvb.hb = bf.dwvb.hb, bf.dwvb.lp = bf.dwvb.lp,
                bf.bqcv.hb = bf.bqcv.hb, bf.bqcv.bl = bf.bqcv.bl,
                bf.abpv.bl = bf.abpv.bl, bf.abpv.am = bf.abpv.am) 

save(bf.list, file = paste0("Data/Results/",date,"_bayes_factor.RData"))

#########################################
####  SENSITIVITY ANALYSIS SOCIALITY ####
#########################################

sociality_param = c(1,2,5,10)


### DWVB

## 2021


sens.R0.dwvb.2021 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "dwvb", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Apis mellifera",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Apis mellifera")), sociality_param) %>%
  bind_rows(.id = "Sociality")

## 2022

sens.R0.dwvb.2022 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "dwvb", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Apis mellifera",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Apis mellifera")), sociality_param) %>%
  bind_rows(.id = "Sociality")

### BQCV

## 2021

sens.R0.bqcv.2021 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "bqcv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Apis mellifera",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Apis mellifera")), sociality_param) %>%
  bind_rows(.id = "Sociality")

## 2022

sens.R0.bqcv.2022 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "bqcv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Apis mellifera",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Apis mellifera")), sociality_param) %>%
  bind_rows(.id = "Sociality")


### ABPV

## 2021

sens.R0.abpv.2021 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2021,
                                true.prev = true.prev.2021,
                                data = data2021,
                                virus = "abpv", 
                                is_subset = TRUE, 
                                morisita = morisita.2021, 
                                species_out = "Andrena minutula",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2021, main.host = "Andrena minutula")), sociality_param) %>%
  bind_rows(.id = "Sociality")

## 2022

sens.R0.abpv.2022 <- set_names(map(sociality_param, ~ r0_wrapper(abundance = abundance.2022,
                                true.prev = true.prev.2022,
                                data = data2022,
                                virus = "abpv", 
                                is_subset = F, 
                                morisita = morisita.2022, 
                                species_out = "Andrena minutula",
                                sociality = .x) %>% 
  bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% mutate(Year = 2022, main.host = "Andrena minutula")), sociality_param) %>%
  bind_rows(.id = "Sociality")


#######

r0.results.sensitivity <- list(dwvb = bind_rows(sens.R0.dwvb.2021, sens.R0.dwvb.2022), 
                               bqcv = bind_rows(sens.R0.bqcv.2021, sens.R0.bqcv.2022), 
                               abpv = bind_rows(sens.R0.abpv.2021, sens.R0.abpv.2022))

save(r0.results.sensitivity, file = paste0("Data/Results/", date, "_r0_sensitivity.RData"))
