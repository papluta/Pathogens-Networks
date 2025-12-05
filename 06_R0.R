### CALCULATING SPECIES R0 BASED ON FENTON ET AL. 2015###


source('05_main_data_file.R') # reruns also previous scripts (01 - 05)

# re-load data if available
load('Data/251126_true_prev_R0.RData')
#load('Data/morisita_all.RData')

library(dplyr)
library(tidyr)
library(tidyverse)
library(glue)


data.both2 = data.both %>% 
  mutate(across(DWVB.abs:ABPV.abs, function(x) x * BUFFER))

data2021 <- data.both2 %>% filter(Year == 2021)
data2022 <- data.both2 %>% filter(Year == 2022)

abundance.2021 <- abundance.both %>% filter(Year == 2021) %>% select(-Year)
abundance.2022 <- abundance.both %>% filter(Year == 2022) %>% select(-Year)


###### ESTIMATING ADJUSTED PREVALENCE
library(brms)
library(tidyverse)
library(tidybayes)
library(tidyr)
library(ggplot2)
library(rstan)
library(bridgesampling)
set.seed(99)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% 
  mutate(prev.inv = round(prev.inv, 2)) %>% 
  mutate(Group = factor(Group, levels = c('hb', 'bl', 'bp', 'bt', 'wb'))) %>%
  arrange(Group)

op1 = as.data.frame(obs.prev.2021)

op1[op1$Virus == 'dwvb',]

# cannot assign prior values dynamicaly, adding them manually
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

prior.abpv <- c(prior(normal(-2.45,3), class = 'b', coef = 'SpeciesApismellifera'),
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
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% 
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

save(true.prev.2021, true.prev.2022, file = 'Data/251126_true_prev_R0.RData')

###################
#### FUNCTIONS ####
###################

dat_into_matrix = function(abundance, true.prev, data, virus, sociality = 2) {
  names(data) = tolower(names(data))
  names(true.prev) = tolower(names(true.prev))
  names(abundance) = tolower(names(abundance))
  t.prev <- true.prev %>% 
    group_by(site, species) %>% 
    rename(vir = paste0(virus, '.true')) %>% 
    summarise(prev = mean(vir), n = n()) %>% 
    filter(n > 2) %>% 
    select(-n)
  load <- data %>% 
    rename(vir = virus, vir.load = paste0(virus, '.abs')) %>% 
    filter(vir > 0) %>% 
    group_by(site, species) %>% 
    summarise(load = mean(vir.load))
  
  mat.dat <- full_join(abundance %>% 
                         select(site, species, bee_abundance), t.prev, by = c('site','species')) %>% 
    full_join(load, by = c('site','species')) %>% 
    group_by(species) %>%
    mutate(prev = ifelse(is.na(prev), mean(prev, na.rm = T), prev),
           load = ifelse(is.na(load), mean(load, na.rm = T), load)) %>%
    ungroup() %>% 
    filter(bee_abundance > 0) %>%
    mutate(load = ifelse(is.na(load), 100, load)) %>% 
    drop_na(prev) %>% 
    left_join(data %>% distinct(species, social), by = 'species') %>%
    mutate(social = as.numeric(ifelse(social == 0, 1, sociality)))
  return(mat.dat)
}

niche_fun = function(mat.dat.list, morisita) {
  inter.niche <- list()
  inter.niche2 <- list()
  for (i in 1:length(mat.dat.list)) {
    inter.niche[[i]] <- as.data.frame(morisita[[i]]) %>% rownames_to_column() %>% mutate(rowname = sub(' agg.', '', rowname)) %>% 
      column_to_rownames()
    diag(inter.niche[[i]]) = 1
    colnames(inter.niche[[i]]) <- rownames(inter.niche[[i]])
    inter.niche2[[i]] <- inter.niche[[i]] %>% filter(row.names(inter.niche[[i]]) %in% mat.dat.list[[i]]$species) 
    inter.niche2[[i]] <- inter.niche2[[i]] %>% select(rownames(inter.niche2[[i]]))
    excluded_species <- setdiff(colnames(inter.niche[[i]]), colnames(inter.niche2[[i]])) 
    message(glue("Species removed from site {i}: {paste0(excluded_species, collapse = ',')}"))
  }
  return(inter.niche2)
}

# niche_null_fun = function(mat.dat.list, morisita.zscore) {
#   inter.niche.null <- list()
#   for (i in 1:length(mat.dat.list)) {
#     inter.niche.null[[i]] <- as.data.frame(morisita.zscore[[i]]) %>% rownames_to_column() %>% mutate(rowname = sub(' agg.', '', rowname)) %>% 
#       column_to_rownames() %>% mutate(across(everything(), function(x) abs(x)))
#     inter.niche.null[[i]][is.na(inter.niche.null[[i]])] <- 1
#     colnames(inter.niche.null[[i]]) <- rownames(inter.niche.null[[i]])
#     inter.niche.null[[i]] <- inter.niche.null[[i]] %>% filter(row.names(inter.niche.null[[i]]) %in% mat.dat.list[[i]]$species) 
#     inter.niche.null[[i]] <- inter.niche.null[[i]] %>% select(rownames(inter.niche.null[[i]]))
#   }
#   return(inter.niche.null)
# }


r0_func_fac <- function(data.matrix.raw, inter.niche.ex, species_out = 'Apis mellifera') {
  result <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = 1)
  r0 <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  transmission <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  data.matrix <- as.matrix(data.matrix.raw[,3:6])
  
  # Calculate niche overlap for each pair of species
  for (j in seq_len(nrow(r0))) {
    for (i in seq_len(nrow(r0))) {
      transmission[j,i] <- inter.niche.ex[j,i]
      transmission[j,j] <- data.matrix[j, 4]
      colnames(transmission) = colnames(inter.niche.ex)
      rownames(transmission) = rownames(inter.niche.ex)
      r0[j,i] <- log10(data.matrix[i, 3])/log10(data.matrix[j, 3]) * data.matrix[i, 1]/data.matrix[j, 1] * data.matrix[i, 2] * 
        (transmission[j,i]/data.matrix[j, 4])
      #r0[j,j] <- 1
    }
    result[j,1] <- 1 / (((1 - data.matrix[j, 2])/data.matrix[j, 2]) * sum(r0[j,], na.rm = T))
    rownames(result) <- rownames(inter.niche.ex)
    colnames(result) <- 'r0'
    #result_com = list(r0 = result, transmission = transmission)
  }
  # remove main host
  data.matrix2 = data.matrix.raw %>% filter(!species %in% species_out)
  transmission2 = as.data.frame(transmission) %>% filter(row.names(transmission) %in% data.matrix2$species)
  transmission2 = as.data.frame(transmission2) %>% select(rownames(transmission2))
  r02 = as.data.frame(result)  %>% filter(row.names(result) %in% data.matrix2$species)
  pred.prev <- matrix(NA, nrow = nrow(data.matrix2), ncol = 1)
  mat <- matrix(NA, nrow = nrow(data.matrix2), ncol = nrow(data.matrix2))
  data.matrix2 <- as.matrix(data.matrix2[,3:6])
  ## run the code solving for prevalence
  for (o in seq_len(nrow(mat))) {
    for (p in seq_len(nrow(mat))) {
      mat[o,p] <- log10(data.matrix2[p, 3])/log10(data.matrix2[o, 3]) * data.matrix2[p, 1]/data.matrix2[o, 1] * data.matrix2[p, 2] * 
        (transmission2[o,p]/data.matrix2[o, 4])
    }
    pred.prev[o,1] <- r02[o,1] * sum(mat[o,], na.rm = T) / (1 + r02[o,1] * sum(mat[o,], na.rm = T))
    rownames(pred.prev) <- rownames(transmission2)
    #colnames(pred.prev) <- 'prev_no_hb'
  }
  result_com = as.data.frame(result) %>% rownames_to_column() %>% left_join(as.data.frame(pred.prev) %>% rownames_to_column(), by = 'rowname') %>% rename(Prev.no.mh = V1) %>%
    left_join(data.matrix.raw %>% select(species, prev), by = join_by('rowname' == 'species'))
  return(result_com)
}


subset <- c('Goe1392', 'Goe1425', 'Goe235', 'Goe288', 'Goe47', 'Goe595',
           'Gos1', 'Gos2', 'Nor1', 'Nor1070', 'Nor1145', 'Nor264', 'Nor508', 'Nor918', 'WM1249', 'WM630')

### DWV-B

### 2021 (subset)

r0_wrapper <- function(abundance, true.prev, data, virus, is_subset, morisita, species_out, sociality = 2) {
  
  mat.dat = dat_into_matrix(abundance, true.prev, data, virus, sociality)
  
  if (is_subset) {
    
    mat.dat = mat.dat %>% filter(site %in% subset)
    
  }
  
  mat.dat.list <- split(mat.dat, mat.dat$site)
  
  inter.niche = niche_fun(mat.dat.list, morisita)
  
  for (i in 1:length(mat.dat.list)) {
    mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche[[i]]))
  }
  
  r0 <- list()
  for (s in 1:length(mat.dat.list)){
    r0[[s]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], species_out = species_out)
  }
  
  return(r0)
}

### DWVB

## 2021
  
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

r0.dwvb <- r0.results.long$dwvb %>% filter(main.host == "Apis mellifera") %>% #removing duplicated values for R0
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


r0.bqcv <- r0.results.long$bqcv %>% filter(main.host == "Apis mellifera") %>% group_by(Species) %>% 
  summarise(r0_mean = mean(r0), r0_sd = sd(r0), n_networks = n())

r0.sim.prev.bqcv <- r0.results.long$bqcv %>% 
  pivot_wider(names_from = main.host, values_from = Prev.no.mh) %>%
  rename(sim_prev_noAM = `Apis mellifera`, sim_prev_noBL = `Bombus lapidarius`) %>%
  group_by(Species) %>% 
  summarise(obs_prev_mean = mean(prev), obs_prev_sd = sd(prev),
            sim_prev_noAM_mean = mean(sim_prev_noAM), sim_prev_noAM_sd = sd(sim_prev_noAM),
            sim_prev_noBL_mean = mean(sim_prev_noBL), sim_prev_noBL_sd = sd(sim_prev_noBL))%>%
  left_join(r0.bqcv, by = "Species")


r0.abpv <- r0.results.long$abpv %>% filter(main.host == "Bombus lapidarius") %>% group_by(Species) %>% 
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

