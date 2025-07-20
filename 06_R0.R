### CALCULATING SPECIES R0 BASED ON FENTON ET AL. 2015###


source('05_main_data_file.R')
load('Data/true_prev_R0.RData')
load('Data/morisita_all.RData')

library(dplyr)
library(tidyr)
library(tidyverse)


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
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


### 2021

obs.prev.2021 = data2021 %>% select(Sample, Site, Year, Species, dwvb, bqcv, abpv, sbv) %>%
  mutate(Group = case_match(Species, 'Apis mellifera' ~ 'hb', 'Bombus lapidarius' ~ 'bl', 'Bombus terrestris' ~ 'bt', 'Bombus pascuorum' ~ 'bp', .default = 'wb')) %>%
  pivot_longer(dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% group_by(Group, Virus) %>%
  summarise(prev = mean(Presence), prev.sd = sd(Presence)) %>% mutate(prev.inv = logit_scaled(prev)) %>% 
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% mutate(prev.inv = round(prev.inv, 2)) %>% mutate(Group = factor(Group, levels = c('hb', 'bl', 'bp', 'bt', 'wb'))) %>%
  arrange(Group)

op1 = as.data.frame(obs.prev.2021)

op1[op1$Virus == 'dwvb',]

prior.dwvb <- c(prior(normal(-0.47,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-1.23,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-3.61,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.56,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-3.10,3), class = 'b'))

op1[op1$Virus == 'bqcv',]

prior.bqcv <- c(prior(normal(3.50,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(2.02,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-0.32,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(0.44,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-2.04,3), class = 'b'))

op1[op1$Virus == 'abpv',]

prior.abpv <- c(prior(normal(-2.45,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-0.37,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-2.14,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.37,3), class = 'b', coef = 'SpeciesBombusterrestris'),
                prior(normal(-2.67,3), class = 'b'))

# op1[op1$Virus == 'sbv',]
# 
# prior.sbv <- c(prior(normal(2.88,3), class = 'b', coef = 'SpeciesApismellifera'),
#                 prior(normal(3.50,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
#                 prior(normal(1.67,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
#                 prior(normal(2.55,3), class = 'b', coef = 'SpeciesBombusterrestris'),
#                 prior(normal(1.17,3), class = 'b'))



prev.pred.dwvb1 <- brm(dwvb ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.dwvb, sample_prior = T, data2021)
prev.pred.abpv1 <- brm(abpv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.abpv, sample_prior = T, data2021)
prev.pred.bqcv1 <- brm(bqcv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.bqcv, sample_prior = T, data2021)
# prev.pred.sbv1 <- brm(sbv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
#                       prior = prior.sbv, sample_prior = T, data2021)


true.prev.2021 <- data.frame(Species = data2021$Species, Site = data2021$Site, dwvb.true = fitted(prev.pred.dwvb1)[,1],
                             abpv.true = fitted(prev.pred.abpv1)[,1],
                             bqcv.true = fitted(prev.pred.bqcv1)[,1]) %>% group_by(Species, Site) %>% mutate(n = n())


### 2022
obs.prev.2022 = data2022 %>% select(Sample, Site, Year, Species, dwvb, bqcv, abpv, sbv) %>%
  mutate(Group = case_match(Species, 'Apis mellifera' ~ 'hb', 'Bombus lapidarius' ~ 'bl', 'Bombus pascuorum' ~ 'bp', .default = 'wb')) %>%
  pivot_longer(dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% group_by(Group, Virus) %>%
  summarise(prev = mean(Presence), prev.sd = sd(Presence)) %>% mutate(prev.inv = logit_scaled(prev)) %>% 
  mutate(prev.inv = ifelse(prev.inv == Inf, 3.5, prev.inv)) %>% mutate(prev.inv = round(prev.inv, 2)) %>% mutate(Group = factor(Group, levels = c('hb', 'bl', 'bp', 'bt', 'wb'))) %>%
  arrange(Group)

op2 = as.data.frame(obs.prev.2022)

op2[op2$Virus == 'dwvb',]

prior.dwvb <- c(prior(normal(-0.47,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(-1.04,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-2.85,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.16,3), class = 'b'))

op2[op2$Virus == 'bqcv',]

prior.bqcv <- c(prior(normal(3.50,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(3.12,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-0.19,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.27,3), class = 'b'))

op2[op2$Virus == 'abpv',]

prior.abpv <- c(prior(normal(-1.00,3), class = 'b', coef = 'SpeciesApismellifera'),
                prior(normal(0.70,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
                prior(normal(-1.34,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
                prior(normal(-1.39,3), class = 'b'))

# op2[op2$Virus == 'sbv',]
# 
# prior.sbv <- c(prior(normal(2.93,3), class = 'b', coef = 'SpeciesApismellifera'),
#                 prior(normal(3.50,3), class = 'b', coef = 'SpeciesBombuslapidarius'),
#                 prior(normal(1.59,3), class = 'b', coef = 'SpeciesBombuspascuorum'),
#                 prior(normal(1.72,3), class = 'b'))


prev.pred.dwvb2 <- brm(dwvb ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.dwvb, sample_prior = T, data2022)
prev.pred.abpv2 <- brm(abpv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.abpv, sample_prior = T, data2022)
prev.pred.bqcv2 <- brm(bqcv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
                      prior = prior.bqcv, sample_prior = T, data2022)
# prev.pred.sbv2 <- brm(sbv ~ Species - 1 + (1|Site), family = bernoulli(link = 'logit'), 
#                       prior = prior.sbv, sample_prior = T, data2022)

true.prev.2022 <- data.frame(Species = data2022$Species, Site = data2022$Site, dwvb.true = fitted(prev.pred.dwvb2)[,1],
                        abpv.true = fitted(prev.pred.abpv2)[,1],
                        bqcv.true = fitted(prev.pred.bqcv2)[,1]) %>% group_by(Species, Site) %>% mutate(n = n())

save(true.prev.2021, true.prev.2022, file = 'Data/true_prev_R0.RData')

###################
#### FUNCTIONS ####
###################

dat_into_matrix = function(abundance, true.prev, data, virus) {
  names(data) = tolower(names(data))
  names(true.prev) = tolower(names(true.prev))
  names(abundance) = tolower(names(abundance))
  t.prev <- true.prev %>% group_by(site, species) %>% rename(vir = paste0(virus, '.true')) %>% summarise(prev = mean(vir), n = n()) %>% filter(n > 2) %>% select(-n)
  load <- data %>% rename(vir = virus, vir.load = paste0(virus, '.abs')) %>% filter(vir > 0) %>% group_by(site, species)  %>% summarise(load = mean(vir.load))
  
  mat.dat <- full_join(abundance %>% select(site, species, sum), t.prev, by = c('site','species')) %>% 
    full_join(load, by = c('site','species')) %>% 
    group_by(species) %>%
    mutate(prev = ifelse(is.na(prev), mean(prev, na.rm = T), prev),
           load = ifelse(is.na(load), mean(load, na.rm = T), load)) %>%
    ungroup() %>% 
    filter(sum > 0) %>%
    mutate(load = ifelse(is.na(load), 100, load)) %>% drop_na(prev) %>% left_join(data %>% distinct(species, social), by = 'species') %>%
    mutate(social = as.numeric(ifelse(social == 0, 1, 2)))
  return(mat.dat)
}

niche_fun = function(mat.dat.list, morisita) {
  inter.niche <- list()
  for (i in 1:length(mat.dat.list)) {
    inter.niche[[i]] <- as.data.frame(morisita[[i]]) %>% rownames_to_column() %>% mutate(rowname = sub(' agg.', '', rowname)) %>% 
      column_to_rownames()
    diag(inter.niche[[i]]) = 1
    colnames(inter.niche[[i]]) <- rownames(inter.niche[[i]])
    inter.niche[[i]] <- inter.niche[[i]] %>% filter(row.names(inter.niche[[i]]) %in% mat.dat.list[[i]]$species) 
    inter.niche[[i]] <- inter.niche[[i]] %>% select(rownames(inter.niche[[i]]))
  }
  return(inter.niche)
}

niche_null_fun = function(mat.dat.list, morisita.zscore) {
  inter.niche.null <- list()
  for (i in 1:length(mat.dat.list)) {
    inter.niche.null[[i]] <- as.data.frame(morisita.zscore[[i]]) %>% rownames_to_column() %>% mutate(rowname = sub(' agg.', '', rowname)) %>% 
      column_to_rownames() %>% mutate(across(everything(), function(x) abs(x)))
    inter.niche.null[[i]][is.na(inter.niche.null[[i]])] <- 1
    colnames(inter.niche.null[[i]]) <- rownames(inter.niche.null[[i]])
    inter.niche.null[[i]] <- inter.niche.null[[i]] %>% filter(row.names(inter.niche.null[[i]]) %in% mat.dat.list[[i]]$species) 
    inter.niche.null[[i]] <- inter.niche.null[[i]] %>% select(rownames(inter.niche.null[[i]]))
  }
  return(inter.niche.null)
}


r0_func_fac <- function(data.matrix.raw, inter.niche.ex, inter.niche.null.ex, species_out = 'Apis mellifera') {
  result <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = 1)
  r0 <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  transmission <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  data.matrix <- as.matrix(data.matrix.raw[,3:6])
  # Calculate niche overlap for each pair of species
  for (j in seq_len(nrow(r0))) {
    for (i in seq_len(nrow(r0))) {
      transmission[j,i] <- rnorm(1, mean = inter.niche.ex[j,i], sd = 1/abs(inter.niche.null.ex[j,i]))
      transmission[j,j] <- data.matrix[j, 4]
      colnames(transmission) = colnames(inter.niche.ex)
      rownames(transmission) = rownames(inter.niche.ex)
      r0[j,i] <- log10(data.matrix[i, 3])/log10(data.matrix[j, 3]) * data.matrix[i, 1]/data.matrix[j, 1] * data.matrix[i, 2] * 
        (ifelse(transmission[j,i] < 0, 0, transmission[j,i])/data.matrix[j, 4])
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
        (ifelse(transmission2[o,p] < 0, 0, transmission2[o,p])/data.matrix2[o, 4])
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

mat.dat = dat_into_matrix(abundance.2021, true.prev.2021, data2021, virus = 'dwvb') %>% filter(site %in% subset)
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Apis mellifera')
  }
}


r0.dwvb.main.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.dwvb.main.2021.sp = r0.dwvb.main.2021 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.dwvb.alt.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.dwvb.alt.2021.sp = r0.dwvb.alt.2021 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))


### 2022 (full set)

mat.dat = dat_into_matrix(abundance.2022, true.prev.2022, data2022, virus = 'dwvb')
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
  mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Apis mellifera')
  }
}

r0.dwvb.main.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.dwvb.main.2022.sp = r0.dwvb.main.2022 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.dwvb.alt.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.dwvb.alt.2022.sp = r0.dwvb.alt.2022 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))


#### BQCV ####

### 2021 (subset)
mat.dat = dat_into_matrix(abundance.2021, true.prev.2021, data2021, virus = 'bqcv') %>% filter(site %in% subset)
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
  mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Apis mellifera')
  }
}


r0.bqcv.main.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.bqcv.main.2021.sp = r0.bqcv.main.2021 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.bqcv.alt.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.bqcv.alt.2021.sp = r0.bqcv.alt.2021 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))


### 2022 (full set)
mat.dat = dat_into_matrix(abundance.2022, true.prev.2022, data2022, virus = 'bqcv')
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
  mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Apis mellifera')
  }
}


r0.bqcv.main.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.bqcv.main.2022.sp = r0.bqcv.main.2022 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.bqcv.alt.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.bqcv.alt.2022.sp = r0.bqcv.alt.2022 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))


#### ABPV

### 2021 (subset)
mat.dat = dat_into_matrix(abundance.2021, true.prev.2021, data2021, virus = 'abpv') %>% filter(site %in% subset)
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
  mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Bombus lapidarius')
  }
}


r0.abpv.main.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.abpv.main.2021.sp = r0.abpv.main.2021 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.abpv.alt.2021 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.abpv.alt.2021.sp = r0.abpv.alt.2021 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

### 2022 (full set)
mat.dat = dat_into_matrix(abundance.2022, true.prev.2022, data2022, virus = 'abpv')
unique(mat.dat$species)
mat.dat.list <- split(mat.dat, mat.dat$site)

inter.niche = niche_fun(mat.dat.list, morisita.2021)
inter.niche.null = niche_null_fun(mat.dat.list, morisita.zscore2021)


for (i in 1:length(mat.dat.list)) {
  mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche.null[[i]]))
}


boot <- list()
for (s in 1:length(mat.dat.list)){
  boot[[s]] <- list()
  for (i in 1:1000) {
    boot[[s]][[i]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], inter.niche.null[[s]], species_out = 'Bombus lapidarius')
  }
}

r0.abpv.main.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
  rename(Species = rowname) %>% group_by(Site, Species) %>% 
  summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
            mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

r0.abpv.main.2022.sp = r0.abpv.main.2022 %>% group_by(Species) %>% 
  summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
            prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))

# r0.abpv.alt.2022 <- boot %>% lapply(function(x) x %>% bind_rows()) %>% bind_rows(.id = 'Site') %>% 
#   rename(Species = rowname) %>% group_by(Site, Species) %>% 
#   summarise(mean.r0 = mean(r0, na.rm = T), sd.r0 = sd(r0, na.rm = T), 
#             mean.prev.no.mh = mean(Prev.no.mh, na.rm = T), mean.prev = mean(prev, na.rm = T), sd.prev = sd(prev, na.rm = T))
# 
# r0.abpv.alt.2022.sp = r0.abpv.alt.2022 %>% group_by(Species) %>% 
#   summarise(r0 = mean(mean.r0), sd = sd(mean.r0), prev.no.mh = mean(mean.prev.no.mh, na.rm = T), 
#             prev.sd.no.mh = sd(mean.prev.no.mh, na.rm = T), prev = mean(mean.prev, na.rm = T), sd.prev = sd(prev, na.rm = T))


#######

r0.main.host <- list(r0.dwvb.main.2021 = r0.dwvb.main.2021, r0.dwvb.main.2022 = r0.dwvb.main.2022, 
                r0.bqcv.main.2021 = r0.bqcv.main.2021, r0.bqcv.main.2022 = r0.bqcv.main.2022, 
                r0.abpv.main.2021 = r0.abpv.main.2021, r0.abpv.main.2022 = r0.abpv.main.2022)

r0.main.host.sp <- list(r0.dwvb.main.2021.sp = r0.dwvb.main.2021.sp, r0.dwvb.main.2022.sp = r0.dwvb.main.2022.sp, 
                        r0.bqcv.main.2021.sp = r0.bqcv.main.2021.sp, r0.bqcv.main.2022.sp = r0.bqcv.main.2022.sp, 
                        r0.abpv.main.2021.sp = r0.abpv.main.2021.sp, r0.abpv.main.2022.sp = r0.abpv.main.2022.sp)

# r0.alt.host = list(r0.dwvb.alt.2021 = r0.dwvb.alt.2021, r0.dwvb.alt.2022 = r0.dwvb.alt.2022, 
#                    r0.bqcv.alt.2021 = r0.bqcv.alt.2021, r0.bqcv.alt.2022 = r0.bqcv.alt.2022, 
#                    r0.abpv.alt.2021 = r0.abpv.alt.2021, r0.abpv.alt.2022 = r0.abpv.alt.2022)
# 
# r0.alt.host.sp <- list(r0.dwvb.alt.2021.sp = r0.dwvb.alt.2021.sp, r0.dwvb.alt.2022.sp = r0.dwvb.alt.2022.sp, 
#                         r0.bqcv.alt.2021.sp = r0.bqcv.alt.2021.sp, r0.bqcv.alt.2022.sp = r0.bqcv.alt.2022.sp, 
#                         r0.abpv.alt.2021.sp = r0.abpv.alt.2021.sp, r0.abpv.alt.2022.sp = r0.abpv.alt.2022.sp)

#save(r0.main.host,r0.main.host.sp, r0.alt.host, r0.alt.host.sp, file = 'Data/250260R0_results.RData')



### comparing simulated prevalence
# data frames from 8_plotting.R

dwvb.hb.b = sim.ms %>% select(Species, sim.m_dwvb, prev.m_dwvb) %>% pivot_longer(cols = sim.m_dwvb:prev.m_dwvb, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Apis mellifera')

mod.dwvb.hb = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                               prior(normal(0,1), class = 'Intercept')), data = dwvb.hb.b)
mod.dwvb.hb.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = dwvb.hb.b)
bayes_factor(mod.dwvb.hb.null , mod.dwvb.hb)
pp_check(mod.dwvb.hb, ndraws = 100)
hypothesis(mod.dwvb.hb, "Typesim.m_dwvb = 0")

dwvb.bl.b = sim.ms.alt %>% select(Species, sim.m_dwvb) %>% left_join(sim.ms %>% select(Species, prev.m_dwvb) %>%
                                                                     mutate(Species = sub('\n',' ', Species)), by = 'Species') %>%
  pivot_longer(cols = sim.m_dwvb:prev.m_dwvb, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Bombus lapidarius')

mod.dwvb.bl = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                                            prior(normal(0,1), class = 'Intercept')), data = dwvb.bl.b)
mod.dwvb.bl.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = dwvb.bl.b)
bayes_factor(mod.dwvb.bl.null, mod.dwvb.bl)
bayes_factor(mod.dwvb.hb, mod.dwvb.bl)
p_direction(mod.dwvb.hb)

####

bqcv.hb.b = sim.ms %>% select(Species, sim.m_bqcv, prev.m_bqcv) %>% pivot_longer(cols = sim.m_bqcv:prev.m_bqcv, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Apis mellifera')

mod.bqcv.hb = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                                            prior(normal(0,1), class = 'Intercept')), data = bqcv.hb.b)
mod.bqcv.hb.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = bqcv.hb.b)
bayes_factor(mod.bqcv.hb.null , mod.bqcv.hb)
pp_check(mod.bqcv.hb, ndraws = 100)
hypothesis(mod.bqcv.hb, "Typesim.m_bqcv = 0")

bqcv.bl.b = sim.ms.alt %>% select(Species, sim.m_bqcv) %>% left_join(sim.ms %>% select(Species, prev.m_bqcv) %>%
                                                                       mutate(Species = sub('\n',' ', Species)), by = 'Species') %>%
  pivot_longer(cols = sim.m_bqcv:prev.m_bqcv, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Bombus lapidarius')

mod.bqcv.bl = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                                            prior(normal(0,1), class = 'Intercept')), data = bqcv.bl.b)
mod.bqcv.bl.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = bqcv.bl.b)
bayes_factor(mod.bqcv.bl.null, mod.bqcv.bl)
bayes_factor(mod.bqcv.hb, mod.bqcv.bl)
bayes_factor(mod.bqcv.bl, mod.bqcv.hb)
p_direction(mod.bqcv.hb)

####

abpv.hb.b = sim.ms.alt %>% select(Species, sim.m_abpv) %>% left_join(sim.ms %>% select(Species, prev.m_abpv) %>%
                                                                       mutate(Species = sub('\n',' ', Species)), by = 'Species') %>%
  pivot_longer(cols = sim.m_abpv:prev.m_abpv, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Apis mellifera')

mod.abpv.hb = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                                            prior(normal(0,1), class = 'Intercept')), data = abpv.hb.b)
mod.abpv.hb.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = abpv.hb.b)
bayes_factor(mod.abpv.hb.null , mod.abpv.hb)
pp_check(mod.abpv.hb, ndraws = 100)
hypothesis(mod.abpv.hb, "Typesim.m_abpv = 0")

abpv.bl.b = sim.ms %>% select(Species, sim.m_abpv, prev.m_abpv) %>%
  pivot_longer(cols = sim.m_abpv:prev.m_abpv, names_to = 'Type', values_to = 'Prev') %>%
  mutate(Species = sub('\n',' ', Species)) %>% filter(Species != 'Bombus lapidarius')

mod.abpv.bl = brm(Prev ~ Type, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'b'),
                                                                            prior(normal(0,1), class = 'Intercept')), data = abpv.bl.b)
mod.abpv.bl.null = brm(Prev ~ 1, family = Beta(), sample_prior = T, prior = c(prior(normal(0,1), class = 'Intercept')), data = abpv.bl.b)
bayes_factor(mod.abpv.bl.null, mod.abpv.bl)
bayes_factor(mod.abpv.bl, mod.abpv.hb)
p_direction(mod.abpv.hb)
