source('02_land_extrapolation_both_years.R')

library(bipartite)
library(purrr)

trans2021.2 <- trans2021 |> filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') |> filter(Bee_species != 'no_identificiation')
trans2022.2 <- trans2022 |> filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee')

trans2021.2 %>% filter(Site %in% c('Goe595')) %>% distinct(Bee_species)
trans2022.2 %>% filter(Site %in% c('Nor174')) %>% distinct(Bee_species)

levels(as.factor(trans2021.2$Bee_group))
levels(as.factor(trans2021.2$Bee_species))
levels(as.factor(trans2022.2$Bee_group))
levels(as.factor(trans2022.2$Bee_species))


adj.matrix.2021 <- frame2webs(trans2021.2, varnames = c("visited_flower_species", "Bee_species",
                                                        "Site"), emptylist = TRUE)
adj.matrix.2022 <- frame2webs(trans2022.2, varnames = c("visited_flower_species", "Bee_species",
                                            "Site"), emptylist = TRUE)


net.met2021 <- adj.matrix.2021 |> lapply(networklevel, index=c("connectance", 
                                                                  "weighted NODF", 
                                                                  "niche overlap"), level = 'higher', 
                                            weighted = TRUE)
net.met2022 <- adj.matrix.2022 |> lapply(networklevel, index=c("connectance", 
                                                                "weighted NODF", 
                                                                "niche overlap"), level = 'higher', 
                                          weighted = TRUE)

net.met2021_species <- adj.matrix.2021 |> lapply(specieslevel, index=c("d", 
                                                                       "normalised degree", 
                                                                       "closeness", "betweenness"), level = 'higher') 

net.met2022_species <- adj.matrix.2022 |> lapply(specieslevel, index=c("d", 
                                                                       "normalised degree", 
                                                                       "closeness", "betweenness"), level = 'higher')
d2021 <- lapply(net.met2021_species, function(x) x %>% rownames_to_column()) %>% bind_rows(.id = 'Site') %>% select(-closeness, -betweenness)
d2022 <- lapply(net.met2022_species, function(x) x %>% rownames_to_column()) %>% bind_rows(.id = 'Site') %>% select(-closeness, -betweenness)


rbind(d2021, d2022) %>% group_by(rowname) %>% summarise(d.m = mean(d), degree.m = mean(normalised.degree), cl.m = mean(weighted.closeness)) %>% filter(rowname == 'Apis mellifera'| rowname == 'Bombus lapidarius')

###################
### NULL MODELS ###
###################


# creating 1000 null models from vazull (connectance is constrained)
null.2021 <- lapply(adj.matrix.2021, function(x) vaznull(1000, x))
null.2022 <- lapply(adj.matrix.2022, function(x) vaznull(1000, x))

for (i in 1:length(null.2021)) {
  for (k in 1:length(null.2021[[i]])) {
    colnames(null.2021[[i]][[k]]) <- colnames(adj.matrix.2021[[i]])
    rownames(null.2021[[i]][[k]]) <- rownames(adj.matrix.2021[[i]])
  }
}

for (i in 1:length(null.2022)) {
  for (k in 1:length(null.2022[[i]])) {
    colnames(null.2022[[i]][[k]]) <- colnames(adj.matrix.2022[[i]])
    rownames(null.2022[[i]][[k]]) <- rownames(adj.matrix.2022[[i]])
  }
}



## looping calculation of metrics for all null models

metrics_null <- function(matrix) {
network_parameters_null <- list()
for (i in 1:length(matrix)) {
  network_parameters_null[[i]] <- matrix[[i]] %>% lapply(networklevel, index=c('connectance',"weighted NODF", 
                                                                             "niche overlap"), level = 'higher', 
                                                       weighted = TRUE) %>% bind_rows()
  }
names(network_parameters_null) <- names(matrix)
return(network_parameters_null)
}

met.null.2021 <- metrics_null(null.2021)
met.null.2022 <- metrics_null(null.2022)

metrics_null_species <- function(matrix) {
  network_parameters_null <- list()
  for (i in 1:length(matrix)) {
    network_parameters_null[[i]] <- matrix[[i]] %>% lapply(specieslevel, index=c("d", 
                                                                                 "normalised degree", 
                                                                                 "closeness", "betweenness"), level = 'higher') %>% 
      lapply(function(x) x %>% rownames_to_column(var = 'Species')) %>% bind_rows()
  }
  names(network_parameters_null) <- names(matrix)
  return(network_parameters_null)
}

met.null.sp.2021 <- metrics_null_species(null.2021)
met.null.sp.2022 <- metrics_null_species(null.2022)


### calculating z-scores, code modified from https://fukamilab.github.io/BIO202/09-B-networks.html

met.zscore <-  function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
}


zscore <- function(met, met.null, metric) {
net.nest.zscore <- list() 
for(i in 1:length(met)){
  net.nest.zscore[[i]] = met.zscore(met[[i]][metric], 
                                    as.vector(met.null[[i]][ ,metric])[[1]])
}
names(net.nest.zscore) <- names(met)
return(net.nest.zscore)
}

zscore.nestedness.2021 <- zscore(net.met2021, met.null.2021, metric = 'weighted NODF') %>% bind_rows(.id = 'Site') %>% mutate(Year = 2021)
zscore.niche.2021 <- zscore(net.met2021, met.null.2021, metric = 'niche.overlap.HL') %>% bind_rows(.id = 'Site') %>% mutate(Year = 2021)
zscore.nestedness.2022 <- zscore(net.met2022, met.null.2022, metric = 'weighted NODF') %>% bind_rows(.id = 'Site') %>% mutate(Year = 2022)
zscore.niche.2022 <- zscore(net.met2022, met.null.2022, metric = 'niche.overlap.HL') %>% bind_rows(.id = 'Site') %>% mutate(Year = 2022)

z.nested <- rbind(zscore.nestedness.2021, zscore.nestedness.2022)
z.niche <- rbind(zscore.niche.2021, zscore.niche.2022)

connectance2021 <- net.met2021 %>% bind_rows(.id = 'Site') %>% select(Site, connectance) %>% mutate(Year = 2021)
connectance2022 <- net.met2022 %>% bind_rows(.id = 'Site') %>% select(Site, connectance) %>% mutate(Year = 2022)
conn <- rbind(connectance2021, connectance2022)

network.metrics.both <- data.frame(Site = z.nested$Site, Nested.z = z.nested$`weighted NODF`, Niche.z = z.niche$niche.overlap.HL, Connectance = conn$connectance, Year = z.nested$Year)

save(network.metrics.both, file = 'Data/network_metrics_both_years.RData')

zscore_sp <- function(met.sp, met.sp.null) {
  net.clos.zscore <- list()
  net.bet.zscore <- list()
  met.sp.null2 <- list()
  for(i in 1:length(met.sp)){
    met.sp.null2[[i]] <- met.sp.null[[i]] %>% group_by(Species) %>% summarise(d.m = mean(d), d.sd = sd(d), closeness.m = mean(weighted.closeness), closeness.sd = sd(weighted.closeness),
                                                         degree.m = mean(normalised.degree), degree.sd = sd(normalised.degree), betweenness.m = mean(betweenness), betweenness.sd = sd(betweenness))
    net.clos.zscore[[i]] = (met.sp[[i]]['weighted.closeness'] - met.sp.null2[[i]]['closeness.m'])/met.sp.null2[[i]]['closeness.sd']
    net.bet.zscore[[i]] = (met.sp[[i]]['weighted.betweenness'] - met.sp.null2[[i]]['betweenness.m'])/met.sp.null2[[i]]['betweenness.sd']
    net.clos.zscore[[i]] = cbind(net.clos.zscore[[i]], net.bet.zscore[[i]])
  }
  names(net.clos.zscore) <- names(met.sp)
  return(net.clos.zscore)
}

closeness2021 <- zscore_sp(net.met2021_species, met.null.sp.2021) %>% lapply(function(x) x %>% rownames_to_column(var = 'Species')) %>% 
  bind_rows(.id = 'Site') %>% mutate(Year = 2021) %>% rename(Closeness.z = weighted.closeness, Betweenness.z = weighted.betweenness)
closeness2022 <- zscore_sp(net.met2022_species, met.null.sp.2022) %>% lapply(function(x) x %>% rownames_to_column(var = 'Species')) %>% 
  bind_rows(.id = 'Site') %>% mutate(Year = 2022) %>% rename(Closeness.z = weighted.closeness, Betweenness.z = weighted.betweenness)

closeness.both <- rbind(closeness2021, closeness2022)

save(closeness.both, file = 'Data/closeness.RData')


##### alternatively calculating 95% CI for null models and comparing to observed values

a = t.test(met.null.2021[[1]]$`weighted NODF`)[["conf.int"]]
net.met2021[[1]][3]

#### CALCULATING INDIVIDUAL NICHE OVERLAP
library(spaa)

morisita.2021 <- map(adj.matrix.2021, function(x) as.matrix(niche.overlap(x, method = 'morisita')))
morisita.2022 <- map(adj.matrix.2022, function(x) as.matrix(niche.overlap(x, method = 'morisita')))
#plotnetwork(morisita.2021[[1]])

## null
morisita.null.2021 <- map(null.2021, function(x) map(x, function(y) as.matrix(niche.overlap(y, method = 'morisita'))))
morisita.null.2022 <- map(null.2022, function(x) map(x, function(y) as.matrix(niche.overlap(y, method = 'morisita'))))


morisita_null_function <- function(morisita.null, morisita) {
  hm <- map(morisita.null, function(x) {
    Y <- do.call(cbind, x)
    Y <- array(Y, dim=c(dim(x[[1]]), length(x)))
    Y.m <- apply(Y, c(1, 2), mean, na.rm = TRUE)
    Y.sd <- apply(Y, c(1, 2), sd, na.rm = TRUE)
    Y2 <- list(mean = Y.m, sd = Y.sd)
  })
  zscore <- list()
  for (i in 1:length(hm)) {
    zscore[[i]] <- (morisita[[i]] - hm[[i]][['mean']])/hm[[i]][['sd']]
  }
  return(zscore)
}
morisita.zscore2021 <- morisita_null_function(morisita.null.2021, morisita.2021)
names(morisita.zscore2021) <- names(adj.matrix.2021)
morisita.zscore2022 <- morisita_null_function(morisita.null.2022, morisita.2022)
names(morisita.zscore2022) <- names(adj.matrix.2022)

#save(morisita.zscore2021, morisita.zscore2022, morisita.2021, morisita.2022, file = 'Data/niche_for_r0.RData')

morisita.zscore2021.2 <- map(morisita.zscore2021, function(x) {
 y <- as.data.frame(x) %>% filter(row.names(x) == 'Apis mellifera') %>% pivot_longer(cols = everything())
 hb <- mean(y$value, na.rm = T)
 y2 <- y %>% mutate(value = ifelse(name == 'Apis mellifera', hb, value)) %>% rename(Morisita.z = value, Species = name) %>%
   mutate(Species = sub(' agg.', '', Species))
}
) %>% bind_rows(.id = 'Site') %>% mutate(Year = 2021)

morisita.zscore2022.2 <- map(morisita.zscore2022, function(x) {
  y <- as.data.frame(x) %>% filter(row.names(x) == 'Apis mellifera') %>% pivot_longer(cols = everything())
  hb <- mean(y$value, na.rm = T)
  y2 <- y %>% mutate(value = ifelse(name == 'Apis mellifera', hb, value)) %>% rename(Morisita.z = value, Species = name) %>%
    mutate(Species = sub(' agg.', '', Species))
}
) %>% bind_rows(.id = 'Site') %>% mutate(Year = 2022)


morisita.zscore.both <- rbind(morisita.zscore2021.2, morisita.zscore2022.2)

save(morisita.zscore.both, morisita.zscore2021, morisita.zscore2022, morisita.zscore.bl.both, file = 'Data/morisita_zscore_both.RData')

morisita.zscore2021.bl <- map(morisita.zscore2021, function(x) {
  y <- as.data.frame(x) %>% filter(row.names(x) == 'Bombus lapidarius') %>% pivot_longer(cols = everything())
  hb <- mean(y$value, na.rm = T)
  y2 <- y %>% mutate(value = ifelse(name == 'Bombus lapidarius', hb, value)) %>% rename(Morisita.z = value, Species = name) %>%
    mutate(Species = sub(' agg.', '', Species))
}
) %>% bind_rows(.id = 'Site') %>% mutate(Year = 2021)

morisita.zscore2022.bl <- map(morisita.zscore2022, function(x) {
  y <- as.data.frame(x) %>% filter(row.names(x) == 'Bombus lapidarius') %>% pivot_longer(cols = everything())
  hb <- mean(y$value, na.rm = T)
  y2 <- y %>% mutate(value = ifelse(name == 'Bombus lapidarius', hb, value)) %>% rename(Morisita.z = value, Species = name) %>%
    mutate(Species = sub(' agg.', '', Species))
}
) %>% bind_rows(.id = 'Site') %>% mutate(Year = 2022)

morisita.zscore.bl.both <- rbind(morisita.zscore2021.bl, morisita.zscore2022.bl)
