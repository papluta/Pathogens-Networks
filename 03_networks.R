source('02_land_extrapolation.R')

library(bipartite)
library(purrr)

trans2021.2 <- trans2021 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') %>% filter(Bee_species != 'no_identificiation')
trans2022.2 <- trans2022 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee')

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


net.met2021 <- adj.matrix.2021 %>% lapply(networklevel, index=c("connectance"), level = 'higher', 
                                            weighted = TRUE)
net.met2022 <- adj.matrix.2022 %>% lapply(networklevel, index=c("connectance"), level = 'higher', 
                                          weighted = TRUE)

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
  network_parameters_null[[i]] <- matrix[[i]] %>% lapply(networklevel, index=c('connectance'), level = 'higher', 
                                                       weighted = TRUE) %>% bind_rows()
  }
names(network_parameters_null) <- names(matrix)
return(network_parameters_null)
}

met.null.2021 <- metrics_null(null.2021)
met.null.2022 <- metrics_null(null.2022)


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


connectance2021 <- net.met2021 %>% bind_rows(.id = 'Site') %>% select(Site, connectance) %>% mutate(Year = 2021)
connectance2022 <- net.met2022 %>% bind_rows(.id = 'Site') %>% select(Site, connectance) %>% mutate(Year = 2022)
conn <- rbind(connectance2021, connectance2022)

network.metrics.both <- data.frame(Site = conn$Site, Connectance = conn$connectance, Year = conn$Year)

#save(network.metrics.both, file = 'Data/network_metrics_both_years.RData')

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

#save(morisita.zscore.both, morisita.zscore2021, morisita.zscore2022, morisita.zscore.bl.both, file = 'Data/morisita_zscore_both.RData')

Flower_rich1 <- trans2021.2 %>% group_by(Site) %>% summarise(FL_rich = n_distinct(visited_flower_species), Year = 2021)
Flower_rich2 <- trans2022.2 %>% group_by(Site) %>% summarise(FL_rich = n_distinct(visited_flower_species), Year = 2022)

Flower_rich <- rbind(Flower_rich1, Flower_rich2)

bee_rich1 <- trans2021.2 %>% group_by(Site) %>% summarise(Bee_rich = n_distinct(Bee_species), Year = 2021)
bee_rich2 <- trans2022.2 %>% group_by(Site) %>% summarise(Bee_rich = n_distinct(Bee_species), Year = 2022)

bee_rich <- rbind(bee_rich1, bee_rich2)
