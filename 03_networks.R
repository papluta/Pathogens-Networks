source('02_land_extrapolation_both_years.R')

library(bipartite)
library(purrr)

trans2021.2 <- trans2021 |> filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') |> filter(Bee_species != 'no_identificiation')
trans2022.2 <- trans2022 |> filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee')

#Chao2
Se = So + (L^2/2 * M)

chao2 = function(network) {
  network$int = paste0(network$Bee_species, network$visited_flower_species)
  So = nrow(network %>% distinct(int))
  L = nrow(network %>% group_by(int) %>% summarise(n = n()) %>% filter(n == 1))
  M = nrow(network %>% group_by(int) %>% summarise(n = n()) %>% filter(n == 2))
  Se = So + (L^2/2 * M)
  result = data.frame(Completeness = 100 * So / Se, Se = Se, So = So, L = L, M = M)
  return(result)
}

tran.split1 = split(trans2021.2, trans2021.2$Site)
tran.split2 = split(trans2022.2, trans2022.2$Site)
completness2021 = lapply(tran.split1, function(x) chao2(x)) %>% bind_rows(.id = 'Site')
completness2022 = lapply(tran.split2, function(x) chao2(x)) %>% bind_rows(.id = 'Site')




net = tran.split[[1]] %>% mutate(int = paste0(Bee_species, visited_flower_species))

So = nrow(net %>% distinct(int))
L = nrow(net %>% group_by(int) %>% summarise(n = n()) %>% filter(n == 1))
M = nrow(net %>% group_by(int) %>% summarise(n = n()) %>% filter(n == 2))


Se = So + (L^2/2 * M)



# n1 = trans2021.2 %>% group_by(Bee_species) %>% summarise(n = n())
# n2 = trans2022.2 %>% group_by(Bee_species) %>% summarise(n = n())
# mm = rbind(n1,n2) %>% group_by(Bee_species) %>% summarise(n2 = sum(n))
# 
# bee_rich2021 =  trans2021.2 %>% group_by(Site) %>% distinct(Bee_species) %>% summarise(Bee_rich = n()) %>% mutate(Year = 2021)
# bee_rich2022 =  trans2022.2 %>% group_by(Site) %>% distinct(Bee_species) %>% summarise(Bee_rich = n()) %>% mutate(Year = 2022)
# 
# bee_rich = rbind(bee_rich2021, bee_rich2022)
#save(bee_rich, file = 'Data/bee_richness.RData')
# Flower_rich2021 =  trans2021.2 %>% group_by(Site) %>% distinct(visited_flower_species) %>% summarise(FL_rich = n()) %>% mutate(Year = 2021)
# Flower_rich2022 =  trans2022.2 %>% group_by(Site) %>% distinct(visited_flower_species) %>% summarise(FL_rich = n()) %>% mutate(Year = 2022)
# 
# Flower_rich = rbind(Flower_rich2021, Flower_rich2022)
# save(Flower_rich, file = 'Data/flower_richness.RData')

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


#################################
####### NOT USED IN MS ##########
#################################

morisita_null_function_CI <- function(morisita.null) {
  hm <- map(morisita.null, function(x) {
    Y <- do.call(cbind, x)
    Y <- array(Y, dim=c(dim(x[[1]]), length(x)))
    Y.t <- apply(Y, c(1, 2), function(y) t.test(y)[["conf.int"]])
  })
  return(hm)
}

ci.morisita.2021 = morisita_null_function_CI(morisita.null = morisita.null.2021)
#ci.morisita.2022 = morisita_null_function_CI(morisita.null = morisita.null.2022)

morisita_CI_yesno <- function(ci.morisita.null, morisita) {
  mm <- list()
  for (i in 1:length(morisita)) {
    mm[[i]] = matrix(NA, nrow(morisita[[i]]), nrow(morisita[[i]]))
    for (k in 1:nrow(morisita[[i]])) {
      for (j in 1:nrow(morisita[[i]])) {
    mm[[i]][k,j] = ifelse(morisita[[i]][k,j] > ci.morisita.null[[i]][1,j,k] & morisita[[i]][k,j] < ci.morisita.null[[i]][2,j,k], 1, 0)
  }}}
  return(mm)
}

lol = morisita_CI_yesno(ci.morisita.2022,morisita.2022)


for (i in 1:length(lol)) {
  colnames(lol[[i]]) = colnames(morisita.2022[[i]])
  rownames(lol[[i]]) = rownames(morisita.2022[[i]])
}


res = data.frame(NA)
for (i in 1:length(lol)) {
  hm = as.data.frame(lol[[i]]) %>% filter(row.names(lol[[i]]) == 'Apis mellifera') %>% pivot_longer(cols = everything())
  res[i,1] = sum(hm$value, na.rm = T)
}

res = data.frame(NA)
for (i in 1:length(lol)) {
  hm = as.data.frame(lol[[i]]) %>% filter(row.names(lol[[i]]) == 'Bombus lapidarius') %>% pivot_longer(cols = everything())
  res[i,1] = sum(hm$value, na.rm = T)
}

res$Site = names(morisita.2021)


################



#D[jk] = (sum((prop[ij] + prop[ik]) * log(prop[ij] + prop[ik]) - sum(p[ij] * log(p[ij]) - sum(p[ik] * log(p[ik])))))/(2 * log(2))
#(sum((prop.resource[2,] + prop.resource[5,]) * log(prop.resource[2,] + prop.resource[5,]), na.rm = T) - sum(prop.resource[2,] * log(prop.resource[2,]), na.rm = T) - sum(prop.resource[5,] * log(prop.resource[5,]), na.rm = T))/(2 * log(2))

b <- Sys.time()
niche.list.func <- function(adj.matrix) {
  # Create a list to store niche overlap data for each adjacency matrix
  list_niche <- map(adj.matrix, function(a) {
    max.resource <- apply(a, 2, sum)
    prop_resource <- apply(a,1, function(x) x/max.resource)
    niche_overlap <- matrix(NA, nrow = nrow(prop_resource), ncol = nrow(prop_resource))
    # Calculate niche overlap for each pair of species
    for (i in seq_len(nrow(prop_resource))) {
      for (j in seq_len(nrow(prop_resource))) {
          combined <- prop_resource[i,] + prop_resource[j,]  
          niche_overlap[i, j] <- (sum(combined * log(combined), na.rm = TRUE) - 
                                    sum(prop_resource[i, ] * log(prop_resource[i, ]), na.rm = TRUE) - 
                                    sum(prop_resource[j, ] * log(prop_resource[j, ]), na.rm = TRUE)) / (2 * log(2))
        }
      }

    # Set column names for the niche overlap matrix
     colnames(niche_overlap) <- rownames(prop_resource)
     rownames(niche_overlap) <- rownames(prop_resource)
    return(niche_overlap)  # Return the processed data frame for the current matrix
  })
}

# Calculate niche pairs for two different years using the function
niche.pair.2021 <- niche.list.func(adj.matrix.2021)
niche.pair.2022 <- niche.list.func(adj.matrix.2022)



hb.overlap.2021 <- lapply(niche.pair.2021, function(x) x %>% filter(Bee_species == 'Apis mellifera') %>% arrange(Bee_species)) 
hb.overlap.2022 <- lapply(niche.pair.2022, function(x) x %>% filter(Bee_species == 'Apis mellifera') %>% arrange(Bee_species)) 

##### null #######
load('Data/list_null.RData')


niche.null.func <- function(null) {
  list_null <- list()
  prop.resource <- list()
  for (k in 1:length(null)) {
    list_null[[k]] <- list()
    prop.resource[[k]] <- list()
    for (g in 1:10) {
      a <- as.data.frame(null[[k]][[g]])
      max.resource <- apply(a, 2, sum)
      prop.resource[[k]][[g]] <- as.data.frame(apply(a,1, function(x) x/max.resource))
      list_null[[k]][[g]] <- list()
      for (i in 1:nrow(prop.resource[[k]][[g]])) {
        list_null[[k]][[g]][[i]] <- list() 
        for (j in 1:nrow(prop.resource[[k]][[g]])) {
          list_null[[k]][[g]][[i]][j] <- (sum((prop.resource[[k]][[g]][i,] + prop.resource[[k]][[g]][j,]) * log(prop.resource[[k]][[g]][i,] + prop.resource[[k]][[g]][j,]), na.rm = T) - sum(prop.resource[[k]][[g]][i,] * log(prop.resource[[k]][[g]][i,]), na.rm = T) - sum(prop.resource[[k]][[g]][j,] * log(prop.resource[[k]][[g]][j,]), na.rm = T))/(2 * log(2))
        }
      }
    }
  }
  names(list_null) <- names(null)
  return(list_null)
}

niche.null.2021 <- niche.null.func(null.2021)
niche.null.2022 <- niche.null.func(null.2022)

#save(niche.null.2021, file = 'Data/niche_null_2021.RData')


load('Data/niche_null_2021.RData')

naming <- function(list_null, adj.matrix) {
  for (i in 1:length(list_null)) {
    for (k in 1:length(list_null[[i]])) {
      names(list_null[[i]][[k]]) <- row.names(as.data.frame(t(adj.matrix[[i]])))
      for (j in 1:length(list_null[[i]][[k]])) {
        names(list_null[[i]][[k]][[j]]) <- names(list_null[[i]][[k]])
        list_null[[i]][[k]][[j]] <- apply(as.data.frame(list_null[[i]][[k]][[j]]), 2, function(x) round(x, digits = 3))
        list_null[[i]][[k]][[j]] <- as.data.frame(list_null[[i]][[k]][[j]]) %>% rename(Niche_overlap = 1) %>% rownames_to_column('Bee_species') %>% 
          mutate(Bee_species = sub("[.]", " ", Bee_species)) %>%
          mutate(Bee_species = sub("terrestris.agg", "terrestris agg", Bee_species)) %>%
          filter(Bee_species != names(list_null[[i]][[k]][j]))
      }
    }
  }
  return(list_null)
}
list.niche.null.2021 <- naming(niche.null.2021, adj.matrix.2021)

list.niche.null.2021 <- lapply(list_null, function(x) lapply(x, function(y) bind_rows(y, .id = 'Bee')))
list.niche.null2 <- lapply(list.niche.null, function(x) x %>% bind_rows(x, .id = 'iter'))
#list.niche.null2 <- do.call(rbind, list.niche.null)


niche.overlap.null <- list()
#nich.m <- list()
for (i in 1:length(list.niche.null2)) {
  niche.overlap.null[[i]] <- list.niche.null2[[i]] %>% group_by(Bee, Bee_species) %>% summarise(Niche = mean(Niche_overlap), Niche.sd = sd(Niche_overlap))
  niche.overlap.null[[i]] <- niche.overlap.null[[i]] %>% filter(Bee == 'Apis mellifera') %>% arrange(Bee_species)
}

#

names(niche.overlap.null) <- unique(pol2$Landscape_ID)

#hb.overlap.null <- lapply(list.niche.null2, function(x) x %>% filter(Bee == 'Apis mellifera'))
#hb.overlap2 <- hb.overlap %>% bind_rows(.id = 'Site') %>% rename(Species = Bee_species) %>% mutate(Site = sub('Wm','WM', Site))


niche.overlap.null <- niche.overlap.null[-23] # removing Nor264 because there was no HB in the network
hb.overlap <- hb.overlap[-23] # removing Nor264 because there was no HB in the network

net.nest.zscore <- list() 
for(i in 1:length(hb.overlap.null)){
  net.nest.zscore[[i]] <- list()
  for (k in 1:nrow(niche.overlap.null[[i]])) {
    net.nest.zscore[[i]][[k]] = (hb.overlap[[i]][k, 3] - niche.overlap.null[[i]][k, 3])/niche.overlap.null[[i]][k, 4]
  }
}

names(net.nest.zscore) <- unique(pol2$Landscape_ID)[-23]

for (i in 1:length(net.nest.zscore)) {
  names(net.nest.zscore[[i]]) <- niche.overlap.null[[i]]$Bee_species
}

niche.hb.zscore <- lapply(net.nest.zscore, function(x) bind_rows(x, .id = 'Species'))
niche.hb.zscore2 <- bind_rows(niche.hb.zscore, .id = 'Site') %>% rename(Niche.hb.z = Niche)
niche.hb.zscore.mean <- niche.hb.zscore2 %>% group_by(Species) %>% summarise(Niche.m.z = mean(Niche.hb.z))

save(niche.hb.zscore2, niche.hb.zscore.mean, file = 'Data/hb_overlap_zscore.RData')

######## CWM null

pol5[[1]] %>% group_by(Bee_species) %>% reframe(n_sum = n/sum(n), visited_flower_species)

pol5.null <- list()
hb <- list()
pol6.null <- list()
for (i in 1:length(null)) {
  pol5.null[[i]]<- list()
  hb[[i]] <- list()
  pol6.null[[i]] <- list()
  for (j in 1:length(null[[i]])) {
    pol5.null[[i]][[j]] <- as.data.frame(null[[i]][[j]]) %>% rownames_to_column() %>% rename(visited_flower_species = rowname) %>%pivot_longer(cols = contains(' '), names_to = 'Bee_species', values_to = 'n')
    hb[[i]][[j]] <- pol5.null[[i]][[j]] %>% filter(Bee_species == 'Apis mellifera') %>% rename(HB_n = n) %>% ungroup() %>% select(visited_flower_species, HB_n) 
    pol5.null[[i]][[j]] <- pol5.null[[i]][[j]] %>% left_join(hb[[i]][[j]], by = c('visited_flower_species')) %>% left_join(fl.ab3[[i]], by = join_by('visited_flower_species' == 'Flower_species')) %>%
      mutate(HB_n = ifelse(is.na(HB_n), 0, HB_n)) %>% drop_na(mean.ab) %>% filter(Bee_species != 'Apis mellifera') %>% 
      group_by(Bee_species) %>% reframe(prop_visited = n/sum(n), visited_flower_species, HB_n, mean.ab) %>%
      mutate(weigh.ab = prop_visited * HB_n * 1/mean.ab) ### INVERSE OF FLOWER COVER! SO THE LESSER THE ABUNDANCE THE MORE OVERLAP
    pol6.null[[i]][[j]] <- pol5.null[[i]][[j]] %>% group_by(Bee_species) %>% summarise(cwm.niche = mean(weigh.ab)) 
  }
  pol6.null[[i]] <- pol6.null[[i]] %>% bind_rows(.id = 'iter') %>% group_by(Bee_species) %>% summarise(mean.cwm.null = mean(cwm.niche), sd.cwm.null = sd(cwm.niche))
}
names(pol5.null) <- names(pol4)
names(pol6.null) <- names(pol4)

pol7.null <- pol6.null %>% bind_rows(.id = 'iter') %>% group_by(Bee_species) %>% summarise(mean.cwm.null = mean(cwm.niche), sd.cwm.null = sd(cwm.niche))

######### HURLBERT NICHE OVERLAP WITH HB (RESOURCE ABUDANCE)

pol5 <- list()
pol6 <- list()
hb <- list()
for (i in 1:length(pol4)) {
  pol5[[i]] <- pol4[[i]] %>% group_by(Bee_species, visited_flower_species) %>% summarise(n = n())
  pol5[[i]]$Bee_species <- as.character(pol5[[i]]$Bee_species)
  pol5[[i]]$visited_flower_species <- as.character(pol5[[i]]$visited_flower_species)
  hb[[i]] <- pol5[[i]] %>% filter(Bee_species == 'Apis mellifera') %>% 
    mutate(HB_n = n/sum(n), visited_flower_species) %>% ungroup() %>% select(visited_flower_species, HB_n) 
  pol5[[i]] <- pol5[[i]] %>% left_join(hb[[i]], by = c('visited_flower_species')) %>% cbind(fl.cv4[[i]]) %>%
    mutate(HB_n = ifelse(is.na(HB_n), 0, HB_n)) %>% filter(Bee_species != 'Apis mellifera') %>% 
    group_by(Bee_species) %>% reframe(prop_visited = n/sum(n), visited_flower_species, HB_n, mean.fl) %>%
    mutate(sim = prop_visited * HB_n * 1) ### INVERSE OF FLOWER COVER! SO THE LESSER THE ABUNDANCE THE MORE OVERLAP
  pol6[[i]] <- pol5[[i]] %>% group_by(Bee_species) %>% summarise(hurlb.niche = sum(sim)/mean(mean.fl)) 
  #pol6[[i]] <- split(pol5[[i]], pol5[[i]]$Bee_species)
}
names(pol5) <- names(pol4)
names(pol6) <- names(pol4)
hurlbert.niche.hb <- pol6 %>% bind_rows(.id = 'Site') %>% mutate(Site = sub('Wm', 'WM', Site), Bee_species = sub(' agg.', '', Bee_species)) 


### hurlbert null
pol5.null <- list()
hb <- list()
pol6.null <- list()
for (i in 1:length(null)) {
  pol5.null[[i]]<- list()
  hb[[i]] <- list()
  pol6.null[[i]] <- list()
  for (j in 1:length(null[[i]])) {
    pol5.null[[i]][[j]] <- as.data.frame(null[[i]][[j]]) %>% rownames_to_column() %>% rename(visited_flower_species = rowname) %>%
      pivot_longer(cols = contains(' '), names_to = 'Bee_species', values_to = 'n') %>% filter(n >0)
    hb[[i]][[j]] <- pol5.null[[i]][[j]] %>% filter(Bee_species == 'Apis mellifera') %>% 
      mutate(HB_n = n/sum(n), visited_flower_species) %>% ungroup() %>% select(visited_flower_species, HB_n) 
    pol5.null[[i]][[j]] <- pol5.null[[i]][[j]] %>% left_join(hb[[i]][[j]], by = c('visited_flower_species')) %>% 
      cbind(fl.cv4[[i]]) %>%
      mutate(HB_n = ifelse(is.na(HB_n), 0, HB_n)) %>% filter(Bee_species != 'Apis mellifera') %>% 
      group_by(Bee_species) %>% reframe(prop_visited = n/sum(n), visited_flower_species, HB_n, mean.fl) %>%
      mutate(sim = prop_visited * HB_n * 1) ### INVERSE OF FLOWER COVER! SO THE LESSER THE ABUNDANCE THE MORE OVERLAP
    pol6.null[[i]][[j]] <- pol5.null[[i]][[j]] %>% group_by(Bee_species) %>% summarise(hurlb.niche.z = mean(sim)/mean(mean.fl)) 
  }
  pol6.null[[i]] <- pol6.null[[i]] %>% bind_rows(.id = 'iter') %>% group_by(Bee_species) %>% summarise(hurlb.mean.z = mean(hurlb.niche.z), hurlb.sd.z = sd(hurlb.niche.z))
}
names(pol5.null) <- names(pol4)
names(pol6.null) <- names(pol4)

pol7.null <- pol6.null %>% bind_rows(.id = 'Site') %>% mutate(Site = sub('Wm','WM', Site), Bee_species = sub(' agg.', '', Bee_species))

hurlbert.niche.z <- hurlbert.niche.hb %>% left_join(pol7.null, by = join_by('Site','Bee_species')) %>% filter(Site != 'Nor264') %>%
  rowwise() %>% mutate(hurlb.z = (hurlb.niche - hurlb.mean.z)/hurlb.sd.z)

plot(hurlbert.niche.z2$hurlb.z ~ hurlbert.niche.z$hurlb.niche)
hurlbert.niche.hb <- hurlbert.niche.hb %>% filter(Site != 'Nor264')

save(hurlbert.niche.z, file = 'Data/hurlbert_zscore.RData')
#### CENTRALITY AT THE SPECIES LEVEL

# closeness centrality for honey bees
network_parameters_species <- adj.matrix %>% lapply(specieslevel, index = 'closeness', level = 'higher')
network_degree_species <- adj.matrix %>% lapply(specieslevel, index = 'normalised degree', level = 'higher')

centrality_Apis <- lapply(network_parameters_species, function(x) x %>% rownames_to_column(var = 'Species') %>%
                                        filter(Species == 'Apis mellifera'))
degree_Apis <- lapply(network_degree_species, function(x) x %>% rownames_to_column(var = 'Species') %>%
                            filter(Species == 'Apis mellifera'))
#centrality_BL <- lapply(network_parameters_species, function(x) x %>% rownames_to_column(var = 'Species') %>%
#                                    filter(Species == 'Bombus lapidarius'))
#centrality_BP <- lapply(network_parameters_species, function(x) x %>% rownames_to_column(var = 'Species') %>%
#                                    filter(Species == 'Bombus pascuorum'))


for(j in 1:length(null)){
  for(i in 1:length(null[[j]]))
    colnames(null[[j]][[i]]) <- rownames(network_parameters_species[[j]])
}


network_closeness_null <- list()
for (i in 1:length(null)) {
  network_closeness_null[[i]] <- null[[i]] %>% lapply(specieslevel, index=c('closeness'), level = 'higher')
}
names(network_closeness_null) <- unique(pol2$Landscape_ID)


Apis_null <- lapply(network_closeness_null, function(x) lapply(x, function(y) rownames_to_column(y, var = 'Species') %>%
         filter(Species == 'Apis mellifera')) %>% bind_rows())

closeness.zscore <- list() 
for(i in 1:length(centrality_Apis)){
  closeness.zscore[[i]] = net.zscore(centrality_Apis[[i]]['weighted.closeness'], 
                                     as.vector(Apis_null[[i]][ ,'weighted.closeness']))
}
names(closeness.zscore) <- unique(pol2$Landscape_ID)

null.av <- Apis_null %>% bind_rows()
null.mean <- mean(null.av$weighted.closeness)
null.sd <- sd(null.av$weighted.closeness)

closeness.z <- closeness.zscore %>% bind_rows(.id = 'Site') %>% mutate(Site = sub('Wm', 'WM', Site)) 


#### Apis degree

network_degree_null <- list()
for (i in 1:length(null)) {
  network_degree_null[[i]] <- null[[i]] %>% lapply(specieslevel, index=c('normalised degree'), level = 'higher')
}
names(network_degree_null) <- unique(pol2$Landscape_ID)


Apis_degree_null <- lapply(network_degree_null, function(x) lapply(x, function(y) rownames_to_column(y, var = 'Species') %>%
                                                                 filter(Species == 'Apis mellifera')) %>% bind_rows())

degree.zscore <- list() 
for(i in 1:length(degree_Apis)){
  degree.zscore[[i]] = net.zscore(degree_Apis[[i]]['normalised.degree'], 
                                     as.vector(Apis_degree_null[[i]][ ,'normalised.degree']))
}
names(degree.zscore) <- unique(pol2$Landscape_ID)

degree.z <- degree.zscore %>% bind_rows(.id = 'Site') %>% mutate(Site = sub('Wm', 'WM', Site)) %>% rename(Degree.z = normalised.degree)

## how to deal with Nor264 (no Apis)?
Nor264 <- (0 - null.mean)/null.sd
#or
Nor264 <- -3

## final data frame
networks.z <- data.frame(Site = unique(pol$Landscape_ID))
networks.z2 <- networks.z %>% left_join(connectance, by = 'Site') %>% left_join(NODF.z, by = 'Site') %>%
  left_join(niche.overlap.z, by = 'Site') %>% rename(Connectance = connectance, NODF.z = `weighted NODF`, Niche.z = niche.overlap.HL)

networks.z2$Site <- sub('Wm', 'WM', networks.z2$Site)

null.av <- Apis_degree_null %>% bind_rows()
null.mean <- mean(null.av$normalised.degree)
null.sd <- sd(null.av$normalised.degree)

networks.z3 <- left_join(networks.z2, closeness.z, by = 'Site') %>% rename(Closeness.z = weighted.closeness)
networks.z3 <- left_join(networks.z3, degree.z, by = 'Site')
networks.z3[23,5] <- Nor264
networks.z3[23,6] <- Nor264

save(networks.z3, file = 'Data/networks_z.RData')


bl2021 <- lapply(morisita.2021, function(x) as.data.frame(as.matrix(x)) %>% rownames_to_column() %>% 
                   select(rowname, `Bombus lapidarius`)) %>% bind_rows() %>% filter(rowname != 'Bombus lapidarius') %>%
                   summarise(bl.m = mean(`Bombus lapidarius`))
bl2022 <- lapply(morisita.2022, function(x) as.data.frame(as.matrix(x)) %>% 
                   filter(row.names(x) == 'Bombus lapidarius') %>% pivot_longer(cols = everything())) %>% bind_rows() %>% filter(name != 'Bombus lapidarius') %>%
                   summarise(bl.m = mean(`value`))

hb2021 <- lapply(morisita.2021, function(x) as.data.frame(as.matrix(x)) %>% 
                   filter(row.names(x) == 'Apis mellifera') %>% pivot_longer(cols = everything())) %>% bind_rows() %>% filter(name != 'Apis mellifera') %>%
  summarise(bl.m = mean(`value`))
hb2022 <- lapply(morisita.2022, function(x) as.data.frame(as.matrix(x)) %>% 
                   filter(row.names(x) == 'Apis mellifera') %>% pivot_longer(cols = everything())) %>% bind_rows() %>% filter(name != 'Apis mellifera') %>%
  summarise(bl.m = mean(`value`))

rbind(bl2021, bl2022) %>% summarise(m = mean(bl.m))              
rbind(hb2021, hb2022) %>% summarise(m = mean(bl.m))              


library(circlize)
library(wesanderson)

b = data.frame(col = wes_palette("Moonrise2"), group = c('hb', 'bb', 'wb', NA))
p = data.frame(p_col = c(wes_palette("Darjeeling2"), wes_palette("Royal2")), Plants = unique(dat$Plants))
as.data.frame(a)

dat = as.data.frame(adj.matrix.2022[[1]]) %>% rownames_to_column('Plants') %>% 
  relocate(Plants, .before = everything()) %>% pivot_longer(cols = 2:12, names_to = 'Bees', values_to = 'N') %>%
  mutate(group = ifelse(grepl('Apis', Bees), "hb", ifelse(grepl('Bombus', Bees), "bb", "wb"))) %>%
  left_join(b, by = 'group') %>% left_join(p, by = 'Plants') %>% mutate(Plants = sub(' ', '\n', Plants), Bees = sub(' ', '\n', Bees))


pdf('Data/fig/example_network.pdf', height = 20, width = 20)

par(cex = 3, mar = c(0, 0, 0, 0))
chordDiagram(dat, annotationTrack = c("grid", "name"), preAllocateTracks = list(track.height = 0.25), directional = -1,
             col = dat$p_col, grid.col = setNames(dat$p_col, dat$Plants), )
dev.off()

ggsave(a,file = 'Data/fig/network_ex.png', height = 20, width = 20)

