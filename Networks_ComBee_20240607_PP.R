pol <- read.csv("Data/Transects2022_20231110.csv", stringsAsFactors = T)

str(pol)

levels(pol$Bee_group)

pol <- pol[pol$Bee_group=="Bombus"|pol$Bee_group=="Wildbee"|pol$Bee_group=="Honeybee",]

pol2 <- pol[pol$run==2,]

library(bipartite)
library(tidyverse)
library(dplyr)

adj.matrix <- frame2webs(pol2, varnames = c("visited_flower_species", "Bee_species",
                                            "Landscape_ID"), emptylist = TRUE)

adj.matrix

### not using this syntax because the sites needs to be in a list ###

#network_parameters <- adj.matrix |> map_df(~ networklevel(.x, index=c("connectance", 
#                                                                 "weighted NODF", 
#                                                                 "niche overlap") ,
#                                                     level="higher", 
#                                                  weighted = TRUE)) %>% mutate(Landscape_ID = unique(pol2$Landscape_ID))

network_parameters <- adj.matrix %>% lapply(networklevel, index=c("connectance", 
                                                                  "weighted NODF", 
                                                                  "niche overlap"), level = 'higher', 
                                            weighted = TRUE)

## NULL MODELS

# creating 1000 null models from vazull (connectance is constrained)
null <- lapply(adj.matrix, function(x) vaznull(1000, x))

## looping calculation of metrics for all null models
network_parameters_null <- list()
for (i in 1:length(null)) {
  network_parameters_null[[i]] <- null[[i]] %>% lapply(networklevel, index=c("connectance", ## connectance values are always the same!
                                                                             "weighted NODF", 
                                                                             "niche overlap"), level = 'higher', 
                                                       weighted = TRUE) %>% bind_rows()
  }
names(network_parameters_null) <- unique(pol2$Landscape_ID)

### calculating z-scores for NODF, code modified from https://fukamilab.github.io/BIO202/09-B-networks.html

net.zscore = function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
}


net.nest.zscore <- list() 
for(i in 1:length(network_parameters)){
  net.nest.zscore[[i]] = net.zscore(network_parameters[[i]]['weighted NODF'], 
                                    as.vector(network_parameters_null[[i]][ ,'weighted NODF'])[[1]])
}
names(net.nest.zscore) <- unique(pol2$Landscape_ID)

# THE Z-SCORES FOR WEIGHTED NODF
NODF.z <- net.nest.zscore %>% bind_rows(.id = 'Site')

#that's the raw NODF score for comparison:
m <- lapply(network_parameters, function(x) x[2]) %>% bind_rows(.id = 'Site')


### now niche overlap
net.nest.zscore <- list() 
for(i in 1:length(network_parameters)){
  net.nest.zscore[[i]] = net.zscore(network_parameters[[i]]['niche.overlap.HL'], 
                                    as.vector(network_parameters_null[[i]][ ,'niche.overlap.HL'])[[1]])
}
names(net.nest.zscore) <- unique(pol2$Landscape_ID)

niche.overlap.z <- net.nest.zscore %>% bind_rows(.id = 'Site')

connectance <- lapply(network_parameters, function(x) x[1]) %>% bind_rows(.id = 'Site')


#### CENTRALITY AT THE SPECIES LEVEL

# closeness centrality for honey bees
network_parameters_species <- adj.matrix %>% lapply(specieslevel, index = 'closeness', level = 'higher')

centrality_Apis <- lapply(network_parameters_species, function(x) x %>% rownames_to_column(var = 'Species') %>%
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

## how to deal with Nor264?
Nor264 <- (0 - null.mean)/null.sd
#or
Nor264 <- -3

## final data frame
networks.z <- data.frame(Site = unique(pol$Landscape_ID))
networks.z2 <- networks.z %>% left_join(connectance, by = 'Site') %>% left_join(NODF.z, by = 'Site') %>%
  left_join(niche.overlap.z, by = 'Site') %>% rename(Connectance = connectance, NODF.z = `weighted NODF`, Niche.z = niche.overlap.HL)

networks.z2$Site <- sub('Wm', 'WM', networks.z2$Site)



networks.z3 <- left_join(networks.z2, closeness.z, by = 'Site') %>% rename(Closeness.z = weighted.closeness)

networks.z3[23,5] <- Nor264

save(networks.z3, file = 'Data/networks_z.RData')
