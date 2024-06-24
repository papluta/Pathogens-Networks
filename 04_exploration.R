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
  mutate(across(Ann.fl:Edge.dens, ~ scale(.)[,1]))

data.bb <- data %>% filter(Group == 'bb' | Group == 'bp') %>% 
  mutate(across(Ann.fl:Edge.dens, ~ scale(.)[,1]))

data.wb <- data %>% filter(Group == 'wb') %>% 
  mutate(across(Ann.fl:Edge.dens, ~ scale(.)[,1])) %>%
  mutate(Genus = gsub(' .*', '', Species)) %>%
  # filtering out erroneous barcodes
  filter(Genus != '') %>% filter(Genus != ('Sipha')) %>% filter(Genus != ('Trypoxylon')) %>% filter(Genus != ('Orisarma')) %>% 
  filter(Genus != ('Lindenius')) %>% filter(Genus != ('Oedogonium')) %>% filter(Genus != ('Orasema')) %>% filter(Genus != ('Megalocoleus')) %>% filter(Genus != ('Bombus'))

b <- data.wb %>% group_by(Species) %>% summarise(S_n = n()) %>% filter(S_n == 1) # 236 Andrena, 55 Lasioglossum
print(data.wb %>% group_by(Site) %>% summarise(n = n()), n = Inf) # Nor 174 only 6 samples, Wm1316 7

data.wb <- data.wb %>% left_join(b, by = 'Species') %>% filter(is.na(S_n))

### AES ###
# separate models for each virus and bee group
# default (improper) priors

dwvb.h <- brm(bf(DWVB.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.b <- brm(bf(DWVB.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + Species + (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site) + (1|Species),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.h <- brm(bf(BQCV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.b <- brm(bf(BQCV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.h <- brm(bf(ABPV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.b <- brm(bf(ABPV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.b <- brm(bf(SBV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density + (1 |Site) + (1|Genus),
                  hu ~ (Org.farm.p + SNH.p + Ann.fl.p) * Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
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

dens_plot_mu(AES_dens_models[['dwvb.h']], effects = 'SNH.p:Density') # Nor578 is an outlier for org farming
dens_plot_hu(AES_dens_models[[3]], effects = 'Ann.fl.p:Density') # Nor578 is an outlier for org farming
dens_plot_mu(AES_dens_models[[1]], effects = 'SNH.p:Density')


library(bipartite)
ex <- data %>% group_by(Species) %>% summarise(n = n()) %>% filter(n >2)

data1 <- data[c(1:7, 13)] %>% group_by(Density, Species) %>% summarise(dwvb = sum(dwvb), bqcv = sum(bqcv), sbv = sum(sbv), abpv = sum(abpv)) %>%
  mutate(Density = factor(Density, levels = c('0', '1'))) %>% filter(Species %in% ex$Species)



net.data <- split(data1, data1$Density)
net.data[[1]] <- net.data[[1]][-1]
net.data[[2]] <- net.data[[2]][-1]
low <- net.data$`0`
high <- net.data$`1`
low <- low[-1]
high <- high[-1]
low2 <- as.matrix(low)
high2 <- as.matrix(high)
rownames(low2) <- net.data[[1]]$Species
rownames(high2) <- net.data[[2]]$Species
a <- betalinkr(webs2array(low2, high2), partitioning="commondenom")
a

low3 <- t(low2)
plotweb(low2)
#low <- low %>% mutate(across(everything(), function(x) ifelse(x > 0, 1, 0)))


### WORKS
col <- adjustcolor(c('#7f5454','#6a7f54', '#547f7f', '#6a547f'), alpha = 0.6)
low$Species <- net.data[[1]]$Species
g.data.low <- low %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% 
  mutate(e.col = rep(col, times = 17)) %>% filter(Yes != 0) %>% rename(weight = Yes) %>%
  #mutate(Species = sub('(.).* ', '. ', Species)) 
  mutate(Species = sub('Andrena', "A.", Species)) %>% mutate(Species = sub('Apis', "A.", Species)) %>% 
  mutate(Species = sub('Bombus', "B.", Species)) %>% 
  mutate(Species = sub('Lasioglossum', "L.", Species)) %>% 
  mutate(Species = sub('Melitta', "M.", Species)) %>% 
  mutate(Species = sub('Eucera', "E.", Species)) %>% mutate(Species = sub('Colletes', "C.", Species)) %>%
  mutate(Virus = recode_factor(as.factor(Virus), 'dwvb' = 'DWV-B', bqcv = 'BQCV', sbv = 'SBV', abpv = 'ABPV'))
#deg <- low %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% filter(Yes != 0) %>% select(Yes)
g.data.low <- g.data.low[,c(2,1,3,4)]
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
g <- graph.data.frame(g.data.low, directed = T)

V(g)$type <- bipartite.mapping(g)$type
V(g)$shape <- ifelse(V(g)$type == TRUE, "circle", "circle")
V(g)$color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54',  rep('#ecde51', times = 17))
V(g)$frame.color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54', rep('#ecde51', times = 17))
V(g)$size <- degree(g)+25
E(g)$color <- E(g)$e.col

E(g)$width <- E(g)$weight*0.1
plot(g, layout = layout_in_circle(g), edge.arrow.size = 0.5, vertex.label.color = 'black', vertex.label.cex = 0.8, vertex.frame.color = NULL)

## high dens
high$Species <- net.data[[2]]$Species
g.data.high <- high %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% 
  mutate(e.col = rep(col, times = 15)) %>% filter(Yes != 0) %>% rename(weight = Yes) %>%
  #mutate(Species = sub('(.).* ', '. ', Species)) 
  mutate(Species = sub('Andrena', "A.", Species)) %>% mutate(Species = sub('Apis', "A.", Species)) %>% 
  mutate(Species = sub('Bombus', "B.", Species)) %>% 
  mutate(Species = sub('Lasioglossum', "L.", Species)) %>% 
  mutate(Species = sub('Melitta', "M.", Species)) %>% mutate(Species = sub('Seladonia', "S.", Species)) %>% 
  mutate(Species = sub('Eucera', "E.", Species)) %>% mutate(Species = sub('Colletes', "C.", Species)) %>%
  mutate(Virus = recode_factor(as.factor(Virus), 'dwvb' = 'DWV-B', bqcv = 'BQCV', sbv = 'SBV', abpv = 'ABPV'))
#deg <- high %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% filter(Yes != 0) %>% select(Yes)
g.data.high <- g.data.high[,c(2,1,3,4)]
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yelhighgreen"), alpha=.6)
g <- graph.data.frame(g.data.high, directed = T)

V(g)$type <- bipartite.mapping(g)$type
#V(g)$shape <- ifelse(V(g)$type == TRUE, "circle", "circle")
V(g)$color <- ifelse(V(g)$type  == TRUE, "#ecde51", "#ec5a51")
V(g)$color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54',  rep('#ecde51', times = 15))
V(g)$frame.color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54', rep('#ecde51', times = 15))
V(g)$size <- degree(g)+25
E(g)$color <- E(g)$e.col

E(g)$width <- E(g)$weight*0.1
plot(g, layout = layout_in_circle(g), edge.arrow.size = 0.5, vertex.label.color = 'black', vertex.label.cex = 0.8, vertex.frame.color = NULL)


dwvb.h <- brm(bf(DWVB.abs ~ weighted.nested + Density + (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.b <- brm(bf(DWVB.abs ~ weighted.nested + Density + Species + (1 |Site),
                 hu ~ weighted.nested + Density + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ weighted.nested + Density + (1 |Site) + (1|Species),
                 hu ~ weighted.nested + Density +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.h <- brm(bf(BQCV.abs ~ weighted.nested + Density +  (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.b <- brm(bf(BQCV.abs ~ weighted.nested + Density +  (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ weighted.nested + Density + (1 |Site) + (1|Genus),
                 hu ~ weighted.nested + Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.h <- brm(bf(ABPV.abs ~ weighted.nested * Density +  (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.b <- brm(bf(ABPV.abs ~ weighted.nested + Density +  (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ weighted.nested + Density + (1 |Site) + (1|Genus),
                 hu ~ weighted.nested + Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ weighted.nested + Density +  (1 |Site),
                hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.b <- brm(bf(SBV.abs ~ weighted.nested + Density +  (1 |Site),
                hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ weighted.nested + Density + (1 |Site) + (1|Genus),
                hu ~ weighted.nested + Density +  (1 |Site) + (1|Genus)), family = hurdle_lognormal(), data = data.wb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

SH_dens_models <- list(dwvb.h, dwvb.b, dwvb.w, 
                        bqcv.h, bqcv.b, bqcv.w, 
                        abpv.h, abpv.b, abpv.w, 
                        sbv.h, sbv.b, sbv.w)

names(SH_dens_models) <- c('dwvb.h', 'dwvb.b', 'dwvb.w', 
                            'bqcv.h', 'bqcv.b', 'bqcv.w', 
                            'abpv.h', 'abpv.b', 'abpv.w', 
                            'sbv.h', 'sbv.b', 'sbv.w')

cs.all <- lapply(SH_dens_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

pp.all <- lapply(SH_dens_models, pp_hurdle) %>% bind_rows(.id = 'id')

dens_plot_mu(SH_dens_models[[12]], effects = 'SH:Density') 
dens_plot_hu(SH_dens_models[[11]], effects = 'SH') 

conditional_effects(SH_dens_models[[11]], dpar = 'hu')

nm <- posterior_density(SH_dens_models, 'b_Density1')
plot(nm)

plot_forest(SH_dens_models[[1]], 'SH:Density', dpar = 'hu')

plot_forest(SH_dens_models, effect = 'Density')

plot <- lapply(SH_dens_models, custom_summary) 
plot2 <- lapply(plot, function(x) filter(x, Predictor == 'SH')) %>% bind_rows(.id = 'mod') %>%
  ggplot(aes(x = Estimate, y = mod))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
  geom_linerange(aes(xmin = CI_low_90, xmax = CI_high_90), linewidth = 1, color = 'black')+
  #geom_linerange(aes(xmin = CI_l90, xmax = CI_h90), linewidth = 2.1, color = '#a40b0b')+
  geom_point(shape = 21, size = 4, color = 'black', fill = '#dc3a3a')+
  theme_bw(base_size = 16, base_line_size = 16/44)+
  #theme(axis.text.y = element_blank())+
  labs(x = NULL, y = NULL)
plot2

## net

library(lme4)
library(glmmTMB)
library(DHARMa)
library(car)

net <- data %>% distinct(Site, connectance, niche.overlap.HL, `weighted.nested`, Density, SH, flcv.m)

m <- lm(`weighted.nested` ~ Density + SH, net)
m <- lm(connectance ~ Density + flcv.m, net)
summary(m)
DHARMa::testDispersion(m)
plot(so <- simulateResiduals(m))
