source('03_data_processing.R')
source('05_plotting_functions.R')

# choose scale (500/1000m)

data <- data_500

data <- data_1000

# correlation between extrapolated honeybees and flower cover standardized to agricultural area
plot(data$FL_per_agr ~ data$HB_per_agr)
plot(data$Connectance ~ data$Niche.z)
plot(data$Connectance ~ data$niche)
plot(data$Connectance ~ data$NODF.z)
plot(data$Connectance ~ data$Closeness.z)
plot(data$Niche.z ~ data$Closeness.z)
plot(data$NODF.z ~ data$Niche.z)
dotchart(data$HB_per_agr)
dotchart(data$FL_per_agr)
dotchart(data$sum.fl)


## separate data frames for each bee group
data.hb <- data %>% filter(Group == 'hb') %>% ungroup() %>% 
  mutate(across(sum.fl:sbv.f, ~ scale(.)[,1]))

data.bb <- data %>% filter(Group == 'bb' | Group == 'bp') %>% ungroup() %>% #filter(Species == 'Bombus lapidarius') %>%
  mutate(across(sum.fl:sbv.f, ~ scale(.)[,1]))


data.wb <- data %>% filter(Group == 'wb') %>% 
  mutate(Genus = gsub(' .*', '', Species)) %>%
  # filtering out erroneous barcodes
  filter(Genus != '') %>% filter(Genus != ('Sipha')) %>% filter(Genus != ('Trypoxylon')) %>% filter(Genus != ('Orisarma')) %>% 
  filter(Genus != ('Lindenius')) %>% filter(Genus != ('Oedogonium')) %>% filter(Genus != ('Orasema')) %>% filter(Genus != ('Megalocoleus')) %>% 
  filter(Genus != ('Bombus')) %>% ungroup() %>%
  mutate(across(sum.fl:sbv.f, ~ scale(.)[,1])) 

b <- data.wb %>% group_by(Species) %>% summarise(S_n = n()) %>% filter(S_n == 1) # finding species with only 1 datapoint
print(data.wb %>% group_by(Site) %>% summarise(n = n()), n = Inf) # Nor 174 only 6 samples, Wm1316 7

data.wb <- data.wb %>% left_join(b, by = 'Species') %>% filter(is.na(S_n)) %>% select(-S_n) # filtering out singleton species

data.site <- data %>% filter(Species == 'Apis mellifera') %>% distinct(Site, sum.hb, HB_per_agr, sum.fl, FL_per_agr, Density, 
                                                                       dwvb.hb.agr, bqcv.hb.agr, abpv.hb.agr, sbv.hb.agr, 
                                                                       dwvb.f, bqcv.f, abpv.f, sbv.f,
                                                                       Connectance, NODF.z, Niche.z, Closeness.z, Degree.z) %>%
  ungroup() %>% mutate(across(4:5, ~scale(.)[,1]))



### bayesian
library(brms)
library(tidybayes)
library(tidyr)
library(ggplot2)
library(rstan)
library(bayesplot)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(11)

load('Data/Models/240823NET_FORCE_models_500.RData')
load('Data/Models/240823NET_FORCE_models_1000.RData')

hbab <- brm(HB_per_agr ~ Density + FL_per_agr, family = 'Gaussian', prior = prior(normal(0,5), class = 'b'), data.site)
pp_check(hbab, ndraws = 100)
summary(hbab)
p_direction(hbab)
get_prior(Connectance ~ Density + FL_per_agr, family = 'Gaussian', prior = prior(normal(0,5), class = 'b'), data.site)

conn <- brm(Connectance ~ FL_per_agr + HB_per_agr, family = student(), prior = prior(normal(0,3), class = 'b'), data.site)
nich <- brm(Niche.z ~ FL_per_agr + HB_per_agr, family = 'Gaussian', prior = prior(normal(0,5), class = 'b'), data.site)
summary(conn)
summary(nich)
pairs(conn, variable = c('b_FL_per_agr', 'b_HB_per_agr'))
p_direction(nich)
p_direction(conn)
pp_check(conn, ndraws = 100)
pp_check(nich, ndraws = 100)

site_level_mods_1000 <- list(conn = conn, nich = nich, hbab = hbab)
site_level_mods_500 <- list(conn = conn, nich = nich, hbab = hbab)
cs.all <- lapply(site_level_mods_1000, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

closeness <- brm(Closeness.z ~ HB_per_agr + FL_per_agr, family = 'Gaussian', prior = prior(normal(0,5), class = 'b'), data.site)
pp_check(closeness, ndraws = 100)
summary(closeness)
p_direction(closeness)

deg <- brm(Degree.z ~ HB_per_agr + FL_per_agr, family = 'Gaussian', prior = prior(normal(0,5), class = 'b'), data.site)
pp_check(deg, ndraws = 100)
summary(deg)
p_direction(deg)

dwvb.h <- brm(bf(DWVB.abs ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site),
                 hu ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.h <- brm(bf(BQCV.abs ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site),
                 hu ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.h <- brm(bf(ABPV.abs ~ HB_per_agr + Connectance + Niche.z + Density +  (1 |Site),
                 hu ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site),
                hu ~ HB_per_agr + Connectance + Niche.z + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

NET_HB_models_1000 <- list(dwvb.h = dwvb.h, 
                      bqcv.h = bqcv.h, 
                      abpv.h = abpv.h, 
                      sbv.h = sbv.h)
NET_HB_models_500 <- list(dwvb.h = dwvb.h, 
                           bqcv.h = bqcv.h, 
                           abpv.h = abpv.h, 
                           sbv.h = sbv.h)

dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f + Connectance + Niche.z + Species + (1 |Site),
                  hu ~ dwvb.f + Connectance + Niche.z + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f + Connectance + Niche.z + (1 |Site) + (1|Species),
                  hu ~ dwvb.f + Connectance + Niche.z +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f + Connectance + Niche.z +  Species + (1 |Site),
                  hu ~ bqcv.f + Connectance + Niche.z + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f + Connectance + Niche.z + (1 |Site) + (1|Species),
                  hu ~ bqcv.f + Connectance + Niche.z +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.b <- brm(bf(ABPV.abs ~ abpv.f + Connectance + Niche.z +  Species +(1 |Site),
                  hu ~ abpv.f + Connectance + Niche.z + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ abpv.f + Connectance + Niche.z + (1 |Site) + (1|Species),
                  hu ~ abpv.f + Connectance + Niche.z +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.b <- brm(bf(SBV.abs ~ sbv.f + Connectance + Niche.z +  Species + (1 |Site),
                  hu ~ sbv.f + Connectance + Niche.z + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ sbv.f + Connectance + Niche.z + (1 |Site) + (1|Species),
                  hu ~ sbv.f + Connectance + Niche.z +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

NET_WB_models_1000 <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                            bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                            bqcv.w = bqcv.w, abpv.w = abpv.w, 
                            sbv.b = sbv.b, sbv.w = sbv.w)
NET_WB_models_500 <- list(dwvb.b = dwvb.b, dwvb.w = dwvb.w, 
                           bqcv.b = bqcv.b, bqcv.w = bqcv.w, 
                           bqcv.w = bqcv.w, abpv.w = abpv.w, 
                           sbv.b = sbv.b, sbv.w = sbv.w)

save(NET_WB_models_1000, site_level_mods_1000, NET_HB_models_1000, file = 'Data/Models/240823NET_FORCE_models_1000.RData')
save(NET_WB_models_500, site_level_mods_500, NET_HB_models_500, file = 'Data/Models/240823NET_FORCE_models_500.RData')

## plotting and diagnosing

cs.all <- lapply(NET_WB_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

cs.all <- lapply(NET_HB_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL
write.csv(cs.all, file = 'Data/Models/honey_bees_500.csv')

# to see which relationships were likely:
sum <- cs.all %>% filter(pd > 0.9) %>% filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept' & Predictor != 'SpeciesBombuspascuorum' & 
                                         Predictor != 'hu_SpeciesBombuspascuorum') %>% arrange(Predictor)

sum

write.csv(sum, file = 'Data/Models/sum_wb_1000.csv', row.names = F)

# POSTERIOR CHECK
pp.all <- lapply(NET_WB_models_500, pp_hurdle) #check for each model
pp.all[['abpv.b']]
pp.all[['sbv.b']]


pairs(dwvb.b, variable = c('b_dwvb.f', 'b_Connectance', 'b_Niche.z', 'b_SpeciesBombuspascuorum'))
pairs(dwvb.b, variable = c('b_hu_dwvb.f', 'b_hu_Connectance', 'b_hu_Niche.z', 'b_hu_SpeciesBombuspascuorum'))

## PLOTTING EFFECTS
conditional_effects(NET_WB_models_500[['bqcv.w']], effects = 'Connectance', dpar = 'mu')

mu_force <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'mu')[[1]])
mu_conn <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'mu')[[2]])
mu_niche <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'mu')[[3]])
hu_force <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'hu')[[1]])
hu_conn <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'hu')[[2]])
hu_niche <- lapply(NET_WB_models_500, function(x) conditional_effects(x, dpar = 'hu')[[3]])

conditional_effects(abpv.w, dpar = 'mu')


plots <- lapply(mu_force, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee force', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(mu_conn, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Connectance', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(mu_niche, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Niche overlap', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_force, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee force', y = 'Probability of 0')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_conn, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Connectance', y = 'Probability of 0')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_niche, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Niche overlap', y = 'Probability of 0')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plot(density(data.site$HB_per_agr))



## PLOTTING EFFECTS


######################### OTHER STUFF, NOT SORTED OUT ############################


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


