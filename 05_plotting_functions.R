library(grid)
library(ggplot2)

custom_summary <- function(mod) {
  sum <- summary(mod)[['fixed']] %>% select(Estimate, Est.Error, Bulk_ESS, Tail_ESS) %>% 
    cbind(ci(as_draws_array(mod), ci = 0.9, method = 'ETI') %>% 
            filter(grepl('b_', Parameter)) %>% filter(!grepl('prior', Parameter)) %>% left_join(p_direction(mod)[,1:2], by = 'Parameter') %>% select(-c(CI, Parameter)) %>% rename(CI_low_90 = CI_low, CI_high_90 = CI_high)) %>% 
    mutate(across(everything(), ~ round(., 2))) %>%
    relocate(c(CI_low_90, CI_high_90), .after = Est.Error) %>% mutate(Response = mod[["formula"]][["resp"]], Predictor = rownames(.)) # %>% relocate(c(Response, Predictor), .before = Estimate)
  #rownames(sum) <- NULL
  return(sum)
}

custom_forest_data <- function(mod) {
  post <- as_draws_array(mod) 
  ci.95 <- ci(post, method = 'ETI') %>% filter(grepl(":", Parameter)) %>% rowwise() %>% mutate(mu_95 = (CI_low+CI_high)/2) %>%
    rename(CI_l95 = CI_low, CI_h95 = CI_high)
  ci.90 <- ci(post, method = 'ETI',ci = 0.90) %>% filter(grepl(":", Parameter)) %>% rowwise() %>% mutate(mu_90 = (CI_low+CI_high)/2)%>%
    rename(CI_l90 = CI_low, CI_h90 = CI_high)
  ci.50 <- ci(post, method = 'ETI', ci = 0.5) %>% filter(grepl(":", Parameter)) %>% rowwise() %>% mutate(mu_50 = (CI_low+CI_high)/2) %>%
    rename(CI_l50 = CI_low, CI_h50 = CI_high)
  ci <- ci.95 %>% left_join(ci.90, by = 'Parameter') %>% left_join(ci.50, by = 'Parameter')
  return(ci)
}

plot_forest <- function(x) {
  ggplot(x, aes(x = mu_95, y = Radius))+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
    geom_linerange(aes(xmin = CI_l95, xmax = CI_h95), linewidth = 1, color = '#464645')+
    geom_linerange(aes(xmin = CI_l90, xmax = CI_h90), linewidth = 2.1, color = '#a40b0b')+
    geom_point(shape = 21, size = 4, color = 'black', fill = '#dc3a3a')+
    theme_bw(base_size = 16, base_line_size = 16/44)+
    theme(axis.text.y = element_blank())+
    labs(x = NULL, y = NULL)
}


custom_fac_plot <- function(variable, mod, CI = 0.90, dpar = NULL) {
  plot <- conditional_effects(mod, effects = variable, prob = CI, dpar = dpar)[[1]] %>%
    ggplot(aes(x = .[,1], y = `estimate__`, fill = as.factor(Density),  # x = .[,1] if Density is behind the continuous variable and .[,2] if it's before
               color = as.factor(Density)))+
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.3, linewidth = 1)+
    geom_point(size = 3)+
    theme_bw(base_size = 22)+
    theme(legend.position = 'none')+
    scale_fill_manual(labels = c('0' = 'Low', '1' = 'High'), values = c('#f1c40f','#0a81f1'))+
    scale_color_manual(labels = c('0' = 'Low', '1' = 'High'), values = c('#f1c40f','#0a81f1'))+
    xlab(variable)
  return(plot)
}

## mu

dens_plot_mu <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'mu')[[1]]
  dat <- mod$data %>% select(x = gsub(':.*','', effects), y = colnames(plot %>% select(contains('abs'))), Density) %>% filter(y != 0)
  pd <- custom_summary(mod) %>% filter(grepl(gsub(':.*','', effects), Predictor)) %>% filter(grepl(":", Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 0) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, estimate__))+
    geom_point(data = dat, aes(x = x, y = log(y), col = Density),alpha = 0.2, shape = 19) +
    geom_line(aes(col = effect2__), linewidth = 1)+
    geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, fill = effect2__), alpha = 0.3)+
    annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                                               gp=gpar(col="black", fontsize=10, fontface=NULL))))+
    labs(x = NULL, y = NULL)+
    theme_bw(base_size = 16)+
    scale_fill_manual(values = c('#f1c40f', '#0a81f1'))+
    scale_color_manual(values = c('#f1c40f', '#0a81f1'))+
    theme(legend.position = 'none', panel.grid = element_blank(), 
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"), plot.title = element_text(size = 16, hjust = 0.5))
  return(plot2)
}

dens_plot_hu <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'hu')[[1]]
  dat <- mod$data %>% select(x = gsub(':.*','', effects), y = colnames(plot %>% select(contains('abs'))), Density) %>% mutate(y = ifelse( y >0, 1, 0))
  pd <- custom_summary(mod) %>% filter(grepl(gsub(':.*','', effects), Predictor)) %>% filter(grepl(":", Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 1) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, (1 - estimate__)))+
    geom_point(data = dat, aes(x = x, y = y, col = Density, fill = Density),alpha = 0.2, shape = 21) +
    scale_size_area()+
    geom_line(aes(col = effect2__), linewidth = 1)+
    geom_ribbon(aes(ymin = 1 - `lower__`, ymax = 1 - `upper__`, fill = effect2__), alpha = 0.3)+
    labs(x = NULL, y = NULL, subtitle = paste0("pd = ", pd))+
    theme_bw(base_size = 16)+
    scale_fill_manual(values = c('#f1c40f', '#0a81f1'))+
    scale_color_manual(values = c('#f1c40f', '#0a81f1'))+
    theme(legend.position = 'none', panel.grid = element_blank(), 
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"), plot.subtitle = element_text(size = 10, vjust = 0.1),
          plot.title = element_text(size = 16, hjust = 0.5))+
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1))
  return(plot2)
}

force_plot <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'mu')[[1]] %>% mutate(effect2__ = as.numeric(as.character(effect2__)))
  dat <- mod$data %>% select(x = gsub(':.*','', effects), y = colnames(plot %>% select(contains('abs'))), Flower = colnames(plot %>% select(contains('Iji')))) %>% 
    filter(y != 0) %>% mutate(Flower2 = as.factor(cut(Flower, 3, labels = c('1000', '2000', '3000'))))
  pd <- custom_summary(mod) %>% filter(grepl(":", Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 0) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, estimate__))+
    geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
    geom_line(aes(col = as.factor(effect2__)), linewidth = 1)+
    geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, fill = as.factor(effect2__)), alpha = 0.3)+
    annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                                               gp=gpar(col="black", fontsize=10, fontface=NULL))))+
    labs(x = NULL, y = NULL)+
    theme_bw(base_size = 16)+
    scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
    scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
    theme(panel.grid = element_blank(),  legend.position = 'none',
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
          plot.title = element_text(size = 16, hjust = 0.5))
  return(plot2)
}

force_plot_hu <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'hu')[[1]]
  dat <- mod$data %>% select(x = gsub(':.*','', effects), y = colnames(plot %>% select(contains('abs'))), Flower = colnames(plot %>% select(contains(gsub('.*:','', effects))))) %>% 
    mutate(y = ifelse(y > 0, 1,0)) %>% mutate(Flower2 = as.factor(cut(Flower, 3)))
  pd <- custom_summary(mod) %>% filter(grepl(":", Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 1) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, 1 - estimate__))+
    geom_point(data = dat, aes(x = x, y = y),alpha = 0.2, shape = 21) +
    geom_line(aes(col = effect2__), linewidth = 1)+
    geom_ribbon(aes(ymin = 1 - `lower__`, ymax = 1 - `upper__`, fill = effect2__), alpha = 0.3)+
    labs(x = NULL, y = NULL, subtitle = paste0("pd = ", pd))+
    theme_bw(base_size = 16)+
    scale_fill_manual(values = c( '#0a81f1','#f1c40f','#00806A'))+
    scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
    theme(panel.grid = element_blank(), legend.position = 'none',
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"), plot.subtitle = element_text(size = 10, vjust = 0.1),
          plot.title = element_text(size = 16, hjust = 0.5))+
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1))
  return(plot2)
}

pp_hurdle <- function(mod) {
  pred <- posterior_predict(mod)
  y1 <- mod[["data"]][[1]]
  ppc_dens_overlay(y = log1p(y1), 
                   yrep = log1p(pred[1:length(y1),]))
}
