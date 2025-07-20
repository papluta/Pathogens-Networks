library(grid)
library(ggplot2)
library(cowplot)

custom_summary <- function(mod) {
  sum <- summary(mod)[['fixed']] %>% select(Estimate, Est.Error, Bulk_ESS, Tail_ESS) %>% 
    cbind(ci(as_draws_array(mod), ci = 0.95, method = 'ETI') %>% 
            filter(grepl('b_', Parameter)) %>% filter(!grepl('prior', Parameter)) %>% left_join(p_direction(mod)[,1:2], by = 'Parameter') %>% 
            select(-c(CI, Parameter)) %>% rename(CI_low_95 = CI_low, CI_high_95 = CI_high)) %>% 
    mutate(across(everything(), ~ round(., 2))) %>%
    relocate(c(CI_low_95, CI_high_95), .after = Est.Error) %>% mutate(Response = mod[["formula"]][["resp"]], Predictor = rownames(.)) %>% 
    relocate(c(Response, Predictor), .before = Estimate)
  rownames(sum) <- NULL
  return(sum)
}


plot_hu <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'hu')[[1]]
  dat <- mod$data %>% select(x = effects, y = colnames(plot %>% select(contains('abs')))) %>% 
    mutate(y = ifelse(y > 0, 1,0))
  pd <- custom_summary(mod) %>% filter(grepl(effects, Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 1) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, 1 - estimate__))+
    geom_jitter(data = dat, aes(x = x, y = y), alpha = 0.2, shape = 21, height = 0.025) +
    geom_ribbon(aes(ymin = 1 - `lower__`, ymax = 1 - `upper__`), fill = '#A9D18E')+
    geom_line(col = 'black', linewidth = 1)+
    labs(x = effects, y = 'Probability of detection', subtitle = paste0("pd = ", pd))+
    theme_bw(base_size = 16)+
    #scale_fill_manual(values = c( '#0a81f1','#f1c40f','#00806A'))+
    #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
    theme(panel.grid = element_blank(), legend.position = 'none',
          #plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"), plot.subtitle = element_text(size = 10, vjust = 0.1),
          plot.title = element_text(size = 16, hjust = 0.5))#+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1))
  return(plot2)
}

plot_mu <- function(mod, effects) {
  plot <- conditional_effects(mod, effects = effects, dpar = 'mu')[[1]]
  dat <- mod$data %>% select(x = effects, y = colnames(plot %>% select(contains('abs')))) %>% 
    filter(y > 0)
  pd <- custom_summary(mod) %>% filter(grepl(effects, Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 1, 0)) %>%
    filter(hu == 0) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
    ggplot(aes(effect1__, estimate__))+
    geom_jitter(data = dat, aes(x = x, y = log(y)), alpha = 0.2, shape = 21, height = 0.025) +
    geom_line(col = 'black', linewidth = 1)+
    geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`), fill = '#00806A', alpha = 0.3)+
    labs(x = effects, y = 'Viral load', subtitle = paste0("pd = ", pd))+
    theme_bw(base_size = 16)+
    #scale_fill_manual(values = c( '#0a81f1','#f1c40f','#00806A'))+
    #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
    theme(panel.grid = element_blank(), legend.position = 'none',
          #plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"), plot.subtitle = element_text(size = 10, vjust = 0.1),
          plot.title = element_text(size = 16, hjust = 0.5))
  return(plot2)
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


pp_hurdle <- function(mod) {
  pred <- posterior_predict(mod)
  y1 <- mod[["data"]][[1]]
  ppc_dens_overlay(y = log1p(y1), 
                   yrep = log1p(pred[1:length(y1),]))
}

posterior_density <- function(mod_list, variable) {
  post <- lapply(mod_list, as_draws_df) 
  post2 <- lapply(post, function (x) select(x, variable))
  post3 <- post2 %>% bind_rows(.id = 'mod') %>% rename(post = 2)
  dens_plot <- post3 %>%
    ggplot(aes(x = post, col = mod, fill = mod))+
    geom_density(alpha = 0.1, linewidth = 1.2)+
    theme_bw(base_size = 20)+
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 1)+
    scale_color_manual(values = c('#f1c40f','#f1c40f','#f1c40f','#00806A','#00806A','#00806A', '#0a81f1','#0a81f1','#0a81f1', '#dc78c7', '#dc78c7', '#dc78c7'))+
    scale_fill_manual(values = c('#f1c40f','#f1c40f','#f1c40f','#00806A','#00806A','#00806A', '#0a81f1','#0a81f1','#0a81f1', '#dc78c7', '#dc78c7', '#dc78c7'))
  return(dens_plot)
}

plot_forest <- function(mod, effect, dpar) {
  plot <- conditional_effects(mod, effects = effects, dpar = dpar)[[1]]
  dat <- mod$data %>% select(x = gsub(':.*','', effects), y = colnames(plot %>% select(contains('abs'))), Density) %>% mutate(y = ifelse( y >0, 1, 0))
  pd <- custom_summary(mod) %>% filter(grepl(gsub(':.*','', effects), Predictor)) %>% filter(grepl(":", Predictor)) %>% mutate(hu = ifelse(grepl("hu",Predictor), 'hu', 'mu')) %>%
    filter(hu == dpar) %>% select(pd) %>% mutate(pd = round(pd, digits = 2))
  plot2 <- plot %>% 
  ggplot(x, aes(x = Estimate, y = Response))+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
    geom_linerange(aes(xmin = CI_low_90, xmax = CI_high_90), linewidth = 1, color = 'black')+
    #geom_linerange(aes(xmin = CI_l90, xmax = CI_h90), linewidth = 2.1, color = '#a40b0b')+
    geom_point(shape = 21, size = 4, color = 'black', fill = '#dc3a3a')+
    theme_bw(base_size = 16, base_line_size = 16/44)+
    theme(axis.text.y = element_blank())+
    labs(x = NULL, y = NULL)
  return(plot2)
}

plot_forest <- function(mod, effect) {
  plot <- lapply(mod, custom_summary) 
  plot2 <- lapply(plot, function(x) filter(x, Predictor == effect)) %>% bind_rows(.id = 'mod') %>%
    ggplot(aes(x = Estimate, y = mod))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
  geom_linerange(aes(xmin = CI_low_90, xmax = CI_high_90), linewidth = 1, color = 'black')+
  #geom_linerange(aes(xmin = CI_l90, xmax = CI_h90), linewidth = 2.1, color = '#a40b0b')+
  geom_point(shape = 21, size = 4, color = 'black', fill = '#dc3a3a')+
  theme_bw(base_size = 16, base_line_size = 16/44)+
  labs(x = NULL, y = NULL)
  return(plot2)
}

