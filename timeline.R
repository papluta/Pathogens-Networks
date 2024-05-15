library(ggplot2)
library(readr)

fl.cv <- read_csv('flower_cover2022.csv')
hive.data <- read_csv('beekepers_data2022.csv')
path.data <- read_csv('pathogen_data2022.csv')

ggplot(fl.cv  %>% mutate(Date = as.Date(Date, '%m/%d/%Y')), aes(Date, 90))+
  geom_area(aes(fill = as.factor(Run)), alpha = 0.6)+
  geom_area(data = hive.data %>% filter(Round == 'June'), aes(Date, 90), fill = '#ebb93a', col = NA, alpha = 0.7)+
  geom_area(data = hive.data %>% filter(Month == 'July'), aes(Date, 90), fill = '#ebb93a', col = NA, alpha = 0.7)+
  geom_area(data = hive.data %>% filter(Month == 'August'), aes(Date, 90), fill = '#ebb93a', col = NA, alpha = 0.7)+
  geom_area(data = hive.data %>% filter(Month == 'October'), aes(Date, 90), fill = '#ebb93a', col = NA, alpha = 0.7)+
  geom_area(data = hive.data %>% filter(Month == 'September'), aes(Date, 90), fill = 'darkgreen', col = NA, alpha = 0.5)+
  geom_area(data = path.data %>% mutate(Date.p = as.Date(Date.p, '%Y-%m-%d')) %>% filter(Month == 'July'), aes(Date.p, 90), fill = 'darkred', col = NA, alpha = 0.7)+
  geom_area(data = path.data %>% mutate(Date.p = as.Date(Date.p, '%Y-%m-%d')) %>% filter(Month == 'August'), aes(Date.p, 90), fill = 'darkred', col = NA, alpha = 0.7)+
  scale_x_date(date_breaks = '2 weeks', date_labels = '%dth %b')+
  scale_y_continuous(expand = c(0,2))+
  scale_fill_manual(values = c('#8ac9b6', '#8abcc9', '#8a9dc9'))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.9, hjust = 0.1), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.grid = element_blank())+
  labs(fill = 'Flower cover \n estimation', y = NULL, x = NULL)

