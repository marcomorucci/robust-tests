
library(reshape2)
library(ggplot2)
library(ggjoy)
library(latex2exp)
library(plyr)

## Check if simulation data is already present -- runs simulations if not
## Warning: running simulations will take a long time.
if(!file.exists('method_comparison2.csv'))
  system('python method_comparison.py')
if(!file.exists('confounded_method_comparison3.csv'))
  system('python confounded_method_comparison.py')

## Plot results

## Simulation 1
sdta = read.csv('method_comparison2.csv')
sdta = melt(sdta[, c('N', 'ATE', 'pval_avg', 'max_p', 'min_p', 'rnd_p')], id.vars = c('N', 'ATE'))
sdta$hp = factor(sdta$ATE > 0, levels = c(FALSE, TRUE), labels = c('SATE = 0', 'SATE != 0'))
sdta$N = factor(sdta$N, labels= c('N = 20', 'N = 50', 'N = 100'))
sdta$variable = factor(sdta$variable, 
                       labels = c('Robust', 'Largest\n first', 'Smallest\n first', 'Random\n order'))
smeans = ddply(sdta, c('hp', 'N', 'variable'), function(x) mean(x$value))
ggplot(sdta, aes(x=value, y=variable)) +
  geom_density_ridges(color='black', size=0.2, fill=alpha('grey', 0.5)) + 
  geom_text(data=smeans, aes(x=V1, y=variable, label=paste('avg =', round(V1, 3))), nudge_y = 0.2) + 
  xlab('p-value') + ylab('') + 
  facet_grid(hp ~ N) + theme_bw() +
  theme(legend.position = 'none')

ggsave('simulation1.png', width = 10, height=7, dpi = 300)


#Simulation 2
sdta = read.csv('confounded_method_comparison3.csv')
sdta = melt(sdta[, c('dlt', 'ATE', 'pval_avg', 'max_p', 'min_p', 'rnd_p')], id.vars = c('dlt', 'ATE'))
sdta$hp = factor(sdta$ATE > 0, levels = c(FALSE, TRUE), labels = c('SATE = 0', 'SATE != 0'))
sdta$dlt = factor(sdta$dlt, labels= c('Light confounding: delta=0.1', 
                                      'Moderate confounding: delta=0.5', 
                                      'Strong confounding: delta=0.9'))
sdta$variable = factor(sdta$variable, 
                       labels = c('Robust', 'Largest\n first', 'Smallest\n first', 'Random\n order'))
smeans = ddply(sdta, c('hp', 'dlt', 'variable'), function(x) mean(x$value))
ggplot(sdta, aes(x=value, y=variable)) +
  geom_density_ridges(color='black', size=0.2, fill=alpha('grey', 0.5)) + 
  geom_text(data=smeans, aes(x=V1, y=variable, label=paste('avg =', round(V1, 3))), nudge_y = 0.2) + 
  xlab('p-value') + ylab('') + 
  facet_grid(hp ~ dlt) + theme_bw() +
  theme(legend.position = 'none')

ggsave('simulation2.png', width = 10, height=7, dpi = 300)
