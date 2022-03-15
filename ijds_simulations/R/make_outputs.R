library(reshape2)
library(ggplot2)
library(dplyr)

figure_folder = "/home/mm7936/robust_test_replication/simulations/results/figures_final"


################################################################################
## Z test performance comparison: linear DGP, increasing n
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/z_linear"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("Z_linear.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## Z test performance comparison: complex DGP, increasing n
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/z_complex"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}

res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Complex DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("Z_complex.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## Z test performance comparison: stratified DGP, increasing n
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/z_stratified"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}

res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Stratified DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("Z_stratified.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)


################################################################################
## Z test performance comparison: p increases linear dgp
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/z_linear_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('p', 'att', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(p, att) %>% 
    summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('p', 'att'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("p (number of covariates)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("Z_linear_largep.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## Z test performance comparison: p increases complex dgp
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/z_complex_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('p', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(p, tau) %>% 
    summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('p', 'tau'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("p (number of covariates)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Complex DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("Z_complex_largep.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)



################################################################################
## McNemar test performance comparison: simple DGP, n increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mcn_linear"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
} 
res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") +
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mcn_linear.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## McNemar test performance comparison: complex DGP, n increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mcn_complex"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
} 
res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") +
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Complex DGP") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mcn_complex.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## McNemar test performance comparison: stratified DGP, n increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mcn_stratified"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
} 
res <- data.frame(res)[, c('n', 'tau', 'robust', 'true', 'naive', 
                            'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(n, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") +
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Stratified DGP") +
     theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mcn_stratified.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)


################################################################################
## McN test performance comparison: p increases linear dgp
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mcn_linear_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('p', 'tau', 'robust', 'true', 'naive', 
                            'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(p, tau) %>% 
    summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('p', 'tau'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("p (number of covariates)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mcn_linear_largep.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## McN test performance comparison: p increases complex dgp
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mcn_complex_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('p', 'tau', 'robust', 'true', 'naive', 
                            'pscore', 'l2', 'optimal')]
resavg  <- res %>% group_by(p, tau) %>% 
    summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('p', 'tau'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("p (number of covariates)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Complex DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mcn_complex_largep.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)



################################################################################
## Approximation comparison: linear DGP, n increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/approx_linear"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)

resavg  <- res[, c('n', 'tau', 'mip', 'approximation_1', 'approximation_2')] %>% 
            group_by(n, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("approx_linear.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)


resavg  <- res[, c('n', 'tau', 'mip_time.elapsed', 
                    'a1_time.elapsed', 'a2_time.elapsed')] %>% 
            group_by(n, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})
names(resavg) <- c('n', 'tau', 'MIP', 'approximation_1', "approximation_2")
ggdata = melt(resavg, id.vars = c('n', 'tau'))
names(ggdata) = c('n', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=n, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("N (sample size)") + 
    ylab("Runtime (seconds)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("approx_linear_time.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## Approximation comparison: linear DGP, p increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/approx_linear_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)

resavg  <- res[, c('p', 'tau', 'mip', 'approximation_1', 'approximation_2')] %>% 
            group_by(p, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})

ggdata = melt(resavg, id.vars = c('p', 'tau'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + ylab("P (number of covariates)") + 
    ylab("P(Reject at 0.05 lv.)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("approx_linear_largep.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

resavg  <- res[, c('p', 'tau', 'mip_time.elapsed', 
                    'a1_time.elapsed', 'a2_time.elapsed')] %>% 
            group_by(p, tau) %>% 
            summarize_all(function(x){ mean(x < 0.05, na.rm=T)})
names(resavg) <- c('p', 'tau', 'MIP', 'approximation_1', "approximation_2")
ggdata = melt(resavg, id.vars = c('p', 'tau'))
names(ggdata) = c('p', 'att', 'Method', 'value')
ph0 = ggplot(ggdata, aes(x=p, y=value, color=Method, linetype=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) +  facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + xlab("P (number of covariates)") + 
    ylab("Runtime (seconds)") + ggtitle("True ATT (Cohen's d) -- Linear DGP") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("approx_linear_largep_time.png", ph0, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)


################################################################################
## MBs comparison, linear, n increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mb_linear"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c(1, 4, 6:11)]
names(res)[2] <- 'att'

resavg  <- res %>% group_by(n, att) %>% 
    summarize(rt_pval = mean(rt_pval , na.rm=T),
              rt_max_att = mean(rt_max_att, na.rm=T), 
              rt_min_att = mean(rt_min_att, na.rm=T),
              mb_pval = mean(mb_pval , na.rm=T), 
              mb_max_att = mean(mb_max_att, na.rm=T),
              mb_min_att=mean(mb_min_att, na.rm=T) )


ggdata = melt(resavg[, c('n', 'att', 'rt_pval', 'mb_pval')], id.vars = c('n', 'att'))
ggdata$Method = ifelse(grepl('rt', ggdata$variable), 'Robust Test', 'Matching Bounds')

ph0 = ggplot(ggdata, aes(x=n, y=value, linetype=Method, color=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) + facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + ylab("P-value") + 
    xlab('N (sample size)') + 
    ggtitle("True ATT") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mb_pval.png", ph0, width = 15, height = 5, units='in', 
        dpi=300, path = figure_folder)

ggdata = cbind(melt(resavg[, c('n', 'att', 'rt_max_att', 'rt_min_att', 
                            'mb_max_att', 'mb_min_att')],
                    id.vars = c('n', 'att')))
ggdata$Bound = ifelse(grepl('max', ggdata$variable), 'max', 'min')
ggdata$Method = ifelse(grepl('rt', ggdata$variable), 'Robust Test', 'Matching Bounds')

ph1 = ggplot(ggdata, aes(x=n, y=value, color=Bound, linetype=Method, pch=Method)) + 
geom_line(size=1) + geom_point(size = 3) + facet_grid(.~att) + 
  scale_color_brewer(palette='Dark2') + ylab("Estimated ATT") + 
  xlab('N (sample size)') + 
  ggtitle("True ATT") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mb_att.png", ph1, width = 15, height = 5, units='in', 
        dpi=300, path = figure_folder)



################################################################################
## MBs comparison, linear, p increases
################################################################################
resdir <- "/scratch/mm7936/linear_sim/data/mb_linear_largep"

res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c(2, 4, 6:11)]
names(res)[2] <- 'att'

resavg  <- res %>% group_by(p, att) %>% 
    summarize(rt_pval = mean(rt_pval , na.rm=T),
              rt_max_att = mean(rt_max_att, na.rm=T), 
              rt_min_att = mean(rt_min_att, na.rm=T),
              mb_pval = mean(mb_pval , na.rm=T), 
              mb_max_att = mean(mb_max_att, na.rm=T),
              mb_min_att=mean(mb_min_att, na.rm=T) )


ggdata = melt(resavg[, c('p', 'att', 'rt_pval', 'mb_pval')], id.vars = c('p', 'att'))
ggdata$Method = ifelse(grepl('rt', ggdata$variable), 'Robust Test', 'Matching Bounds')

ph0 = ggplot(ggdata, aes(x=p, y=value, linetype=Method, color=Method, pch=Method)) + 
    geom_line(size=1) + geom_point(size=3) + facet_grid(.~att) + 
    scale_color_brewer(palette='Dark2') + ylab("P-value") + 
    xlab('P (number of covariates)') + 
    ggtitle("True ATT") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mb_largep_pval.png", ph0, width = 15, height = 5, units='in', 
        dpi=300, path = figure_folder)

ggdata = cbind(melt(resavg[, c('p', 'att', 'rt_max_att', 'rt_min_att', 
                            'mb_max_att', 'mb_min_att')],
                    id.vars = c('p', 'att')))
ggdata$Bound = ifelse(grepl('max', ggdata$variable), 'max', 'min')
ggdata$Method = ifelse(grepl('rt', ggdata$variable), 'Robust Test', 'Matching Bounds')

ph1 = ggplot(ggdata, aes(x=p, y=value, color=Bound, linetype=Method, pch=Method)) + 
geom_line(size=1) + geom_point(size = 3) + facet_grid(.~att) + 
  scale_color_brewer(palette='Dark2') + ylab("Estimated ATT") + 
  xlab('P (number of covariates)') + 
  ggtitle("True ATT") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("mb_largep_att.png", ph1, width = 15, height = 5, units='in', 
        dpi=300, path = figure_folder)



################################################################################
## Constraint complexity study: n increases
################################################################################
resdir <- '/scratch/mm7936/linear_sim/data/constraints'
res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)
resavg  <- res %>% group_by(n, att) %>% 
    summarize(mom_time.elapsed = mean(mom_time.elapsed),
                mom_pval = mean( mom_pval < 0.05, na.rm=T), 
                cal_time.elapsed=mean(cal_time.elapsed), 
                cal_pval = mean(cal_pval < 0.05, na.rm=T), 
                qua_time.elapsed=mean(qua_time.elapsed),
                qua_pval=mean(qua_pval < 0.05, na.rm=T))

ggdata  <- cbind(melt(resavg[, c(1,2,4,6,8)], id.vars=c('n', 'att')),
                 time=melt(resavg[, c(1,2,3,5,7)], id.vars=c('n', 'att'))[, 4])
names(ggdata)  <-  c('n', 'att', 'constraint', 'pval', 'time')
ggdata$constraint <- factor(ggdata$constraint,
                            labels=c('Moment', 'Caliper', 'Quantile'))
p1  <- ggplot(ggdata, aes(x=n, color=constraint, pch=constraint, linetype=constraint)) + 
    geom_line(aes(y=pval), size=1) + geom_point(aes(y=pval), size=3) +
    scale_color_brewer(palette='Dark2') + ylab("P(reject at 0.05 lv.)") +
    xlab("N (sample size)") + 
    ggtitle("True ATT (Cohen's d) -- Linear DGP") + facet_grid(.~att) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("constraints_pval.png", p1, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

p2  <- ggplot(ggdata, aes(x=n, color=constraint, pch=constraint, linetype=constraint)) + 
    geom_line(aes(y=time), size=1) + geom_point(aes(y=time), size=3) +
    scale_color_brewer(palette='Dark2') + ylab("Runtime (seconds)") +
    xlab("N (sample size)") + 
    ggtitle("True ATT (Cohen's d) -- Linear DGP") + facet_grid(.~att) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("constraints_time.png", p2, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

################################################################################
## Constraint complexity study: p increases
################################################################################
resdir <- '/scratch/mm7936/linear_sim/data/constraint_p'
res <- NULL
for (fpath in list.files(resdir)){
    df <- readRDS(paste0(resdir, '/', fpath))
    if(ncol(df) != ncol(res) && !is.null(res)){
        next
    }
    res <- rbind(res, df)
}
res <- data.frame(res)[, c('p', 'tau', 'mom_time.elapsed', 'mom_pval', 'cal_time.elapsed', 
                            'cal_pval', 'qua_time.elapsed', 'qua_pval')]

resavg  <- res %>% group_by(p, tau) %>% 
    summarize(mom_time.elapsed = mean(mom_time.elapsed),
                mom_pval = mean( mom_pval < 0.05, na.rm=T), 
                cal_time.elapsed=mean(cal_time.elapsed), 
                cal_pval = mean(cal_pval < 0.05, na.rm=T), 
                qua_time.elapsed=mean(qua_time.elapsed),
                qua_pval=mean(qua_pval < 0.05, na.rm=T))

ggdata  <- cbind(melt(resavg[, c(1,2,4,6,8)], id.vars=c('p', 'tau')),
                 time=melt(resavg[, c(1,2,3,5,7)], id.vars=c('p', 'tau'))[, 4])
names(ggdata)  <-  c('p', 'att', 'constraint', 'pval', 'time')
ggdata$constraint <- factor(ggdata$constraint,
                            labels=c('Moment', 'Caliper', 'Quantile'))
p1  <- ggplot(ggdata, aes(x=p, color=constraint, pch=constraint, linetype=constraint)) + 
    geom_line(aes(y=pval), size=1) + geom_point(aes(y=pval), size=3) +
    scale_color_brewer(palette='Dark2') + ylab("P(reject at 0.05 lv.)") +
    xlab("p (number of covariates)") + 
    ggtitle("True ATT (Cohen's d) -- Linear DGP") + facet_grid(.~att) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("constraints_p_pval.png", p1, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)

p2  <- ggplot(ggdata, aes(x=p, color=constraint, pch=constraint, linetype=constraint)) + 
    geom_line(aes(y=time), size=1) + geom_point(aes(y=time), size=3) +
    scale_color_brewer(palette='Dark2') + ylab("Runtime (seconds)") +
    xlab("p (number of covariates)") +
    ggtitle("True ATT (Cohen's d) -- Linear DGP") + facet_grid(.~att) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggsave("constraints_p_time.png", p2, width = 15, height = 5, units='in',
       dpi=300, path = figure_folder)