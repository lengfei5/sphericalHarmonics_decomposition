rm(list =ls())

require(ggplot2)
library(plotly)
library(ggrepel)
library(tidyverse)
library(patchwork)

##########################################
# import the output of SH 
##########################################
dataDir = "res/test_organoid_WTd6_v2"
Dirlist = list.dirs(path = dataDir, full.names = TRUE)

Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0
for(n in 1:length(Dirlist))
{
  if(file.exists(paste0(Dirlist[n], '/Condion_eachOrgnoid_SHspectrum_power.csv'))){
    i = i + 1
    xx = read.csv(paste0(Dirlist[n], '/Condion_eachOrgnoid_SHspectrum_power.csv'))
    cat(' nb cysts --- ', nrow(xx), '\n')
    
    if(i == 1){
      res6 = xx 
    }else{
      res6 = rbind(res6, xx)
    }
  }
}


dataDir = "res/test_organoid_WTd4_v2"
Dirlist = list.dirs(path = dataDir, full.names = TRUE)

Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0
for(n in 1:length(Dirlist))
{
  if(file.exists(paste0(Dirlist[n], '/Condion_eachOrgnoid_SHspectrum_power.csv'))){
    i = i + 1
    xx = read.csv(paste0(Dirlist[n], '/Condion_eachOrgnoid_SHspectrum_power.csv'))
    cat(' nb cysts --- ', nrow(xx), '\n')
    
    if(i == 1){
      res4 = xx 
    }else{
      res4 = rbind(res4, xx)
    }
  }
}

cat(' nb cysts for WT6 --- ', nrow(res6), '\n')
cat(' nb cysts for WT4 --- ', nrow(res4), '\n')

res6$condition = "wt_d6" 
res4$condition = "wt_d4" 
res = rbind(res6, res4)

res = res[, -1]

rm(res4); 
rm(res6)

##
## define which 
spect = as.matrix(res[, grep('power_per_l_', colnames(res))[1:21]])
l_max = apply(spect, 1, which.max) - 1
lmax_pct = apply(spect, 1,  function(x){x[which.max(x)]/sum(x)})

res$lmax = l_max
res$lmax_pct = lmax_pct
res$r2 = res$cyst_r^2
res$lmax_norm = res$lmax/res$r2


res = data.frame(res, stringsAsFactors = TRUE)

## filter the cyst sizes: tiny or segmentation artifacts
hist(log10(res$cyst_size), breaks = 50, xlim = c(2.0, 6.0))
abline(v = c(3, 4.7), col = 'red')

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^4.7), ]

hist(res$lmax)
table(res$lmax)
res = res[which(res$lmax > 0 & res$lmax <5), ]

ggplot(res, aes(x=r2, y=lmax, color=condition)) +
  geom_point()


ggplot(res, aes(x=lmax, y=lmax_pct, color=condition)) +
  geom_point()

ggplot(res, aes(x=r2, y=lmax_pct, color=condition)) +
  geom_point()


ggplot(res, aes(x=condition, y=r2, fill=condition)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) +
  #geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_vline(xintercept=5, col='darkgray') +
  #geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Condition", y= 'Cyst sizes (R2)') 

stat = table(res$lmax, res$condition)
stat[,1] = stat[,1]/sum(stat[,1])
stat[,2] = stat[,2]/sum(stat[,2])

ggplot(res, aes(x = lmax, fill = condition)) + 
  geom_histogram(position = "dodge")


ggplot(res, aes(x=condition, y=lmax, fill=condition)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) +
  #geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_vline(xintercept=5, col='darkgray') +
  #geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Condition", y= 'Cyst dominant SH component') 


ggplot(res, aes(x=condition, y=lmax_pct, fill=condition)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) +
  #geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_vline(xintercept=5, col='darkgray') +
  #geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Condition", y= 'weights of the dominant SH component') 


ggplot(res, aes(x=condition, y=lmax_pct, fill=condition)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size=10), 
      legend.title = element_text(size = 10)) +
  #geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_vline(xintercept=5, col='darkgray') +
  #geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Condition", y= 'weights of the dominant SH component') 
 


data = res[, c(47:50)]


# Create a sample dataset

# Preview the dataset
head(data)
#data = data[which(data$lmax < 5), ]

library(tidyverse)
library(htmlwidgets)
# Convert the ggplot object to a plotly object for 3D plotting
p <- plot_ly(data, x = ~r2, y = ~lmax_pct, z = ~ lmax, type = "scatter3d", mode = "markers",
             split = ~ condition,
             marker = list(size = 5)) %>%
  layout(title = "3D Scatter Plot",
         scene = list(xaxis = list(title = "R2"),
                      yaxis = list(title = "lmax_pct"),
                      zaxis = list(title = "lmax")))

# Display the 3D scatter plot
p
saveWidget(ggplotly(p), file = "myplot.html")




########################################################
########################################################
# Section II: test the comparison WT between day3, d3.5, d4, d5 and day6
# 
########################################################
########################################################
##########################################
# quickly test examples of day 3 and day6
##########################################
dataDir = "res/test_organoid_timepoints"
figureDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/",
                   "projects/RA_competence/results/figures_tables_R13547_10x_mNT_20240522")

filelist = list.files(path = dataDir, pattern = '*.csv', full.names = TRUE)

for(n in 1:length(filelist))
{
  xx = read.csv(file = filelist[n])
  cat(' nb cysts --- ', nrow(xx), '\n')
  
  if(n == 1){
    res = xx 
  }else{
    res = rbind(res, xx)
  }
  
}

aa = as.matrix(res[, grep('power_per_l_', colnames(res))])
plot(aa)

matplot(log(t(aa)), type = "l")

jj = 6

nb_l = 30
ll = c(0:nb_l)

yy = log10(aa[jj, c(1:(nb_l+1))])

plot((ll), yy,  type = 'b')

plot(2*pi/sqrt(ll*(ll+1)), yy, type = 'b')

plot((ll), log10(10^yy*(ll*(ll +1))),  type = 'b')


fit = lm(yy[15:length(yy)] ~ ll[15:length(ll)])
yy_hat = fit$coefficients[1] + ll*fit$coefficients[2]

plot((ll), yy,  type = 'b')
abline(fit$coefficients, col = 'red')

plot(ll, yy -yy_hat, type = 'b')


########################################################
########################################################
# Section II: compare all WT time points 
# 
########################################################
########################################################
dataDir = "res/test_organoid_wt_timepoints_v1"

Dirlist = list.dirs(path = dataDir, full.names = TRUE)

Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0

for(n in 1:length(Dirlist))
{
  filelist = list.files(path = Dirlist[n], pattern = '*.csv', full.names = TRUE)
  if(length(filelist) > 0){
    cat(basename(Dirlist[n]), ' :  nb cysts --- ', length(filelist), '\n')
    
    for(m in 1:length(filelist))
    {
      xx = read.csv(filelist[m])
      
      i = i + 1
      
      if(i == 1){
        res = xx 
      }else{
        res = rbind(res, xx)
      }
    }
  }
  
}

rm(xx); rm(i)

colnames(res)[1] = 'time'
res$time = sapply(res$image, function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$time[which(res$time == 'd3-5')] = 'd3.5'

rownames(res) = paste0(res$image, '_cyst_', res$cyst_index)

saveRDS(res, file = paste0(dataDir, '/saved_res.rds'))

## reload the SH outcome 
res = readRDS(file = paste0(dataDir, '/saved_res.rds'))

hist(log10(res$cyst_size), breaks = 20)

res = res[which(res$cyst_size > 10^3.75 & res$cyst_size < 10^5), ]

aa = as.matrix(res[, grep('power_per_l_', colnames(res))])
#bb = as.matrix(res[, grep('power_per_dlogl_', colnames(res))])
matplot(log(t(aa[c(1:20), c(1:30)])), type = "l")

## convert the power into contribution percentage
xx = aa
for(n in 1:nrow(aa))
{
  xx[n, ] = aa[n, ]/sum(aa[n,], na.rm = FALSE)
}


aa = xx
rm(xx)

hist(log10(aa), breaks = 100);
abline(v = c(-4, -3, -2), lwd = 2.0, col = 'red')


matplot(log10(t(aa[c(1:20), c(1:30)])), type = "l")
abline(h = c(-4, -3, -2), lwd = 2.0, col = 'red')


## example from https://stackoverflow.com/questions/61909908/
## adding-a-shaded-standard-deviation-to-line-plots-on-ggplot2-for-multiple-variabl
library(ggplot2)
library(dplyr)
library(tidyr)

aa = log10(aa)


tt = unique(res$time)

Means = matrix(NA, nrow = length(tt), ncol = 30)
colnames(Means) = colnames(aa)[1:ncol(Means)]
SD = Means

for(n in 1:length(tt))
{
  #kk = which(res$time == tt[n] & res$cyst_size > 10^4 & res$cyst_size < 10^4.2)
  kk = which(res$time == tt[n])
  Means[n, ] = apply(aa[kk, 1:ncol(Means)], 2, mean)
  SD[n, ] = apply(aa[kk, 1:ncol(SD)], 2, sd)
  
}


Means = data.frame(Means, time = tt, stringsAsFactors = FALSE)
SD = data.frame(SD, time = tt, stringsAsFactors = FALSE)

means_long <- pivot_longer(Means, -time, values_to = "mean", names_to = "variable")
sd_long <- pivot_longer(SD, -time, values_to = "sd", names_to = "variable")

df_join <- means_long %>% 
  left_join(sd_long)
#> Joining, by = c("date", "variable")

df_join$l = sapply(df_join$variable, function(x) {as.numeric(unlist(strsplit(as.character(x), '_'))[4])})

ggplot(data = df_join, aes(x = l, group = time, fill = time)) + 
  geom_line(aes(y = mean, color = time), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("harmonics degree l") + 
  ylab("% of variance ") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

ggsave(filename = paste0(figureDir, '/Mean_distribution_powerL_RA_wt.pdf'), height = 6, width = 10)

ggplot(data = df_join, aes(x = l, group = time, fill = time)) + 
  geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("harmonics degree l") + 
  ylab("% of variance ") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

ggsave(filename = paste0(figureDir, '/Mean_distribution_powerL_RA_wt_onlyMean.pdf'), height = 6, width = 10)


##########################################
# define the l_max and l_opt
##########################################
res = readRDS(file = paste0(dataDir, '/saved_res.rds'))

hist(log10(res$cyst_size))

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5), ]
res$r2 = res$cyst_r^2

aa = as.matrix(res[, grep('power_per_l_', colnames(res))])
#bb = as.matrix(res[, grep('power_per_dlogl_', colnames(res))])
matplot(log(t(aa[c(1:20), c(1:30)])), type = "l")

## convert the power into contribution percentage
xx = aa
for(n in 1:nrow(aa))
{
  xx[n, ] = aa[n, ]/sum(aa[n,], na.rm = FALSE)
}
aa = xx
rm(xx)

hist(log10(aa), breaks = 100);
abline(v = c(-4, -3, -2), lwd = 2.0, col = 'red')

aa = log10(aa)

l_max = apply(aa, 1, which.max) - 1
lmax_pct = apply(aa, 1,  function(x){x[which.max(x)]})

res$lmax = l_max
res$lmax_pct = lmax_pct

find_optimal_l = function(yy, pct_cutoff = -2)
{
  peaks = findpeaks(yy, nups = 1) 
  #yy2 = (bb[jj, c(1:(nb_l+1))])
  #yy3 = (10^yy*sqrt(ll*(ll+1)))
  index_peaks = peaks[, 2]
  index_cutoff = which(yy > pct_cutoff)
  
  index_peaks = intersect(index_peaks, index_cutoff)
  
  index_opt = index_peaks[ceiling(length(index_peaks)/2)]
  index_extrem = index_peaks[length(index_peaks)]
  
  return(c(index_opt -1, yy[index_opt], index_extrem - 1, yy[index_extrem]))
  
}

l_opt = t(apply(aa, 1, find_optimal_l))
colnames(l_opt) = c('l_optm', 'l_optm_pct', 'l_extrem', 'l_extrem_pct')

res = data.frame(res, l_opt, stringsAsFactors = FALSE)

res = res[, -grep('power_per_', colnames(res))]

#res$lmax_norm = res$lmax/res$r2
res$r2_log = log10(res$r2)
res$size_log = log10(res$cyst_size)

res$lmax_norm = log10(res$lmax/res$r2)
res$wavelength = log10(res$r2/res$lmax)
res$wavelength_v2 = 2*pi*sqrt(res$r2)/(res$lmax*(res$lmax+1))



ggplot(data = res, aes(x = lmax_norm,  y = lmax_pct, color = time)) + 
  geom_point(size = 1.6) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("normalized l (log10 lmax/surface)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(figureDir, '/normalizedLmax_vs_patternQuality_v1.pdf'), height = 6, width = 10)

## test some other normalization
res$lmax_norm2 = log10(res$lmax/res$cyst_size)
ggplot(data = res, aes(x = lmax_norm2,  y = lmax_pct, color = time)) + 
  geom_point(size = 1.6) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("normalized l (log10 lmax/surface)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


## add image name and cyst index
library(ggrepel)
res$cyst_name = rownames(res)

#xx = res[which(res$image == '211209_d3_RAd2_D2_52_01_isotropic'), ]
#xx$cyst = xx$cyst_index - 1

res$cyst = res$cyst_index

ggplot(data = res, aes(x = lmax_norm,  y = lmax_pct, color = time, label = cyst)) + 
  geom_point(size = 1.6) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("normalized l (log10 lmax/surface)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm")) +
  theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14)) +
  geom_text(aes(label=ifelse(image == '211209_d3_RAd2_D2_52_01_isotropic', as.character(cyst),'')),
            hjust=0,vjust=0)
  # geom_label_repel(aes(label = cyst),
  #                  box.padding   = 0.1, 
  #                  point.padding = 0.1,
  #                  segment.color = 'grey50')

ggsave(filename = paste0(figureDir, '/normalizedLmax_vs_patternQuality_v1.pdf'), height = 6, width = 10)


plot(res$nb_local_max, res$lmax, ylim = range(c(res$l_extrem, res$l_optm, res$lmax)), col = 'darkorange')
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$nb_local_max, res$l_optm , type = 'p', col = 'darkgreen', cex = 1.0, pch = 16)
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$nb_local_max, res$l_extrem , type = 'p', col = 'darkblue', cex = 1.0, pch = 16)
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$cyst_size, res$l_optm, type = 'p', col = 'darkblue', cex = 1.0, pch = 16)
abline(0, 1, lwd = 2.0, col = 'red')


ggplot(data = res, aes(x = size_log,  y = lmax_pct, color = time)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

ggsave(filename = paste0(figureDir, '/patternQuality_cystSize_test_v1.pdf'), height = 6, width = 10)




ggplot(data = res, aes(x = wavelength,  y = lmax_pct, color = time)) + 
  geom_point(size = 1.6) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("normalized wavelength (surface size/lmax)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
                      axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14)) 
ggsave(filename = paste0(figureDir, '/dominantWavelength_patternQuality_v1.pdf'), height = 6, width = 10)






ggplot(data = res, aes(x = size_log,  y = lmax, color = time)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2)") + 
  ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())


#xx = res[which((res$time == 'd5'|res$time == 'd6') & res$lmax <= 3), ]
xx = res

xx$r2_log2 = log2(xx$cyst_size)
ggplot(, aes(x=lmax, y=r2_log, group = lmax)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) +
  labs(x = "lmax", y= 'size in log scale') 

ggsave(filename = paste0(figureDir, '/RA_wt_d5_d6_size_vs_lmax.pdf'), height = 6, width = 10)

p2 = ggplot(data = xx, aes(x = r2_log,  y = lmax, color = time)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())






ggplot(res, aes(x=time, y=l_optm, fill=time)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) +
  labs(x = "Condition", y= 'weights of the dominant SH component') 



## example
check_individualExample = FALSE
if(check_individualExample){
  jj = 10
  head(res[jj, c(1:7)])
  
  nb_l = 30
  ll = c(0:nb_l)
  require(pracma)
 
  
  
  #plot(2*pi/sqrt(ll*(ll+1)), yy, type = 'b')
  
  yy = (aa[jj, c(1:(nb_l+1))])
  peaks = findpeaks(yy, nups = 1) 
  #yy2 = (bb[jj, c(1:(nb_l+1))])
  #yy3 = (10^yy*sqrt(ll*(ll+1)))
  index_peaks = peaks[, 2]
  
  plot(ll, yy,  type = 'b', ylim  = range(c(yy), finite = TRUE))
  points(ll[index_peaks], yy[index_peaks] , type = 'p', col = 'darkred', cex = 2.0, pch = 16)
  #points(ll, yy2, type = 'b', col = 'blue')
  #points(ll, yy3, type = 'b', col = 'red')
  abline(h = c(-3, -2), lwd = 2.0, col = 'red')
  
  
  #plot((ll), log10(10^yy*(ll +1)),  type = 'b')
  
}


##########################################
# ## prepare the data for singlecellExperiment
##########################################
Test_umap_summary = FALSE
if(Test_umap_summary){
  res = readRDS(file = paste0(dataDir, '/saved_res.rds'))
  hist(log10(res$cyst_size))
  
  res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5), ]
  res$r2 = res$cyst_r^2
  
  
  mat = res[, c(1:7)]
  mat$treatment = 'RA_wt'
  mat$condition = mat$time
  metadata = mat
  metadata$marker_class = 'type'
  
  counts = as.matrix(t(res[, grep('power_per_l_', colnames(res))]))
  
  #load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
  
  #counts = as.matrix(t(mat))
  #colnames(counts) = paste0('cell_', c(1:ncol(counts)))
  #rownames(metadata) = colnames(counts)
  metadata = data.frame(metadata, stringsAsFactors = FALSE)
  
  library(SingleCellExperiment)
  
  sce <- SingleCellExperiment(assays=list(counts=counts),
                              colData=metadata, 
                              metadata = metadata)
  
  y <- assay(sce, "counts")
  #y <- asinh(sweep(y, 1, cf, "/"))
  assay(sce, "exprs", FALSE) <- y
  
  sce$sample_id = as.character(sce$condition)
  #sce$condition = gsub('noRA_d2', "beforeRA_d2", sce$condition)
  
  cc.levels = c('d3', 'd3.5', 'd4', 'd5', 'd6')
  
  sce$condition = factor(sce$condition, levels = cc.levels)
  sce@metadata$condition = sce$condition
  
  rowData(sce)$marker_name = rownames(sce)
  rowData(sce)$channel_name = NULL
  rowData(sce)$marker_class = 'type'
  
  saveRDS(sce, file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/data/", 
                             "image_SHout_sce.rds"))
}

########################################################
########################################################
# Section III: compare the WT, KO-KO and TetOn-TetOn
# 
########################################################
########################################################
dataDir = "res/test_organoid_KOKO_TetOnTetOn_timepoints_v1"

figureDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/",
                   "projects/RA_competence/results/figures_tables_R13547_10x_mNT_20240522")

##########################################
# import the WT data 
##########################################
dataDir_WT = "res/test_organoid_wt_timepoints_v1" 

res_wt = readRDS(file = paste0(dataDir_WT, '/saved_res.rds'))
res_wt = data.frame(condition = rep('WT', nrow(res_wt)), res_wt, stringsAsFactors = FALSE)


##########################################
# process SH-analysis in the KO-KO and TetOn-TetOn conditions
##########################################
Dirlist = list.dirs(path = dataDir, full.names = TRUE)
Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0
for(n in 1:length(Dirlist))
{
  filelist = list.files(path = Dirlist[n], pattern = '*.csv', full.names = TRUE)
  if(length(filelist) > 0){
    cat(basename(Dirlist[n]), ' :  nb cysts --- ', length(filelist), '\n')
    
    for(m in 1:length(filelist))
    {
      xx = read.csv(filelist[m])
      
      i = i + 1
      
      if(i == 1){
        res = xx 
      }else{
        res = rbind(res, xx)
      }
    }
  }
  
}

rm(xx); rm(i)

colnames(res)[1] = 'time'
res$time = sapply(res$image, function(x){unlist(strsplit(as.character(x), '_'))[3]})
res$time = gsub('-2umZ', '', res$time)
res$time = gsub('-P6-FA2', '', res$time)

rownames(res) = paste0(res$image, '_cyst_', res$cyst_index)

res = data.frame(condition = rep('KO_KO', nrow(res)), res, stringsAsFactors = FALSE)
res$condition[grep('iF-iP_dox', res$image)] = 'TetOn_TetON_dox'
res$condition[grep('iF-iP_RA', res$image)] = 'TetOn_TetON_RA'
table(res$condition[grep('2umZ', res$image)]) 

saveRDS(res, file = paste0(dataDir, '/saved_res.rds'))

## reload the SH outcome 
res = readRDS(file = paste0(dataDir, '/saved_res.rds'))

## filtering the cyst by sizes
res = data.frame(rbind(res_wt, res), stringsAsFactors = FALSE)

hist(log10(res$cyst_size), breaks = 20)
res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5.5), ]

saveRDS(res, file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn.rds'))


##########################################
# add more features for the organoids 
##########################################
res = readRDS(file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn.rds'))

metaDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/',
                 'images_data/results/featureCollection_organoid_WT')

filelist = list.files(path = metaDir, pattern = '*.csv', full.names = TRUE)

metadata = c()
i = 0
if(length(filelist) > 0){
  #cat(basename(Dirlist[n]), ' :  nb cysts --- ', length(filelist), '\n')
  
  for(m in 1:length(filelist))
  {
    xx = read.csv(filelist[m])
    
    i = i + 1
    
    if(i == 1){
      metadata = xx 
    }else{
      metadata = rbind(metadata, xx)
    }
  }
}

metaDir2 = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/',
                 'images_data/results/featureCollection_organoid_timepoints_KOKO_TetOnTetOn')

filelist = list.files(path = metaDir2, pattern = '*.csv', full.names = TRUE)

metadata2 = c()
i = 0
if(length(filelist) > 0){
  #cat(basename(Dirlist[n]), ' :  nb cysts --- ', length(filelist), '\n')
  
  for(m in 1:length(filelist))
  {
    xx = read.csv(filelist[m])
    
    i = i + 1
    
    if(i == 1){
      metadata2 = xx 
    }else{
      metadata2 = rbind(metadata2, xx)
    }
  }
}

metadata$genotype_foxa2 = 1
metadata$cutoff_genotype_foxa2 = NA
metadata$cutoff_genotype_pax6 = NA
metadata$pct_foxa2 = metadata$nb_foxa2_c3/metadata$cyst_size
metadata$cutoff_foxa2 = metadata$cutoff_foxa2_c3
metadata = metadata[, -c(5:13)]


metadata2$genotype_foxa2 = metadata2$nb_genotype_fx_c2/(metadata2$nb_genotype_fx_c2 + metadata2$nb_genotype_px_c3)
metadata2$cutoff_genotype_foxa2 = metadata2$cutoff_c2
metadata2$cutoff_genotype_pax6 = metadata2$cutoff_c3
metadata2$pct_foxa2 = metadata2$nb_foxa2_c1/(metadata2$nb_foxa2_c1 + metadata2$nb_pax6_c4)
metadata2$cutoff_foxa2 = metadata2$cutoff_c1
metadata2 = metadata2[, -c(5:13)]

metadata = rbind(metadata, metadata2)
rm(metadata2)

metadata$ids = paste0(metadata$image, "_cyst_", metadata$cyst_index)

mm = match(rownames(res), metadata$ids)

xx = data.frame(res[, c(1:8)], metadata[mm, c(5:27)], res[, -c(1:8)])

res = xx; rm(xx)

saveRDS(res, file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection.rds'))

##########################################
# filter cysts with metadata  
##########################################
res = readRDS(paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection.rds'))

## define genotype_pct
res$pct_genotype = res$genotype_foxa2

jj = which(res$condition == 'WT')
res$pct_genotype[jj] = 1 - res$genotype_foxa2[jj]

res$r2 = res$cyst_r^2

res = data.frame(res)

res$outliers = FALSE

## filtered the segmentation with sphericity
ggplot(data = res, aes(x = cyst_size,  y = sphericity_wadell, color = time, group = condition)) + 
  geom_point(size = 1.6, aes(shape = condition)) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  #xlab("normalized l (log10 lmax/surface)") + 
  #ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  scale_x_continuous(trans='log10') +
  #scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

jj = which(res$sphericity_wadell >=0.75)
res = res[jj, ]


## filter cyst with no foxa2 
ggplot(data = res, 
       aes(x = pct_foxa2,  y = intensity_max_foxa2, color = time)) + 
  geom_point(size = 1.6, aes(shape = condition)) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  #xlab("normalized l (log10 lmax/surface)") + 
  #ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  geom_vline(xintercept = 0.01, linetype="solid", color = "blue") +
  geom_hline(yintercept = 100, linetype="solid", color = "blue") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.2,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


jj = which(res$pct_foxa2 > 0.01 & res$intensity_mean_foxa2 > 100)
res = res[jj, ]

## check the genotype pct quantification
ggplot(data = res, 
       aes(x = cutoff_genotype_foxa2,  y = cutoff_genotype_pax6, color = condition)) + 
  geom_point(size = 1.6, aes(shape = condition)) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  #xlab("normalized l (log10 lmax/surface)") + 
  #ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
  theme_bw() +  
  geom_vline(xintercept = 0.01, linetype="solid", color = "blue") +
  geom_hline(yintercept = 100, linetype="solid", color = "blue") +
  #scale_x_continuous(trans='log10') +
  #scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.9,.8), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

jj = which(res$condition == "WT" | (res$condition != "WT"  & res$cutoff_genotype_foxa2 < 2000))
res = res[jj, ]

saveRDS(res, file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered.rds'))


##########################################
# ## prepare the data for singlecellExperiment
##########################################
Test_umap_summary = FALSE
if(Test_umap_summary){
  library(SingleCellExperiment)
  #library(CATALYST)
  #res = readRDS(file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn.rds'))
  
  res = readRDS(file = paste0(dataDir, '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered.rds'))
  hist(log10(res$cyst_size))
  
  res = res[which(res$sphericity_wadell > 0.8), ]
  
  table(res$condition, res$time)
  
  ## further filter cysts with sphericity and FoxA2 sd 
  res$intensity_std_foxa2_log = log10(res$intensity_std_foxa2)
  
  ggplot(res, aes(x=intensity_std_foxa2_log, y=sphericity_wadell, color=time, shape = condition)) + 
    geom_point(size = 1.5) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())
  
  res$cyst_size_log = log10(res$cyst_size)
  
  jj1 = which(res$condition == 'WT')
  jj2 = which(res$condition == 'KO_KO')
  jj3 = which(res$condition == 'TetOn_TetON_RA')
  jj4 = which(res$condition == 'TetOn_TetON_dox')
  
  ggplot(res[c(jj1), ], aes(x=cyst_size_log, y=intensity_std_foxa2_log)) + 
    geom_point(size = 1.5) +
    geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 

  keep = c()
  ## WT
  xx = res[which(res$condition == 'WT' & res$cyst_size_log > 3.75 & res$cyst_size_log <5), ]
  ggplot(xx, aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time)) + 
    geom_point(size = 1.5) +
    #geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 
  
  #fit = loess(xx$intensity_std_foxa2_log ~ xx$cyst_size_log)
  #hist(fit$residuals, breaks = 100)
  #xx = xx[which(abs(fit$residuals) < 0.5), ]
  xx = xx[which((xx$intensity_std_foxa2_log > 2.0 & (xx$time == 'd3'|xx$time == 'd3.5'|xx$time == 'd4')) |
                  ((xx$time == 'd5'|xx$time == 'd6') & xx$intensity_std_foxa2_log > 2.25)), ]
  
  
  keep = c(keep, rownames(xx))
  
  
  ## KO-KO
  xx = res[which(res$condition == 'KO_KO' & res$cyst_size_log > 4.0), ]
  ggplot(xx, aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time)) + 
    geom_point(size = 1.5) +
    #geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 
  
  test = xx[which(xx$intensity_std_foxa2_log < 2.5), ]
  xx = xx[which(xx$intensity_std_foxa2_log > 2.5), ]
  keep = c(keep, rownames(xx))
  
  ## TetOn_TetON_RA
  xx = res[which(res$condition == 'TetOn_TetON_RA' & res$cyst_size_log > 3.75), ]
  ggplot(xx, aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time)) + 
    geom_point(size = 1.5) +
    #geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 
  
  test = xx[which(xx$intensity_std_foxa2_log < 2.5), ]
  xx = xx[which(xx$intensity_std_foxa2_log > 2.5), ]
  keep = c(keep, rownames(xx))
  
  
  ## TetOn_TetON_dox
  xx = res[which(res$condition == 'TetOn_TetON_dox' & res$cyst_size_log > 3.75), ]
  ggplot(xx, aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time)) + 
    geom_point(size = 1.5) +
    #geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 
  
  test = xx[which(xx$intensity_std_foxa2_log < 2.5), ]
  xx = xx[which(xx$intensity_std_foxa2_log > 2.5), ]
  keep = c(keep, rownames(xx))
  
  
  res = res[which(!is.na(match(rownames(res), keep))), ]
  
  saveRDS(res, file = paste0(dataDir, 
                             '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered_FoxA2SDfiltered.rds'))
  
  #res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5), ]
  res$r2 = res$cyst_r^2
  
  mat = res[, c(1:31)]
  #mat$treatment = 'RA_wt'
  #mat$condition = mat$time
  metadata = mat
  metadata$marker_class = 'type'
  metadata = data.frame(metadata, stringsAsFactors = FALSE)
  
  counts = as.matrix(t(res[, grep('power_per_l_', colnames(res))]))
  #counts = as.matrix(t(res[, c(grep('power_per_dlogl_', colnames(res))[1:20], 9:15, 17:18, 22:24)]))
  
  sce <- SingleCellExperiment(assays=list(counts=counts),
                              colData=metadata, 
                              metadata = metadata)
  
  y <- assay(sce, "counts")
  #y <- asinh(sweep(y, 1, cf, "/"))
  assay(sce, "exprs", FALSE) <- y
  
  sce$sample_id = as.character(sce$condition)
  #sce$condition = gsub('noRA_d2', "beforeRA_d2", sce$condition)
  
  cc.levels = levels = c("WT", "KO_KO", "TetOn_TetON_RA", 'TetOn_TetON_dox')
  
  sce$condition = factor(sce$condition, levels = cc.levels)
  sce@metadata$condition = sce$condition
  
  rowData(sce)$marker_name = rownames(sce)
  rowData(sce)$channel_name = NULL
  rowData(sce)$marker_class = 'type'
  
  saveRDS(sce, file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/data/", 
                             "image_SHout_dlogl_sce_wt_conditions_metadata_filteringFoxA2SD_v3.rds"))
  
  table(sce$condition, sce$time)
  
}

##########################################
# define the l_max and l_opt
##########################################
ExploraryAnalaysis_by_definingL_opt = FALSE
if(ExploraryAnalaysis_by_definingL_opt)
{
  res = readRDS(file = paste0(dataDir, 
                              '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered_FoxA2SDfiltered.rds'))
  hist(log10(res$cyst_size))
  
  res = res[which(res$time != "d3.5" & res$time != "d5"), ]
  res$r2 = res$cyst_r^2
  
  cc.levels = levels = c("WT", "KO_KO", "TetOn_TetON_RA", 'TetOn_TetON_dox')
  res$condition = factor(res$condition, levels = cc.levels)
  
  jj1 = which(res$condition == 'WT')
  jj2 = which(res$condition == 'KO_KO')
  jj3 = which(res$condition == 'TetOn_TetON_RA')
  jj4 = which(res$condition == 'TetOn_TetON_dox')
  
  jj = which(res$time == 'd3' & res$condition != 'TetOn_TetON_dox')
  ggplot(res[c(jj1, jj2), ],  aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time, shape = condition)) + 
    geom_point(size = 2.0) +
    #geom_smooth(method=loess, aes(x=cyst_size_log, y=intensity_std_foxa2_log))+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    #geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    #geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank()) 
  
  table(res$condition, res$time)
  
  
  aa = as.matrix(res[, grep('power_per_l_', colnames(res))])
  #bb = as.matrix(res[, grep('power_per_dlogl_', colnames(res))])
  matplot(log(t(aa[c(1:100), ])), type = "l")
  
  ## convert the power into contribution percentage
  xx = aa
  for(n in 1:nrow(aa))
  {
    xx[n, ] = aa[n, ]/sum(aa[n,], na.rm = FALSE)
  }
  aa = xx
  rm(xx)
  
  aa = aa[, -1]
  
  xx = aa
  for(n in 1:nrow(xx))
  {
    xx[n, ] = 2*pi*sqrt(res$r2[n])/sqrt(c(1:100) * (c(1:100) +1))
  }
  
  plot(c(1,1), type = 'n', xlim = c(1, 150), ylim = range(aa), log = 'xy')
  
  nb_sample = 50
  res$cc = paste0(res$condition, '_', res$time)
  index = c(which(res$cc == 'WT_d3')[1:nb_sample], 
            which(res$cc == 'KO_KO_d6')[1:nb_sample],
            which(res$cc == 'WT_d6')[1:nb_sample])
  cols = rep(c('gray', 'orange', 'darkblue'), each = nb_sample)
  for(n in 1:length(index))
  {
    points(xx[index[n], ], aa[index[n],], col = cols[n], type = 'l')
  }
  
  #matplot(log(t(xx[c(1:100), ])), type = "l")
  
  table(res$condition, res$time)
  
  ll = c(1:100)
  
  plot(1, 1, type = 'n', col = 1, log = 'y', xlim = c(0, 1.5), ylim = range(aa))
  
  jj = which(res$condition == "WT" & res$time == "d3")
  
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  
  points(xx[-1], yy, type = 'l', col = 'dodgerblue', lwd = 3.0)
  
  
  jj = which(res$condition == "WT" & res$time == "d4")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'blue', lwd = 3.0)
  
  
  jj = which(res$condition == "WT" & res$time == "d6")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'darkblue', lwd = 3.0)
  
  
  jj = which(res$condition == "KO_KO" & res$time == "d3")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'lightpink', lwd = 3.0)
  
  jj = which(res$condition == "KO_KO" & res$time == "d4")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'deeppink', lwd = 3.0)
  
  
  jj = which(res$condition == "KO_KO" & res$time == "d6")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'darkred', lwd = 3.0)
  
  
  ## TetOn_TetOn RA
  cc = "TetOn_TetON_RA"
  jj = which(res$condition == "TetOn_TetON_RA" & res$time == "d3")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'chartreuse', lwd = 3.0)
  
  jj = which(res$condition == "TetOn_TetON_RA" & res$time == "d4")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'forestgreen', lwd = 3.0)
  
  
  jj = which(res$condition == "TetOn_TetON_RA" & res$time == "d6")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'darkgreen', lwd = 3.0)
  
  ## TetOn-TetOn dox
  cc = "TetOn_TetON_RA"
  jj = which(res$condition == "TetOn_TetON_dox" & res$time == "d3")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'blueviolet', lwd = 3.0)
  
  jj = which(res$condition == "TetOn_TetON_dox" & res$time == "d4")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'magenta', lwd = 3.0)
  
  
  jj = which(res$condition == "TetOn_TetON_dox" & res$time == "d6")
  ctb = aa[jj, ]
  lambda = matrix(0, nrow = nrow(ctb), ncol = ncol(ctb))
  for(n in 1:nrow(lambda))
  {
    lambda[n, ] =  1/sqrt(ll*(ll+1)) * res$cyst_r[jj[n]]
  }
  xx = seq(log10(min(lambda, na.rm = TRUE)), log10(max(lambda, na.rm = TRUE)), length.out = 101)
  yy = rep(NA, length(xx)-1)
  for(n in 1:length(yy))
  {
    # n = 1
    kk = which(lambda >= 10^xx[n] & lambda < 10^xx[n+1])
    yy[n] = mean(ctb[kk])
  }
  points(xx[-1], yy, type = 'l', col = 'maroon4', lwd = 3.0)
  
  
  hist(log10(aa), breaks = 100);
  abline(v = c(-4, -3, -2), lwd = 2.0, col = 'red')
  
  aa = log10(aa)
  
  l_max = apply(aa, 1, which.max) - 1
  lmax_pct = apply(aa, 1,  function(x){x[which.max(x)]})
  
  res$lmax = l_max
  res$lmax_pct = lmax_pct
  
  find_optimal_l = function(yy, pct_cutoff = -2)
  {
    peaks = pracma::findpeaks(yy, nups = 1)
    #yy2 = (bb[jj, c(1:(nb_l+1))])
    #yy3 = (10^yy*sqrt(ll*(ll+1)))
    index_peaks = peaks[, 2]
    index_cutoff = which(yy > pct_cutoff)
    
    index_peaks = intersect(index_peaks, index_cutoff)
    
    index_opt = index_peaks[ceiling(length(index_peaks)/2)]
    index_extrem = index_peaks[length(index_peaks)]
    
    return(c(index_opt -1, yy[index_opt], index_extrem - 1, yy[index_extrem]))
    
  }
  
  l_opt = t(apply(aa, 1, find_optimal_l))
  colnames(l_opt) = c('l_optm', 'l_optm_pct', 'l_extrem', 'l_extrem_pct')
  
  res = data.frame(res, l_opt, stringsAsFactors = FALSE)
  
  res = res[, -grep('power_per_', colnames(res))]
  
  #res$lmax_norm = res$lmax/res$r2
  res$r2_log = log10(res$r2)
  res$size_log = log10(res$cyst_size)
  
  #res$lmax_norm = log10(res$lmax/res$r2)
  #res$wavelength = log10(res$r2/res$lmax)
  res$wavelength_v2 = 2*pi*sqrt(res$r2)/(res$lmax*(res$lmax+1))
  
  res$condition = factor(res$condition, levels = c("WT", "KO_KO", "TetOn_TetON_RA",  'TetOn_TetON_dox'))
  
  res = res[grep('TetOn_TetON', res$condition), ]
  
  ggplot(data = res[c(jj3, jj4), ], aes(x = lmax_norm,  y = lmax_pct, color = time, group = condition)) + 
    geom_point(size = 1.6, aes(shape = condition)) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("normalized l (log10 lmax/surface)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  ggsave(filename = paste0(figureDir, '/normalizedLmax_vs_patternQuality_WTvsConditions_v1.pdf'), height = 6, width = 10)
  
  
  ## test some other normalization
  res$lmax_norm2 = log10(res$lmax/res$cyst_size)
  ggplot(data = res, aes(x = lmax_norm2,  y = lmax_pct, color = time, group = condition)) + 
    geom_point(size = 1.6, aes(shape = condition)) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("normalized l (log10 lmax/surface)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  
  ## add image name and cyst index
  library(ggrepel)
  res$cyst_name = rownames(res)
  
  #xx = res[which(res$image == '211209_d3_RAd2_D2_52_01_isotropic'), ]
  #xx$cyst = xx$cyst_index - 1
  
  res$cyst = res$cyst_index
  
  ggplot(data = res, aes(x = lmax_norm,  y = lmax_pct, color = time, label = cyst)) + 
    geom_point(size = 1.6) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("normalized l (log10 lmax/surface)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm")) +
    theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14)) +
    geom_text(aes(label=ifelse(image == '211209_d3_RAd2_D2_52_01_isotropic', as.character(cyst),'')),
              hjust=0,vjust=0)
  # geom_label_repel(aes(label = cyst),
  #                  box.padding   = 0.1, 
  #                  point.padding = 0.1,
  #                  segment.color = 'grey50')
  
  ggsave(filename = paste0(figureDir, '/normalizedLmax_vs_patternQuality_v1.pdf'), height = 6, width = 10)
  
  
  plot(res$nb_local_max, res$lmax, ylim = range(c(res$l_extrem, res$l_optm, res$lmax)), col = 'darkorange')
  abline(0, 1, lwd = 2.0, col = 'red')
  
  plot(res$nb_local_max, res$l_optm , type = 'p', col = 'darkgreen', cex = 1.0, pch = 16)
  abline(0, 1, lwd = 2.0, col = 'red')
  
  plot(res$nb_local_max, res$l_extrem , type = 'p', col = 'darkblue', cex = 1.0, pch = 16)
  abline(0, 1, lwd = 2.0, col = 'red')
  
  plot(res$cyst_size, res$l_optm, type = 'p', col = 'darkblue', cex = 1.0, pch = 16)
  abline(0, 1, lwd = 2.0, col = 'red')
  
  
  ggplot(data = res, aes(x = size_log,  y = lmax_pct, color = time)) + 
    geom_point(size = 1.5) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("cyst surface (r2)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())
  
  ggsave(filename = paste0(figureDir, '/patternQuality_cystSize_test_v1.pdf'), height = 6, width = 10)
  
  
  
  
  ggplot(data = res, aes(x = wavelength,  y = lmax_pct, color = time)) + 
    geom_point(size = 1.6) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("normalized wavelength (surface size/lmax)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.8), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14)) 
  
  ggsave(filename = paste0(figureDir, '/dominantWavelength_patternQuality_v1.pdf'), height = 6, width = 10)
  
  
  
  
  
  
  ggplot(data = res, aes(x = size_log,  y = lmax, color = time)) + 
    geom_point(size = 1.5) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("cyst surface (r2)") + 
    ylab("pattern quality (log10 % variance by dominant SH degree) ") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())
  
  
  #xx = res[which((res$time == 'd5'|res$time == 'd6') & res$lmax <= 3), ]
  xx = res
  
  xx$r2_log2 = log2(xx$cyst_size)
  ggplot(, aes(x=lmax, y=r2_log, group = lmax)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=10), 
          legend.title = element_text(size = 10)) +
    labs(x = "lmax", y= 'size in log scale') 
  
  ggsave(filename = paste0(figureDir, '/RA_wt_d5_d6_size_vs_lmax.pdf'), height = 6, width = 10)
  
  p2 = ggplot(data = xx, aes(x = r2_log,  y = lmax, color = time)) + 
    geom_point(size = 1.5) +
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    xlab("cyst surface (r2 log2)") + 
    ylab("lmax") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())
  
  
  
  
  
  
  ggplot(res, aes(x=time, y=l_optm, fill=time)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=10), 
          legend.title = element_text(size = 10)) +
    labs(x = "Condition", y= 'weights of the dominant SH component') 
  
  
  
  ## example
  check_individualExample = FALSE
  if(check_individualExample){
    jj = 10
    head(res[jj, c(1:7)])
    
    nb_l = 30
    ll = c(0:nb_l)
    require(pracma)
    
    
    
    #plot(2*pi/sqrt(ll*(ll+1)), yy, type = 'b')
    
    yy = (aa[jj, c(1:(nb_l+1))])
    peaks = findpeaks(yy, nups = 1) 
    #yy2 = (bb[jj, c(1:(nb_l+1))])
    #yy3 = (10^yy*sqrt(ll*(ll+1)))
    index_peaks = peaks[, 2]
    
    plot(ll, yy,  type = 'b', ylim  = range(c(yy), finite = TRUE))
    points(ll[index_peaks], yy[index_peaks] , type = 'p', col = 'darkred', cex = 2.0, pch = 16)
    #points(ll, yy2, type = 'b', col = 'blue')
    #points(ll, yy3, type = 'b', col = 'red')
    abline(h = c(-3, -2), lwd = 2.0, col = 'red')
    
    
    #plot((ll), log10(10^yy*(ll +1)),  type = 'b')
    
  }
  
}



##########################################
# analyze the UMAP structure of cyst across conditions
##########################################
sce = readRDS(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/data/", 
                            "image_SHout_dlogl_sce_wt_conditions_metadata_pca_umap_v6.rds"))


rds = reducedDim(sce, 'UMAP')
colnames(rds) = c('umap_1', 'umap_2')
pcs = reducedDim(sce, 'PCA')[, c(1:2)]
rds = data.frame(rds, pcs, colData(sce))


ggplot(data = rds, aes(x = PC1,  y = PC2, color = time, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("PC1") + 
  ylab("PC2") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())



jj = which(rds$condition == 'WT')
p1 = ggplot(data = rds[jj, ], aes(x = PC1,  y = PC2, color = time, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

jj = which(rds$condition == 'KO_KO')
p2 = ggplot(data = rds[jj, ], aes(x = umap_1,  y = umap_2, color = time, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

jj = which(rds$condition == 'TetOn_TetON_RA')
p3 = ggplot(data = rds[jj, ], aes(x = umap_1,  y = umap_2, color = time, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())


jj = which(rds$condition == 'TetOn_TetON_dox')
p4 = ggplot(data = rds[jj, ], aes(x = umap_1,  y = umap_2, color = time, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

(p1 + p2)/(p3 + p4)

jj1 = which(rds$condition == 'WT')
jj2 = which(rds$condition == 'KO_KO')
jj3 = which(rds$condition == 'TetOn_TetON_RA')
jj4 = which(rds$condition == 'TetOn_TetON_dox')
rds$intensity_std_foxa2 = log10(rds$intensity_std_foxa2)

library(RColorBrewer)
#colfunc<-colorRampPalette(c("#4CC9F0", "yellow","#C61010"))
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "RdYlBu")))


ggplot(data = rds[which(rds$intensity_std_foxa2 >2.5), ], 
       aes(x = umap_1,  y = umap_2, color = intensity_std_foxa2, shape = condition)) + 
  geom_point(size = 1.5) +
  #geom_raster()  + 
  #geom_line(aes(y = mean, color = time), size = 1) + 
  #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
  xlab("cyst surface (r2 log2)") + 
  ylab("lmax") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank()) +
  scale_colour_gradientn(colours = SpatialColors(n = 21))

  #scale_colour_gradientn(colours = c("purple", "orange")) 

xx = (rds[which(rds$umap_1 > 7.5 & rds$time == "d6"), ])


