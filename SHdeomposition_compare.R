rm(list =ls())

require(ggplot2)
library(plotly)

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




