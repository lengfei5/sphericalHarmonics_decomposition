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



