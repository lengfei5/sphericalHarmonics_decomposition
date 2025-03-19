rm(list =ls())

dataDir = "res/test_organoid_WTd6"
Dirlist = list.dirs(path = dataDir, full.names = TRUE)

Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0
for(n in 1:length(Dirlist))
{
  if(file.exists(paste0(Dirlist[n], '/Condion_echOrgnoid_SHspectrum_power.csv'))){
    i = i + 1
    xx = read.csv(paste0(Dirlist[n], '/Condion_echOrgnoid_SHspectrum_power.csv'))
    if(i == 1){
      res6 = xx 
    }else{
      res6 = rbind(res6, xx)
    }
  }
  
}

dataDir = "res/test_organoid_WTd4"
Dirlist = list.dirs(path = dataDir, full.names = TRUE)

Dirlist = Dirlist[grep('isotropic', Dirlist)]

i = 0
for(n in 1:length(Dirlist))
{
  if(file.exists(paste0(Dirlist[n], '/Condion_echOrgnoid_SHspectrum_power.csv'))){
    i = i + 1
    xx = read.csv(paste0(Dirlist[n], '/Condion_echOrgnoid_SHspectrum_power.csv'))
    if(i == 1){
      res4 = xx 
    }else{
      res4 = rbind(res4, xx)
    }
  }
}

res6$condition = "wt_d6" 
res4$condition = "wt_d4" 
res = rbind(res6, res4)

res = res[, -1]

spect = as.matrix(res[, grep('power_per_lm', colnames(res))[1:21]])
l_max = apply(spect, 1, which.max) - 1
lmax_pct = apply(spect, 1,  function(x){x[which.max(x)]/sum(x)})

res$lmax = l_max
res$lmax_pct = lmax_pct
res$r2 = res$cyst_r^2
res$lmax_norm = res$lmax/res$r2

require(ggplot2)

res = data.frame(res, stringsAsFactors = TRUE)
res = res[which(res$cyst_size > 10^3 & res$lmax <=5), ]

ggplot(res, aes(x=r2, y=lmax, color=condition)) +
  geom_point()


ggplot(res, aes(x=lmax, y=lmax_pct, color=condition)) +
  geom_point()

data = res[, c(47:50)]

library(plotly)
# Create a sample dataset

# Preview the dataset
head(data)
data = data[which(data$lmax < 5), ]


# Convert the ggplot object to a plotly object for 3D plotting
p <- plot_ly(data, x = ~lmax, y = ~lmax_pct, z = ~ r2, type = "scatter3d", mode = "markers",
             split = ~ condition,
             marker = list(size = 7)) %>%
  layout(title = "3D Scatter Plot",
         scene = list(xaxis = list(title = "lmax"),
                      yaxis = list(title = "lmax_pct"),
                      zaxis = list(title = "R2")))

# Display the 3D scatter plot
p



