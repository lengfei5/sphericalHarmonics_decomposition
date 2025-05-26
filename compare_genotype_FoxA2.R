##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: compare the genotype and FoxA2+ cells after image analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr 22 14:02:11 2025
##########################################################################
##########################################################################

rm(list = ls())

outDir = "./results/"
if(!dir.exists(outDir)) dir.create(outDir)

library(ggplot2)
library(plotly)
library(patchwork)

########################################################
########################################################
# Section I : distribution of FoxA2+ in WT (figure 1K and 1L) 
# 
########################################################
########################################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                             "images_data/results/test_WT_FoxA2_pct/",
                             "wt_FoxA2_Pax6_voxel_counts_global_cyst_otsu_li_mean_localThreshold.csv"), 
               header = TRUE, row.names = c(1))

res$pct_foxa2_otsuCyst = res$nb_foxa2_otsu_cyst/(res$nb_foxa2_otsu_cyst + res$nb_pax6_otsu_cyst + res$nb_double_otsu_cyst)
res$pct_foxa2_otsuGlobal = res$nb_foxa2_otsu_global/(res$nb_foxa2_otsu_global + 
                                                       res$nb_pax6_otsu_global + res$nb_double_otsu_global)

res$pct_foxa2_liCyst = res$nb_foxa2_li_cyst/(res$nb_foxa2_li_cyst + res$nb_pax6_li_cyst + res$nb_double_li_cyst)
res$pct_foxa2_liGlobal = res$nb_foxa2_li_global/(res$nb_foxa2_li_global + res$nb_pax6_li_global +
                                                   res$nb_double_li_global)

res$pct_foxa2_local = res$nb_foxa2_localThreshold/(res$nb_foxa2_localThreshold + res$nb_pax6_localThreshold +
                                                     res$nb_double_localThreshold)

res$pct_foxa2_mean = res$nb_foxa2_mean_cyst/(res$nb_foxa2_mean_cyst + res$nb_pax6_mean_cyst + 
                                               res$nb_double_mean_cyst)


USE_OTSU_cyst = FALSE
if(USE_OTSU_cyst){
  res$pct_foxa2 = res$pct_foxa2_otsuCyst
  res$pct_double = res$nb_double_otsu_cyst/(res$nb_foxa2_otsu_cyst + res$nb_pax6_otsu_cyst + res$nb_double_otsu_cyst)
  res$cutoff_foxa2 = res$cutoff_otsu_cyst_foxa2
  res$cutoff_pax6 = res$cutoff_otsu_cyst_pax6
  
  plot(res$cutoff_foxa2, res$cutoff_pax6)
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  
  hist(res$cutoff_otsu_cyst_foxa2, breaks = 40)
  abline(v = c(600, 1200), lwd = 2.0, col = 'red')
  
  res = res[which(res$cutoff_foxa2 > 600 & res$cutoff_foxa2 < 1200), ]
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  hist(res$cutoff_pax6, breaks = 40)
  abline(v = c(1500, 5500), lwd = 2.0, col = 'red')
  
  res = res[which(res$cutoff_pax6 > 1500 & res$cutoff_foxa2 < 5500), ]
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  plot(res$cyst_size, res$pct_foxa2)
  
  plot(res$cyst_size, res$pct_double) ## here total is the number of FoxA2+, Pax6+ and FoxA2+&Pax6+ 
  abline(a = 0, b = 0.5, lwd = 2.0, col = 'red')
  
  res = res[which(res$cyst_size > 10000), ]
  
  plot(density(res$pct_foxa2, adjust = 1.2), col = 'darkgreen', lwd = 3.0,  xlim = c(0, 1), ylim=c(0, 3.5))
  abline(v = 0.3, lwd = 2.0)
  
  
}else{
  
  res$pct_foxa2 = res$pct_foxa2_mean
  res$pct_double = res$nb_double_mean_cyst/(res$nb_foxa2_mean_cyst + res$nb_pax6_mean_cyst + 
                                              res$nb_double_mean_cyst)
  
  res$cutoff_foxa2 = res$cutoff_mean_cyst_foxa2
  res$cutoff_pax6 = res$cutoff_mean_cyst_pax6
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  res = res[which(res$cyst_size > 10000), ]
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  kk = which(res$pct_foxa2 < 0.2)
  plot(res$cutoff_foxa2, res$cutoff_pax6)
  abline(v = c(450, 1000), lwd = 2.0, col = 'red')
  abline(h = c(2500, 5000), lwd = 2.0, col = 'red')
  points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'orange', pch = 16)
  
  hist(res$cutoff_otsu_cyst_foxa2, breaks = 40)
  abline(v = c(500, 1200), lwd = 2.0, col = 'red')
  
  plot(res$cyst_size, res$pct_foxa2)
  kk = which(res$cutoff_foxa2 < 450 | res$cutoff_foxa2 > 1000 | res$cutoff_pax6 < 2500 | res$cutoff_pax6 > 5000)
  points(res$cyst_size[kk], res$pct_foxa2[kk], pch = 16, col = 'red')
  
  res = res[which(res$cutoff_foxa2 > 450 & res$cutoff_foxa2 < 1000 & res$cutoff_pax6 > 2500 & res$cutoff_pax6 < 5000), ]
  
  # plot(res$cyst_size, res$pct_foxa2)
  # abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  # 
  # hist(res$cutoff_pax6, breaks = 40)
  # abline(v = c(1500, 5000), lwd = 2.0, col = 'red')
  # 
  # res = res[which(res$cutoff_pax6 > 1500 & res$cutoff_foxa2 < 5000), ]
  # 
  # plot(res$cyst_size, res$pct_foxa2)
  # abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  # 
  plot(res$cyst_size, res$pct_foxa2)
  
  plot(res$cyst_size, res$pct_double) ## here total is the number of FoxA2+, Pax6+ and FoxA2+&Pax6+ 
  abline(a = 0, b = 0.5, lwd = 2.0, col = 'red')
  
  plot(density(res$pct_foxa2, adjust = 1.5), col = 'darkgreen', lwd = 3.0,  xlim = c(0, 1), ylim=c(0, 6))
  abline(v = 0.28, lwd = 2.0)
  
  df = data.frame(pct = res$pct_foxa2, condition = rep('wt', nrow(res)))
  
  ggplot(df, aes(x=pct, fill = condition)) +
    geom_density(alpha=0.7, adjust = 1.5) + 
    xlim(0, 1) + 
    scale_fill_manual(values=c("darkgreen")) + 
    xlab("% FoxA2+ ") + 
    ylab("Density") + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  ggsave(filename = paste0(outDir, 'RA_WT_day4_FoxA2pct_thresholdMean.pdf'), height = 6, width = 6)
  
  
}


hist(res$pct_double, breaks = 40)
res = res[which(res$pct_double <0.3), ]

hist(res$pct_foxa2, breaks = 20)
plot(density(res$pct_foxa2, adjust = 1.5))


xx = log10(res$cyst_size[which(res$cyst_size>10^4)])
hist(xx,  probability = TRUE)
lines(density(xx), lwd = 2, col = "chocolate3")

rd = rnorm(n = 500, mean = mean(xx), sd = sd(xx))
lines(density(rd), lwd = 2, col = "red")

nb_cyst = 5000;

set.seed(2025)
a = mean(xx)/log10(200)
rd2 = rnorm(n = nb_cyst, mean = mean(xx)/a, sd = sd(xx)/a)
rd2 = 10^rd2
rd_size = floor(rd2)
hist(rd_size, breaks = 20)

set.seed(2025)
rd_cells = rep(0, sum(rd_size))
rd_cells[sample(c(1:length(rd_cells)), size = floor(length(rd_cells)*0.3), replace = FALSE)] = 1

pct = c()

for(size in rd_size)
{
  #prob = runif(n = 1, min = 0, max = 1);
  #pct = c(pct, rbinom(n = 1, size = size, prob = prob)/size)
  #pct2 = c(pct2, rbinom(n = 1, size = size, prob = 0.3)/size)
  cyst = rd_cells[1:size]
  rd_cells = rd_cells[-c(1:size)]
  pct = c(pct, sum(cyst)/length(cyst))
  
}


plot(density(res$pct_foxa2, adjust = 1.5), col = 'darkgreen', lwd = 3.0,  xlim = c(0, 1), ylim=c(0, 12))
lines(density(pct, adjust = 1.2), col = 'black', lwd = 2.0)
#lines(density(pct2), col = 'darkorange', lwd = 2.0)

df = rbind(cbind(rep('wt', nrow(res)), res$pct_foxa2), 
           cbind(rep('random', length(pct)), pct))
df = data.frame(df, stringsAsFactors = FALSE)
colnames(df) = c('condition', 'pct')
df$pct = as.numeric(df$pct)
df$condition = factor(df$condition, levels = c('wt', 'random'))

ggplot(df, aes(x=pct, fill=condition)) +
  geom_density(alpha=0.7, adjust = 1.5) +
  scale_fill_manual(values=c("darkgreen", "#999999")) + 
  xlab("% FoxA2+ ") + 
  ylab("Density") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_WT_day4_FoxA2pct.pdf'), height = 6, width = 6)


df = rbind(cbind(rep('wt', nrow(res)), res$pct_foxa2), 
           cbind(rep('random', length(pct)), pct),
           cbind(rep('perfectReg', length(pct2)), pct2)
           )

df = data.frame(df, stringsAsFactors = FALSE)
colnames(df) = c('condition', 'pct')
df$pct = as.numeric(df$pct)
df$condition = factor(df$condition, levels = c('wt', 'random', 'perfectReg'))

ggplot(df, aes(x=pct, fill=condition)) +
  geom_density(alpha=0.5, adjust = 1.5) +
  scale_fill_manual(values=c("darkgreen", "#999999", "darkorange")) + 
  xlab("% FoxA2+ ") + 
  ylab("Density") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_WT_day4_FoxA2pct_random_perfectRegulation.pdf'), height = 6, width = 6)


plot(res$cyst_size, res$total_postives)
abline(0, 0.5, lwd = 2.0, col = 'red')

plot(res$cyst_size * (res$nb_foxa2_postive + res$nb_double_positive)/res$total_postives, 
     (res$nb_foxa2_postive + res$nb_double_positive))
abline(0, 0.5, lwd = 2.0, col = 'red')

plot(res$cyst_size * (res$nb_foxa2_postive)/res$total_postives, 
     res$nb_foxa2_postive)
abline(0, 0.5, lwd = 2.0, col = 'red')


########################################################
########################################################
# Section II: chimeric FoxA2-KO & Pax6-KO
# 
########################################################
########################################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "PKO_FKO_d3_4_5_genotype_FoxA2_Pax6/cyst_size_genotype_FoxA2_Pax6_cystThresholds.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})
res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-2umZ','', res$time)

#res = res[which(res$time == 'd4'), ]
res = res[which(res$time == 'd4' & res$treatment == 'RA'), ]
#res = res[which(res$treatment == 'RA'), ]

res$genotype_pct_pxko = res$nb_pxko/res$genotype_total

res$pct_foxa2 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/res$genotype_total*2

res$pct_foxa22 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/(res$fxko_nb_foxa2+res$pxko_nb_foxa2 + 
                                                          res$fxko_nb_pax6+res$pxko_nb_pax6)

res$pct_foxa2_fxko = res$fxko_nb_foxa2/res$nb_fxko
res$pct_pax6_pxko = res$pxko_nb_pax6/res$nb_pxko


## size filtering 
hist(log10(res$cyst_size), breaks = 50)

res = res[which(res$cyst_size > 10^4), ]

plot(res$cutoff_foxa2, res$cutoff_pax6, xlim = c(0, 2), ylim = c(0, 2))
#jj = which(res$treatment == 'noRA')
#points(res$cutoff_foxa2[jj], res$cutoff_pax6[jj], col = 'red')
abline(v = 0.25, col = 'red')
abline(h = 0.25, col = 'red')

res = res[which(res$cutoff_foxa2 < 2), ]

hist(res$cutoff_foxa2, breaks = 40)
abline(v = c(0.25, 1.7), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_foxa2 > 0.25 & res$cutoff_foxa2 < 1.7), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$genotype_pct_pxko, res$pct_foxa22)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_pax6, breaks = 40)
abline(v = c(0.25, 1.2), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_pax6 > 0.25 & res$cutoff_foxa2 < 1.0), ]

#plot(res$genotype_pct_pxko, res$pct_foxa2)
#abline(0, 1, lwd = 2.0, col = 'red')

plot(res$genotype_pct_pxko, res$pct_foxa22)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$pct_foxa2_fxko, breaks = 40)

res = res[which(res$pct_foxa2_fxko < 0.05), ] 

plot(res$genotype_pct_pxko, res$pct_foxa22)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$pct_pax6_pxko, breaks = 40)

res = res[which(res$pct_pax6_pxko < 0.15), ]


hist(res$cutoff_wt, breaks = 40)
abline(v = c(0.1, 1.5), lwd = 2.0, col = 'red')
res = res[which(res$cutoff_wt > 0.1 & res$cutoff_foxa2 < 1.5), ]

plot(res$genotype_pct_pxko, res$pct_foxa22)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_ko, breaks = 20)
abline(v = c(0.2, 0.9), lwd = 2.0, col = 'red')

#res = res[which(res$cutoff_ko < 0.9), ]


plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

# Add regression lines
# ggplot(res, aes(x=genotype_pct_pxko, y=pct_foxa2, color=time, shape=time)) +
#   geom_point() + 
#   geom_smooth(method=lm, aes(fill=time))+
#   ylab("% FoxA2+ ") + 
#   xlab("genotype % Pax6-/-") + 
#   theme_bw() +  
#   theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
#         axis.text.y = element_text(angle = 0, size = 12)) +
#   theme(legend.key = element_blank()) + 
#   theme(plot.margin=unit(c(1,3,1,1),"cm"))+
#   #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
#   theme(legend.title = element_blank(), 
#         legend.text = element_text(size = 14))
# 
# ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day3_4_6_FoxA2pct.pdf'), height = 6, width = 10)
# 
# res = res[which(res$time == 'd4'), ]
ggplot(res, aes(x=genotype_pct_pxko, y=pct_foxa22)) +
  geom_point() + 
  geom_smooth(method=lm) +
  ylab("% FoxA2+ ") + 
  xlab("genotype % Pax6-/-") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day4_FoxA2pct_final.pdf'), height = 6, width = 10)


##########################################
# Pax6KO and FoxA2KO quantifying the genotypes and marker genes  
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "PKO_FKO_d6_stain_DVpatterning/",
                             "cyst_size_genotype_FoxA2_markers_genotype_FoxA_ostu_markerMultiotsu_manualOlig2global.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})
res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-2umZ','', res$time)
res$marker = sapply(res$image, function(x){unlist(strsplit(x, '_'))[6]})
res$marker = gsub('FA2-', '', res$marker)

res$genotype_pct_cyst = res$nb_pxko_cyst/(res$nb_pxko_cyst + res$nb_fxko_cyst)
res$genotype_pct = res$nb_pxko_global/(res$nb_pxko_global + res$nb_fxko_global)
#res$pct_foxa2 = (res$nb_foxa2_fxko + res$nb_foxa2_pxko)/res$genotype_total
#res$pct_marker = (res$nb_marker_fxko + res$nb_marker_pxko)/res$genotype_total

#res$pct_foxa2_cyst = res$nb_foxa2_cyst/res$cyst_size
res$pct_marker_cyst = res$nb_pax6_cyst/res$cyst_size
res$pct_marker = res$nb_pax6_global/res$cyst_size

table(res$treatment, res$marker)
#res$pct_foxa2_global = res$nb_foxa2_global/res$cyst_size
#res$pct_marker_global = res$/res$cyst_size



## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(3.5, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5.5), ]

table(res$treatment, res$marker)

saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))

## filtering for genotype quantification
res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))

plot(res$cutoff_fxko, res$cutoff_pxko)
kk = which(res$treatment == 'noRA')
points(res$cutoff_fxko[kk], res$cutoff_pxko[kk], col = 'red')
abline(v =c(3), col = 'red')
abline(h =c(6), col = 'red')

hist(res$cutoff_fxko, breaks = 100)
abline(v =c(2.5), col = 'red')

hist(res$cutoff_pxko, breaks = 200)
abline(v =c(6), col = 'red')

jj = which(res$cutoff_fxko < 3 &  res$cutoff_pxko <6)

res = res[jj, ]

saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_filtering_size_genotype.rds'))

### consider separate marker Olig2
res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_filtering_size_genotype.rds'))

res = res[which(res$marker == "Og2"), ]

aa = res[which(res$marker == "Og2" & res$treatment == 'RA'), ]

ggplot(res, aes(x=genotype_pct_cyst, y=pct_marker, color = treatment)) +
  geom_point() + 
  #geom_smooth(method=lm) +
  ylab("% Marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  geom_hline(yintercept=0.03, linetype="dashed", color = "red") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

hist(res$pct_marker[which(res$treatment == 'noRA')], breaks = 50)

df = c()
cutoffs = seq(0, 1, by = 0.2)
for(n in 1:(length(cutoffs)-1))
{
  index_group = which(res$treatment == 'RA'& res$genotype_pct_cyst >=cutoffs[n] & res$genotype_pct_cyst < cutoffs[n+1])
  df = c(df, length(which(res$pct_marker[index_group] > 0.03))/length(index_group))
}

df = data.frame(group = c('WT',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'), 
                pct_cystOlig2 = c(0.6733333, df))

df$group = factor(df$group, levels = c('WT',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'))

ggplot(data=df, aes(x=group, y=pct_cystOlig2, fill=group)) +
  geom_bar(stat="identity")+
  #theme_minimal() +
  ylab("% cyst with Olig2") + 
  xlab("% genotype Pax6-/-") +
  ylim(0, 0.7) + 
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") + 
  theme_bw() +  
  scale_fill_manual(values=c("darkgreen", rep("steelblue", 5))) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_KOKO_day6_DVpatterning_Olig2_vsWT.pdf'), height = 6, width = 8)


####### quantifying the Olig2 in WT with RA
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "WT_d6_stain_DVpatterning/",
                             "WT_d6_cyst_size_manualOlig2global_v2.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})
#res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
#res$time = gsub('-2umZ','', res$time)
#res$marker = sapply(res$image, function(x){unlist(strsplit(x, '_'))[6]})
res$pct_marker_cyst = res$nb_marker_cyst/res$cyst_size
res$pct_marker = res$nb_marker_global/res$cyst_size

plot(res$cyst_size, res$pct_marker)

## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(4.0, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^4.0 & res$cyst_size < 10^5.), ]

jj = which(res$treatment == 'NoRA')
hist(res$pct_marker[jj], breaks = 100)


jj2 = which(res$treatment == 'RA')
hist(res$pct_marker)
jj2 = which(res$treatment == 'RA')
plot(log10(res$cyst_size[jj2]), res$pct_marker[jj2])


length(which(res$pct_marker[jj2] > 0.005))/length(jj2)



### marker NKX22
res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_filtering_size_genotype_otsu.rds'))

res = res[which(res$marker == "N22"), ]

jj2 = which(res$treatment == 'RA')
jj3 = which(res$treatment == 'noRA')

plot(res$genotype_pct_cyst, res$pct_marker_cyst)
points(res$genotype_pct_cyst[jj3], res$pct_marker_cyst[jj3], pch = 16, col = 'red')
abline(h =c(0.15), col = 'red')


### Marker Nkx61
res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_filtering_size_genotype_otsu.rds'))

res = res[which(res$marker == "N61"), ]

jj2 = which(res$treatment == 'RA')
jj3 = which(res$treatment == 'noRA')

plot(res$genotype_pct_cyst, res$pct_marker_cyst)
points(res$genotype_pct_cyst[jj3], res$pct_marker_cyst[jj3], pch = 16, col = 'red')
abline(h =c(0.15), col = 'red')


aa = res[which(res$treatment == 'RA'), ]
## plot the pct FoxA2 in function of genotype pct
ggplot(aa, aes(x=genotype_pct_cyst, y=pct_marker_cyst, color = marker)) +
  geom_point() + 
  geom_smooth(method=lm) +
  ylab("% Marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

aa = res[which(res$marker == "Og2" & res$treatment == 'RA'), ]

ggplot(aa, aes(x=genotype_pct_cyst, y=pct_marker, color = marker)) +
  geom_point() + 
  #geom_smooth(method=lm) +
  ylab("% Marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

aa$cutoff_pax6_image[match(unique(aa$image[which(aa$marker == "Og2")]), aa$image)]

ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day6_DVpatterning_FoxA2pct.pdf'), height = 6, width = 10)

#saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))
#res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))
#res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_all_global.rds'))
res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_all_nomarkerNormalization.rds'))


hist(log10(res$cyst_size), breaks = 50)
abline(v = c(4, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^4 & res$cyst_size < 10^5.5), ]

plot(res$pct_foxa2_global, res$pct_foxa2_cyst)

jj = which(res$marker == 'Og2')
df = data.frame(cutoff = res$cutoff_marker[jj], quantile = res$quantile95_marker[jj], condition = res$treatment[jj])
df$condition = factor(df$condition, levels = c('RA', 'noRA'))

ggplot(df, aes(x=cutoff, fill=condition)) +
  geom_histogram(alpha = 0.7, position = 'identity') + 
  #geom_histogram(alpha=0.7, adjust = 1.5) +
  scale_fill_manual(values=c("darkgreen", "#999999")) + 
  xlab("% FoxA2+ ") + 
  ylab("Density") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


ggplot(df, aes(x=quantile, fill=condition)) +
  geom_histogram(alpha = 0.7, position = 'identity') + 
  #geom_histogram(alpha=0.7, adjust = 1.5) +
  scale_fill_manual(values=c("darkgreen", "#999999")) + 
  xlab("% FoxA2+ ") + 
  ylab("Density") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

abline(v = c(1.5, 2.5), lwd = 2.0, col = 'red')


kk1 = which(res$marker == 'Og2' & res$treatment == 'RA' & res$cutoff_marker < 1.5)
res$pct_marker[kk1] = 0

hist(res$cutoff_marker[which(res$marker == 'N22' & res$treatment == "noRA")], breaks = 30)
abline(v = c(1.5, 2.5), lwd = 2.0, col = 'red')

kk2 = which(res$marker == 'N22' & res$treatment == 'RA' & res$cutoff_marker < 1.0)
res$pct_marker[kk2] = 0

hist(res$cutoff_marker[which(res$marker == 'N61' & res$treatment == "noRA")], breaks = 30)
abline(v = c(1.5, 2.5), lwd = 2.0, col = 'red')

kk2 = which(res$marker == 'N22' & res$treatment == 'RA' & res$cutoff_marker < 1.0)
res$pct_marker[kk2] = 0

res = res[which(res$treatment == 'RA'), ]

#jj = which(res$marker == 'Og2')
#plot(res$cutoff_marker[jj], res$pct_marker[jj])


#jj = which(res$marker == 'Og2' & res$pct_marker >0.5 & res$genotype_pct_pxko <0.5)
jj = which(res$marker == 'Og2')
aa = res[jj, ]

# Add regression lines
ggplot(aa, aes(x=genotype_pct_pxko, y=pct_marker_cyst, color=treatment, shape=marker)) +
  geom_point(size = 2.5) + 
  #geom_smooth(method=lm, aes(fill=marker))+
  ylab("% marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  #xlim(0, 2) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))
#res = res[which(res$pct_marker >0), ]


# Add regression lines
ggplot(aa, aes(x=genotype_pct_pxko, y=pct_marker, color=treatment, shape=marker)) +
  geom_point(size = 2.5) + 
  #geom_smooth(method=lm, aes(fill=marker))+
  ylab("% marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  #xlim(0, 2) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'PXKO_FXKO_day6_Olig2pct_genotype_RA.noRA.pdf'), height = 6, width = 10)

jj = which(res$marker == 'Og2' & res$treatment == 'RA')
aa = res[jj, ]

ggplot(aa, aes(x=genotype_pct_pxko, y=pct_marker_cyst, color = cutoff_marker,   shape=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=lm, aes(fill=marker))+
  ylab("% marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  #xlim(0, 2) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

jj = which(res$marker == 'Og2' & res$treatment == 'RA' & res$cutoff_marker > 170)
aa = res[jj, ]


ggplot(aa, aes(x=genotype_pct_pxko, y=pct_marker,  shape=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=lm, aes(fill=marker))+
  ylab("% marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  #xlim(0, 2) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'PXKO_FXKO_day6_Olig2pct_genotype_RA.pdf'), height = 6, width = 10)


########################################################
########################################################
# Section III TetOn-FoxA2-TetOnPax6:
# genotype and FoxA2
########################################################
########################################################

##########################################
# genotype and FoxA2 at d4
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "TetOnF_TetOnP_chim_d3_4_5_genotype_FoxA2_Pax6/",
              "cyst_size_genotype_FoxA2_Pax6_cystThresholds_globalThreshold_quantile_normalization_refineGenotype_thresholdMean_v2.csv"), 
               header = TRUE, row.names = c(1))

# res2 = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
#                               "TetOnF_TetOnP_chim_d3_4_5_genotype_FoxA2_Pax6/",
#                               "cyst_size_genotype_FoxA2_Pax6_cystThresholds_globalThreshold_quantile.csv"), 
#                 header = TRUE, row.names = c(1))
# 
# rownames(res) = paste0(res$image, '_', res$cyst_index) 
# rownames(res2) = paste0(res2$image, '_', res2$cyst_index) 
# res2 = res2[match(rownames(res), rownames(res2)), ]

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})

res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-P6-FA2','', res$time)

res$genotype_pct_global = res$nb_fxko_global/(res$nb_fxko_global + res$nb_pxko_global)
res$pct_foxa2_global = (res$nb_foxa2_global)/(res$nb_foxa2_global + res$nb_pax6_global)

#res$genotype_pct_global2 = res$nb_fxko_global/res$cyst_size
#res$pct_foxa2_global2 = (res$nb_foxa2_global)/res$cyst_size

#res$genotype_pct_cyst2 = res$nb_fxko_cyst/res$cyst_size
#res$pct_foxa2_cyst2 = res$nb_foxa2_cyst/res$cyst_size

res$genotype_pct_cyst = res$nb_fxko_cyst/(res$nb_fxko_cyst + res$nb_pxko_cyst)
#res$genotype_pct_cyst2 = res$nb_fxko/(res$nb_fxko + res$nb_pxko) 
res$pct_foxa2_cyst = res$nb_foxa2_cyst/(res$nb_foxa2_cyst + res$nb_pax6_cyst + res$nb_double_cyst)
res$pct_foxa2_cyst2 = (res$fxko_nb_foxa2 + res$pxko_nb_foxa2)/(res$fxko_nb_foxa2 + res$pxko_nb_foxa2 + 
                                                                 res$fxko_nb_pax6 + res$pxko_nb_pax6)

res$pct_foxa2_cyst3 = res$pct_foxa2_cyst
kk = which(res$treatment == 'dox')

xx1 = res$pxko_nb_foxa2[kk]/(res$pxko_nb_foxa2[kk] + res$pxko_nb_pax6[kk])
hist(xx1, breaks = 100)

xx2 = res$fxko_nb_foxa2[kk]/(res$fxko_nb_foxa2[kk] + res$fxko_nb_pax6[kk])
hist(xx2, breaks = 100)

res$pct_foxa2_cyst3[kk] = (res$fxko_nb_foxa2[kk])/(res$fxko_nb_foxa2[kk] + res$pxko_nb_pax6[kk])
#res$pct_foxa2_cyst3[kk] = (res$fxko_nb_foxa2[kk])/(res$fxko_nb_foxa2[kk] + res$fxko_nb_pax6[kk] + res$pxko_nb_pax6[kk])
#res$pct_foxa2_cyst3 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/res$genotype_total*2
#res$genotype_pct_pxko = res$nb_fxko/res$genotype_total
#res$pct_foxa2 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/res$genotype_total*2


## select only the d4 and discard RAdox
# res = res[which(res$time == 'd4' & res$treatment != "RAdox"), ]
# rownames(res) = paste0(res$image, '_', res$cyst_index) 
# 
# saveRDS(res, file = paste0(outDir, 'TetOn_TetOn_genotype_FoxA2pct_Threshold_otsu.rds'))
# 
# res = readRDS(file = paste0(outDir, 'TetOn_TetOn_genotype_FoxA2pct_Threshold_otsu.rds'))
# res = res[which(res$treatment != 'noRA'), ]
# res2 = readRDS(file = paste0(outDir, 'TetOn_TetOn_genotype_FoxA2pct_Threshold_mean.rds'))
# res2 = res2[match(rownames(res), rownames(res2)), ]
# 
# plot(res$genotype_pct_cyst, res2$genotype_pct_cyst);
# abline(0, 1, lwd = 2.0, col = 'red')
# kk = which(res$treatment == 'dox')
# points(res$genotype_pct_cyst[kk], res2$genotype_pct_cyst[kk], col = 'darkblue')
# 
# plot(res$pct_foxa2_cyst3, res2$pct_foxa2_cyst3);
# abline(0, 1, lwd = 2.0, col = 'red')
# kk = which(res$treatment == 'dox')
# points(res$pct_foxa2_cyst3[kk], res2$pct_foxa2_cyst3[kk], col = 'darkblue')
# 
# kk = which(res$treatment == 'dox')
# 
# res$genotype_pct_cyst[kk] = res2$genotype_pct_cyst[kk]
# res$pct_foxa2_cyst3[kk] = res2$pct_foxa2_cyst3[kk]

plot(res$genotype_pct_cyst, res$genotype_pct_global);abline(0, 1, col = 'red', lwd = 2.0)

plot(res$pct_foxa2_cyst, res$pct_foxa2_global);abline(0, 1, col = 'red', lwd = 2.0)

plot(res$pct_foxa2_cyst, res$pct_foxa2_cyst2);abline(0, 1, col = 'red', lwd = 2.0)

plot(res$pct_foxa2_cyst2, res$pct_foxa2_cyst3);abline(0, 1, col = 'red', lwd = 2.0)

plot(res$pct_foxa2_cyst, res$pct_foxa2_cyst3);abline(0, 1, col = 'red', lwd = 2.0)

#plot(res$genotype_pct_cyst, res$genotype_pct_cyst2);
#abline(0, 1, lwd = 2.0, col = 'red')
aa = res[which(res$treatment != 'noRA'), ]
p1 = ggplot(aa, aes(x=genotype_pct_cyst, y=pct_foxa2_cyst3, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=treatment)) + 
  geom_hline(yintercept=c(0.3), linetype="solid", linewidth = 1) +  
  ylab("% FoxA2+ ") + 
  xlab("genotype % TetOn-FoxA2") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


########################
## filtering steps 
## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(3.5, 5), col = 'red', lwd = 2.0)

hist(log10(res$cyst_size[which(res$treatment == 'dox')]), breaks = 100)

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5), ]

table(res$treatment)

aa = res[which(res$treatment != 'noRA'), ]
p2 = ggplot(aa, aes(x=genotype_pct_cyst, y=pct_foxa2_cyst3, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=treatment))+ 
  geom_hline(yintercept=c(0.3), linetype="solid", linewidth = 1) +  
  ylab("% FoxA2+ ") + 
  xlab("genotype % TetOn-FoxA2") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

p1 /p2

kk1 = which(res$treatment == 'dox' & res$pct_foxa2_cyst3 > 0.75 & res$genotype_pct_cyst > 0.8)

plot(res$cutoff_fxko, res$cutoff_pxko, cex = 0.7, xlim = c(0, 5), ylim = c(0, 4))
abline(v = c(0.4, 5), col = 'red')
abline(h = 0.2, col = 'red')
kk = which(res$treatment == 'dox')
points(res$cutoff_fxko[kk], res$cutoff_pxko[kk], col = 'blue', cex = 0.7)
kk = which(res$treatment == 'noRA')
points(res$cutoff_fxko[kk], res$cutoff_pxko[kk], col = 'red', cex = 0.7, pch =4)
points(res$cutoff_fxko[kk1], res$cutoff_pxko[kk1], col = 'orange', cex = 1.5, pch =16)

#hist(res$cutoff_fxko, breaks = 100);abline(v = c(0.5, 5))

kk1 = which(res$treatment == 'dox')
plot(res$cutoff_fxko[kk1], res$cutoff_pxko[kk1], cex = 0.7, xlim = c(0, 5), ylim = c(0, 4))

jj1 = kk1[which(res$cutoff_fxko[kk1] < 4.5 & res$cutoff_pxko[kk1] <2.5)]

kk2 = which(res$treatment == 'RA')
plot(res$cutoff_fxko[kk2], res$cutoff_pxko[kk2], cex = 0.7, xlim = c(0, 5), ylim = c(0, 5))

jj2 = kk2[which(res$cutoff_fxko[kk2] < 4 & res$cutoff_pxko[kk2] <4.5)]

res = res[c(jj1, jj2), ]

#hist(res$cutoff_fxko[)], breaks = 100);abline(v = c(0.5, 5))
#hist(res$cutoff_fxko[which(res$treatment == 'RA')], breaks = 100);abline(v = c(0.5, 5))
#hist(res$cutoff_pxko, breaks = 100);abline(v = c(0.4, 4))

#jj1 = which(res$cutoff_fxko > 0.5  & res$cutoff_fxko < 5 & res$cutoff_pxko > 0.4 &res$cutoff_pxko < 4 )
#res = res[jj1, ]
table(res$treatment)

aa = res[which(res$treatment != 'noRA'), ]
ggplot(aa, aes(x=genotype_pct_cyst, y=pct_foxa2_cyst3, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=treatment))+
  ylab("% FoxA2+ ") + 
  xlab("genotype % TetOn-FoxA2") + 
  theme_bw() +  
  geom_hline(yintercept=c(0.3), linetype="solid", linewidth = 1) +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


plot(res$cutoff_foxa2, res$cutoff_pax6, ylim = c(0, 3), xlim = c(0, 6),  type = 'n')
#kk = which(res$treatment == 'noRA')
#points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'red', pch = 2)
kk = which(res$treatment == 'dox')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'blue', pch = 16)

kk = which(res$treatment == 'RA')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'orange', pch = 4)
abline(v = c(4.0, 1.0), col = 'red')
abline(h = c(0.3, 2.0), col = 'red')


kk1 = which(res$treatment == 'dox')
plot(res$cutoff_foxa2[kk1], res$cutoff_pax6[kk1], col = 'blue', pch = 16)
abline(v = c(0.7, 3.5))

ii1 = kk1[which(res$cutoff_foxa2[kk1]< 0.7)]
points(res$cutoff_foxa2[ii1], res$cutoff_pax6[ii1], col = 'red', pch = 4)

plot(res$genotype_pct_cyst[kk1], res$pct_foxa2_cyst3[kk1])
points(res$genotype_pct_cyst[ii1], res$pct_foxa2_cyst3[ii1], col = 'red', pch = 4)

jj1 = kk1[which(res$cutoff_foxa2[kk1]> 0.7 & res$cutoff_foxa2[kk1] < 3.5)]

kk2 = which(res$treatment == 'RA')

plot(res$cutoff_foxa2[kk2], res$cutoff_pax6[kk2], col = 'blue', pch = 16)
abline(v = c(1., 4))

ii2 = kk2[which(res$cutoff_pax6[kk2] > 1.7)]
points(res$cutoff_foxa2[ii2], res$cutoff_pax6[ii2], col = 'red', pch = 4)

plot(res$genotype_pct_cyst[kk2], res$pct_foxa2_cyst3[kk2])
points(res$genotype_pct_cyst[ii2], res$pct_foxa2_cyst3[ii2], col = 'red', pch = 4)
abline(h = 0.3, lwd = 2.0, col = 'red')

jj2 = kk2[which(res$cutoff_foxa2[kk2]> 1.0 & res$cutoff_foxa2[kk2] < 3.5 & res$cutoff_pax6[kk2] < 2.0)]


res = res[c(jj1, jj2), ]

# hist(res$cutoff_foxa2[which(res$treatment == 'RA')], breaks = 40)
# abline(v = c(1000, 4000))
# 
# hist(res$cutoff_foxa2[which(res$treatment == 'dox')], breaks = 40)
# abline(v = c(1000, 6000))
# 
# hist(res$cutoff_pax6[which(res$treatment == 'dox')], breaks = 40)
# abline(v = c(1000, 6000))
# 
# #hist(res$cutoff_pax6[which(res$treatment == 'RAdox')], breaks = 40)
# #abline(v = c(450, 1500))
# 
# 
# j1 = which(res$cutoff_foxa2 > 1000 & res$cutoff_foxa2 <4000 & res$treatment == 'RA')
# j1 = which(res$cutoff_foxa2 <3000 & res$treatment == 'RA')
# 
# #j2 = which(res$cutoff_pax6 > 450 & res$cutoff_pax6 < 1500 & res$treatment == 'RAdox')
# j3 = which(res$cutoff_foxa2 > 1000 & res$cutoff_foxa2 <5000 & res$cutoff_pax6 > 450 & res$treatment == 'dox')
# 
# jj2 = c(j1, j3)
# 
# hist(res$cutoff_foxa2[which(res$treatment != 'noRA')], breaks = 40)
# abline(v = c(1500, 5000))
# 
# hist(res$cutoff_pax6[which(res$treatment != 'noRA')], breaks = 40)
# abline(v = c(400, 1200))

# jj = which(res$cutoff_foxa2 < 4.0 & res$cutoff_foxa2 > 1.0 &
#            res$cutoff_pax6 < 2.0 & res$cutoff_pax6 > 0.3)
# 
# res =res[jj, ]

aa = res[which(res$treatment != 'noRA'), ]

ggplot(aa, aes(x=genotype_pct_cyst, y=pct_foxa2_cyst3, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=treatment))+
  ylab("% FoxA2+ ") + 
  xlab("genotype % TetOn-FoxA2") + 
  geom_hline(yintercept=c(0.3), linetype="solid", linewidth = 1) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'TetOnFx_TetOnPx_day4_FoxA2pct_genotypepct_thresholdMean_v3.pdf'), 
       height = 6, width = 10)


##########################################
# d6 stain patterning (markers SHH, NKX22) 
##########################################

## test manually specify global throsholds
aa = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                            "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning/",
                            "image_gloablThresholds_multiotsu_Nkx_Shh.csv"), 
              header = TRUE, row.names = c(1))

aa = aa[, c(1, 3, 5, 6, 8, 10, 11)]

cutoff_nkx = c()
for(n in 1:nrow(aa))
{
  cuts = as.numeric(aa[n, c(2:4)])
  diff = abs(cuts - 1.0)
  cutoff_nkx = c(cutoff_nkx, cuts[which.min(diff)])
}

aa = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                            "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning/",
                            "image_gloablThresholds_multiotsu_Nkx_ShhfilteringLargeValues.csv"), 
              header = TRUE, row.names = c(1))

aa = aa[, c(1,  8, 10, 11)]

cutoff_shh = c()
for(n in 1:nrow(aa))
{
  cuts = as.numeric(aa[n, c(2:4)])
  diff = abs(cuts - 1.25)
  cutoff_shh = c(cutoff_shh, cuts[which.min(diff)])
}




res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning/",
                             "cyst_size_genotype_cyst_Nkx22_Shh_image_manualThresholds.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})

res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-Nk22-Shh','', res$time)

## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(3.7, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^4.0 & res$cyst_size < 10^5.5), ]

#res = res[which(res$treatment != 'noRA'), ]
### check the genotype
saveRDS(res, file = paste0(outDir, 'd6_TetOn_TetOn_Shh_Nkx22_sizeFiltering.rds'))


### check the SHH
res = readRDS(file = paste0(outDir, 'd6_TetOn_TetOn_Shh_Nkx22_sizeFiltering.rds'))

#res$genotype_pct_global = res$nb_fxko_global_otsu/(res$nb_fxko_global_otsu + res$nb_pxko_global_otsu)
res$genotype_pct_cyst = res$nb_fxko_cyst_otsu/(res$nb_fxko_cyst_otsu + res$nb_pxko_cyst_otsu)
#res$genotype_pct_cyst = res$nb_fxko_cyst_mean/(res$nb_fxko_cyst_mean + res$nb_pxko_cyst_mean)

#res$genotype_pct_cyst = res$nb_fxko_cyst_li/(res$nb_fxko_cyst_li + res$nb_pxko_cyst_li)
#res$genotype_pct_cyst = res$nb_fxko_cyst_yen/(res$nb_fxko_cyst_yen + res$nb_pxko_cyst_yen)
#res$genotype_pct_cyst = res$nb_fxko_cyst_isodata/(res$nb_fxko_cyst_isodata + res$nb_pxko_cyst_isodata)


#res$pct_nkx22_cyst = res$nb_foxa2_cyst/res$cyst_size
#res$pct_shh_global = res$nb_shh_global_multiotsu/res$cyst_size
#res$pct_shh_cyst = res$nb_shh_cyst_otsu/res$cyst_size

res$pct_shh_global = res$nb_shh_global_multiotsu/res$cyst_size
res$pct_nkx22_global = res$nb_nkx_global_multiotsu/res$cyst_size

res$pct_shh_global_2 = res$nb_shh_fxko/res$nb_total_genotype
res$pct_nkx22_global_2 = res$nb_nkx_fxko/res$nb_total_genotype

#res$pct_shh_global = res$nb_shh_global_isodata/res$cyst_size
#res$pct_nkx22_global = res$nb_nkx_global_isodata/res$cyst_size

## filtering the genotype pct
#jj1 = which(res$treatment != 'noRA')
#res = res[jj1, ]

plot((res$cutoff_fxko_cyst_otsu/res$cutoff_fxko_global_otsu), 
     (res$cutoff_pxko_cyst_otsu/res$cutoff_pxko_global_otsu))
jj1 = which(res$treatment == 'noRA')
points(res$cutoff_fxko_cyst_otsu[jj1]/res$cutoff_fxko_global_otsu[jj1], 
       res$cutoff_pxko_cyst_otsu[jj1]/res$cutoff_pxko_global_otsu[jj1], col = 'red', pch = 2)

abline(v = 2.1)
abline(h = 1.7)

res$outliers = FALSE
res$outliers[which(res$cutoff_fxko_cyst_otsu/res$cutoff_fxko_global_otsu > 2.1 |
                     res$cutoff_pxko_cyst_otsu/res$cutoff_pxko_global_otsu > 1.7)] = TRUE

res = res[!res$outliers, ]

plot(res$cutoff_fxko_cyst_otsu, res$cutoff_pxko_cyst_otsu)
jj1 = which(res$treatment == 'dox')
points(res$cutoff_fxko_cyst_otsu[jj1], res$cutoff_pxko_cyst_otsu[jj1], col = 'red', pch = 2)

jj1 = which(res$treatment == 'RA')
points(res$cutoff_fxko_cyst_otsu[jj1], res$cutoff_pxko_cyst_otsu[jj1], col = 'blue', pch = 2)


saveRDS(res, file = paste0(outDir, 'd6_TetOn_TetOn_Shh_Nkx22_Filtering_size_genotype.rds'))

# jj1 = which(res$treatment == 'dox')
# res = res[jj1, ]
# res$pct_shh_false = res$nb_shh_pxko/res$nb_pxko_cyst_otsu
# res$pct_nkx_false = res$nb_nkx_pxko/res$nb_pxko_cyst_otsu
# 
# plot(res$genotype_pct_cyst, res$pct_shh_false)
# 
# res$outliers = FALSE
# res$outliers[which(res$cutoff_fxko_cyst_otsu/res$cutoff_fxko_global_otsu > 2.1 |
#                      res$cutoff_pxko_cyst_otsu/res$cutoff_pxko_global_otsu > 1.7)] = TRUE
# 
# res = res[!res$outliers, ]



### check the SHH
res = readRDS(file = paste0(outDir, 'd6_TetOn_TetOn_Shh_Nkx22_Filtering_size_genotype.rds'))

ggplot(res, aes(x=genotype_pct_cyst, y=pct_shh_global, color = treatment, fill=treatment)) +
  geom_point() + 
  #geom_smooth(method=lm) +
  ylab("% Marker+ ") + 
  xlab("genotype % Pax6-/-") + 
  geom_hline(yintercept=0.03, linetype="dashed", color = "red") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

hist(res$pct_shh_global[which(res$treatment == 'noRA')], breaks = 50)

df = c()
cutoffs = seq(0, 1, by = 0.2)

for(n in 1:(length(cutoffs)-1))
{
  index_group = which(res$treatment == 'dox'& 
                        res$genotype_pct_cyst >=cutoffs[n] & res$genotype_pct_cyst < cutoffs[n+1])
  df = c(df, length(which(res$pct_shh_global[index_group] > 0.05))/length(index_group))
}

index_wt = which(res$treatment == 'RA')
pct_wt = length(which(res$pct_shh_global[index_wt] > 0.05))/length(index_wt)

index_noRA = which(res$treatment == 'noRA')
pct_noRA = length(which(res$pct_shh_global[index_noRA] > 0.05))/length(index_noRA)

df = data.frame(group = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'), 
                pct_shh = c(pct_wt, df))

df$group = factor(df$group, levels = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'))

ggplot(data=df, aes(x=group, y=pct_shh, fill=group)) +
  geom_bar(stat="identity")+
  #theme_minimal() +
  ylab("% cyst with Shh") + 
  xlab("% genotype TetOn-FoxA2") +
  ylim(0, 0.4) + 
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") + 
  theme_bw() +  
  scale_fill_manual(values=c("darkgreen", rep("deepskyblue", 5))) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'TetOn_TetOn_day6_DVpatterning_Shh_vsWT.pdf'), height = 6, width = 8)



res = res[which(res$treatment != "noRA"), ]
p1 = ggplot(res, aes(x=genotype_pct_cyst, y = pct_shh_global, color = treatment, fill=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=loess, aes(fill=treatment))+
  ylab("% Shh+ ") + 
  xlab("genotype % FoxA2-TetOn") + 
  #ylim(0, 0.5) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

p2 = ggplot(res, aes(x=genotype_pct_cyst, y = pct_nkx22_global, color = treatment, fill=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=loess, aes(fill=treatment))+
  ylab("% Nkx2.2 ") + 
  xlab("genotype % FoxA2-TetOn") + 
  #xlim(0, 2) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

p1 / p2


########################################################
########################################################
# Section IV : WT_FoxA2 chimeras and WT_Pax6 chimeras
# Fig 4B
########################################################
########################################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "FoxAKO_WTchim_d4/",
                             "cyst_size_genotype_FoxA2_Pax6_cystThresholds_test_otsuThresholdGenotype_meanThresholdFoxA2_originalGenetype_Normallization_v6.csv"), 
               header = TRUE, row.names = c(1))


res$pct_foxa2 = res$nb_foxa2/(res$nb_foxa2 + res$nb_pax6)

plot(res$pct_ko, res$pct_foxa2)


## size filtering 
hist(log10(res$cyst_size), breaks = 50)

res = res[which(res$cyst_size > 10^4), ]
plot(res$pct_ko, res$pct_foxa2)

plot(res$cutoff_wt, res$cutoff_ko)
abline(v = 17)
abline(h = 3.5)


res = res[which(res$cutoff_wt < 17 & res$cutoff_ko <3.5), ]

plot(res$pct_ko, res$pct_foxa2)

plot(res$cutoff_foxa2, res$cutoff_pax6)
abline(v = 17)
abline(h = 3.5)

plot(res$pct_ko, res$pct_foxa2)
jj1 = which(res$cutoff_foxa2 > 6)
points(res$pct_ko[jj1], res$pct_foxa2[jj1], pch = 16, col = 'red')

res = res[which(res$cutoff_foxa2 < 6), ]

plot(res$pct_ko, (res$nb_foxa2 + res$nb_pax6))


ggplot(res, aes(x=pct_ko, y=pct_foxa2)) +
  geom_point() + 
  geom_smooth(method=loess) +
  ylab("% FoxA2+ in WT ") + 
  xlab("genotype % FoxA2-/-") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'FoxA2KO_WTchimeras_day4.pdf'), height = 6, width = 10)


##########################################
# Pax6KO-WT chimeras
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "Pax6KO_WTchim_d4/",
                             "cyst_size_genotype_FoxA2_Pax6_cystThresholds_testMeanThreshodling_v2.csv"), 
               header = TRUE, row.names = c(1))

Compare_compressed_version = FALSE
if(Compare_compressed_version){
  
  xx = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                              "Pax6KO_WTchim_d4/",
                              "cyst_size_genotype_FoxA2_Pax6_cystThresholds_testMeanThreshodling_test_jetraw_compression.csv"), 
                header = TRUE, row.names = c(1))
  
  rownames(res) = paste0(res$image, "_", res$cyst_index)
  rownames(xx) = paste0(xx$image, "_", xx$cyst_index)
  
  xx = xx[match(rownames(res), rownames(xx)), ]
  
  plot(res$cyst_size, xx$cyst_size, xlab = 'original', ylab = 'compressed', main = 'cyst size')
  abline(0, 1, col = 'red', lwd =2.0)
  
  
  plot(res$pct_wt, xx$pct_wt, xlab = 'original', ylab = 'compressed', main = '% of wt genotype')
  abline(0, 1, col = 'red', lwd =2.0)
  
  plot(res$pct_foxa2, xx$pct_foxa2, xlab = 'original', ylab = 'compressed', main = '% of FoxA2+')
  abline(0, 1, col = 'red', lwd =2.0)
  
  #res = xx
    
}


#res$pct_foxa2 = res$nb_foxa2/(res$nb_foxa2 + res$nb_pax6)
res$pct_ko = 1.0- res$pct_wt 
plot(res$pct_ko, res$pct_foxa2)


## size filtering 
hist(log10(res$cyst_size), breaks = 50)

res = res[which(res$cyst_size > 10^4), ]
plot(res$pct_ko, res$pct_foxa2)

plot(res$cutoff_wt, res$cutoff_ko, log = 'xy')
abline(v = 20)
abline(h = 4)


res = res[which(res$cutoff_wt < 20 & res$cutoff_ko < 4 & res$cutoff_wt > 1.5), ]

plot(res$pct_ko, res$pct_foxa2)

plot(res$cutoff_foxa2, res$cutoff_pax6, log = '')
abline(v = 17)
abline(h = 3.5)

plot(res$pct_ko, res$pct_foxa2)
jj1 = which(res$cutoff_foxa2 > 6)
points(res$pct_ko[jj1], res$pct_foxa2[jj1], pch = 16, col = 'red')

#res = res[which(res$cutoff_foxa2 < 4 & res$cutoff_foxa2 > 0.5 & res$cutoff_pax6 > 1 & res$cutoff_pax6 <3.5), ]

#plot(res$pct_ko, (res$nb_foxa2 + res$nb_pax6))

ggplot(res, aes(x=pct_ko, y=pct_foxa2)) +
  geom_point() + 
  geom_smooth(method=loess) +
  ylab("% FoxA2+ in WT ") + 
  xlab("genotype % Pax6-/-") + 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'PAX6KO_WTchimeras_day4_compressedVerion.pdf'), height = 6, width = 10)

