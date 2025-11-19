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

outDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                "results/figures_tables_R13547_10x_mNT_20240522/")


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


## here use the mean threshold for FoxA2 proportions
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
  
  res$pct_foxa2 = res$nb_foxa2_localThreshold/(res$nb_double_localThreshold + res$nb_foxa2_localThreshold +
                                                 res$nb_pax6_localThreshold)
  res$pct_pax6 = res$nb_pax6_localThreshold/(res$nb_double_localThreshold + res$nb_foxa2_localThreshold +
                                                 res$nb_pax6_localThreshold)
  
  res$pct_double = res$nb_double_mean_cyst/(res$nb_foxa2_mean_cyst + res$nb_pax6_mean_cyst + 
                                              res$nb_double_mean_cyst)
  
  
  res$pct_pax6 = res$nb_pax6_otsu_cyst/(res$nb_pax6_otsu_cyst + res$nb_foxa2_otsu_cyst + 
                                              res$nb_double_otsu_cyst)
  
  res$pct_pax6 = res$nb_pax6_otsu_cyst/(res$nb_pax6_otsu_cyst+ res$nb_foxa2_mean_cyst)
  hist(res$pct_pax6)
  abline(v = 0.6, col = 'blue')
  
  hist(res$pct_foxa2)
  abline(v = 0.3, col = 'blue')
  
  res$pct_pax6 = res$nb_pax6_otsu_cyst/(res$nb_pax6_otsu_cyst+ res$nb_foxa2_mean_cyst)
  
  plot(res$nb_pax6_mean_cyst, res$nb_)
  abline(0, 1, lwd =2.0, col = 'red')
  
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
  
  
  plot(res$cyst_size, res$pct_foxa2)
  
  plot(res$cyst_size, res$pct_double) ## here total is the number of FoxA2+, Pax6+ and FoxA2+&Pax6+ 
  abline(a = 0, b = 0.5, lwd = 2.0, col = 'red')
  
  plot(density(res$pct_foxa2, adjust = 1.5), col = 'darkgreen', lwd = 3.0,  xlim = c(0, 1), ylim=c(0, 6))
  abline(v = 0.28, lwd = 2.0)
  
  df = data.frame(pct = res$pct_foxa2, condition = rep('wt', nrow(res)))
  
  saveRDS(df, file = paste0(outDir, 'FoxA2_pct_wt.rds'))
  
  df = readRDS(file = paste0(outDir, 'FoxA2_pct_wt.rds'))
  
  ggplot(df, aes(x=pct, fill = condition)) +
    geom_density(alpha=0.7, adjust = 1.5) + 
    xlim(0, 1) + 
    scale_fill_manual(values=c("#337f01")) + 
    xlab("% FoxA2+ ") + 
    ylab("Density") + 
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  ggsave(filename = paste0(outDir, 'RA_WT_day4_FoxA2pct_thresholdMean.pdf'), height = 6, width = 6)
  
  
  res$pct_pax6 = 1.0 - res$pct_foxa2 - res$pct_double
  
  df = data.frame(pct = c(res$pct_foxa2, res$pct_double, res$pct_pax6), 
                  gene = c(rep('FoxA2+', nrow(res)), rep('double+', nrow(res)), rep('Pax6+', nrow(res)))
  )
  
  ggplot(df, aes(x=pct, fill = gene)) +
    geom_density(alpha=0.7, adjust = 1.5) + 
    xlim(0, 1) + 
    #scale_fill_manual(values=c("darkgreen")) + 
    xlab("% FoxA2+ ") + 
    ylab("Density") + 
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  ggsave(filename = paste0(outDir, 'RA_WT_day4_FoxA2pct_thresholdMean.pdf'), height = 6, width = 6)
  
  
}

Use_Ilastik_pixelProb = FALSE
if(Use_Ilastik_pixelProb)
{
  res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                               "images_data/results/test_WT_FoxA2_pct_ilastik/",
                               "wt_FoxA2_Pax6_voxel_counts_ilastikPixelProb_v2.csv"), 
                 header = TRUE, row.names = c(1))
  
  res$total = res$nb_pax6_cyst + res$nb_foxa2_cyst + res$nb_double_cyst
  res$pct_foxa2 = res$nb_foxa2_cyst/res$total
  res$pct_pax6 = res$nb_pax6_cyst/res$total
  res$pct_double = res$nb_double_cyst/res$total
  
  plot(res$cyst_size, res$pct_foxa2)
  abline(v = c(5000, 10000), lwd = 2.0, col = 'red')
  
  res = res[which(res$cyst_size > 5000), ]
  
  df = data.frame(pct = c(res$pct_foxa2, res$pct_double, res$pct_pax6), 
                  gene = c(rep('FoxA2+', nrow(res)), rep('double+', nrow(res)), rep('Pax6+', nrow(res)))
  )
  
  ggplot(df, aes(x=pct, fill = gene)) +
    geom_density(alpha=0.7, adjust = 1.5) + 
    xlim(0, 1) + 
    #scale_fill_manual(values=c("darkgreen")) + 
    xlab("% FoxA2+ ") + 
    ylab("Density") + 
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14))
  
  
  ggsave(filename = paste0(outDir, 'RA_WTd4_pct_FoxA2_Pax6_double_v2.pdf'), height = 6, width = 8)
  
  
}


##########################################
# simulation of FoxA2 proportion 
##########################################
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

saveRDS(df, file = paste0(outDir, ))

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
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day4_FoxA2pct_final.pdf'), height = 6, width = 8)


## compare the FoxA2 percentages with the WT 
df = readRDS(file = paste0(outDir, 'FoxA2_pct_wt.rds'))

df = rbind(df, data.frame(pct = res$pct_foxa22, condition = rep('KO_KO', nrow(res))))

ggplot(df, aes(x=pct, fill=condition)) +
  geom_density(alpha=0.5, adjust = 1.5) +
  scale_fill_manual(values=c("#007BB7", "#337f01")) + 
  xlab("% FoxA2+ ") + 
  ylab("Density") + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14)) +
  ggtitle('KS-test : p < 1e-10')

ggsave(filename = paste0(outDir, 'RA_day4_FoxA2pct_compare_KOKOvsWT.pdf'), height = 4, width = 8)

ks.test(x = df$pct[which(df$condition == 'wt')], 
        y = df$pct[which(df$condition == 'KO_KO')])

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
pvals = c()

for(n in 1:(length(cutoffs)-1))
{
  index_group = which(res$treatment == 'RA'& res$genotype_pct_cyst >=cutoffs[n] & res$genotype_pct_cyst < cutoffs[n+1])
  
  p = binom.test(x = length(which(res$pct_marker[index_group] > 0.03)), 
             n = length(index_group), 
             p = 0.6733333, 
             alternative = c('less'))
  
  print(p$p.value)
  pvals = c(pvals, p$p.value)
  df = c(df, length(which(res$pct_marker[index_group] > 0.03))/length(index_group))
  
}

df = data.frame(group = c('WT',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'), 
                pct_cystOlig2 = c(0.6733333, df),
                pvals = c(NA, pvals))

df$group = factor(df$group, levels = c('WT',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'))

ggplot(data=df, aes(x=group, y=pct_cystOlig2, fill=group)) +
  geom_bar(stat="identity")+
  #theme_minimal() +
  ylab("% cyst with Olig2") + 
  xlab("% genotype Pax6-/-") +
  ylim(0, 0.7) + 
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") + 
  theme_classic() +  
  scale_fill_manual(values=c("darkgreen", rep("steelblue", 5))) +
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'RA_KOKO_day6_DVpatterning_Olig2_vsWT.pdf'), height = 6, width = 8)

write.table(df, file = paste0(outDir, 'RA_KOKO_day6_DVpatterning_Olig2_vsWT_pvals.txt'), 
            sep = '\t', col.names = TRUE, quote = FALSE, row.names = FALSE)

##########################################
# ####### quantifying the Olig2 in WT with RA
##########################################
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
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0),
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title=element_text(size=14,face="bold")) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 18))

ggsave(filename = paste0(outDir, 'TetOnFx_TetOnPx_day4_FoxA2pct_genotypepct_thresholdMean_v3.pdf'), 
       height = 6, width = 9)


##########################################
## d6 stain patterning (markers SHH, NKX22)
## after many test, the automatic thresholds did not work
## the final analysis was done with manual scoring by Elena
##########################################
# res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
#                              "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning/",
#                              #"cyst_size_genotype_Nkx22_Shh_global_manual.csv"),
#                              #"cyst_size_genotype_cyst_Nkx22_Shh_image_manualThresholds.csv"),
#                              # "cyst_size_genotype_cyst_Nkx22_Shh_image_manualThresholds_v2.csv"), 
#                              #"cyst_size_genotype_cyst_Nkx22_Shh_image_manualThresholds_shh250_v3.csv"),
#                              #"cyst_size_genotype_cyst_Nkx22_Shh_image_manualThresholds_shh225_v4.csv"), 
#                              "cyst_size_genotype_cystThresholdsMultiple_Nkx22_Shh_image_manualThresholds_shh300_nkx1000_v5.csv"), 
#                header = TRUE, row.names = c(1))

res = read.csv2(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning_SHH_NKX_manual/", 
                             "d6_TetOn_TetOn_Shh_Nkx22_Filtering_size_pctGenotype.otsu.li.mean.csv"), 
               header = TRUE, sep = "\t")

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})

res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-Nk22-Shh','', res$time)

scores = readxl::read_xlsx(path = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                                         "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning_SHH_NKX_manual/", 
                                         "TetOnTetOn_manual_scoring.xlsx"), 
                           sheet = 1)

scores = data.frame(scores, stringsAsFactors = FALSE)
scores$cyst = paste0(scores$image, '_', scores$cyst_index)

res$cyst = paste0(res$image, '_', res$cyst_index)

mm = match(res$cyst, scores$cyst)

res = data.frame(res, scores[, c(1:5)])

res = res[is.na(res$note), ]

res = res[which(res$treatment != 'noRA'), ]


## first select the size 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(3.7, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5.5), ]

## secondly use three cyst-level thresholds to determine the genotype proportions
res$genotype_pct_otsu = res$nb_fxko_cyst_otsu/(res$nb_fxko_cyst_otsu + res$nb_pxko_cyst_otsu)
res$genotype_pct_mean = res$nb_fxko_cyst_mean/(res$nb_fxko_cyst_mean + res$nb_pxko_cyst_mean)
res$genotype_pct_li = res$nb_fxko_cyst_li/(res$nb_fxko_cyst_li + res$nb_pxko_cyst_li)
#res$genotype_pct_cyst2 = res$nb_total_fxko/(res$nb_total_fxko + res$nb_total_pxko)
#res$genotype_pct_cyst = res$nb_fxko_cyst_mean/(res$nb_fxko_cyst_mean + res$nb_pxko_cyst_mean)
res$treatment = factor(res$treatment, levels = c('dox', 'RA'))

## clean again the manual scores
res$manual_score_nkx22[grep('weak', res$manual_score_nkx22)] = 1

res = res[which(res$manual_score_shh == 1 | res$manual_score_shh == 0), ]
res = res[which(res$manual_score_nkx22 == 1 | res$manual_score_nkx22 == 0), ]


### cleaning a bit the genotype pct 
plot(res$genotype_pct_otsu, res$genotype_pct_mean)
abline(0, 1, lwd = 2.0, col = 'red')
abline(v = c(0.05, 0.95))

plot(res$genotype_pct_otsu, res$genotype_pct_li)
abline(0, 1, lwd = 2.0, col = 'red')


jj = which(res$treatment == 'RA')
length(which(res$manual_score_nkx22[jj] == 1))/length(jj)
length(which(res$manual_score_shh[jj] == 1))/length(jj)


## filtering the genotype pct

#jj1 = which(res$treatment != 'noRA')
#res = res[jj1, ]
res$cutoff_fxko_cyst_otsu = as.numeric(res$cutoff_fxko_cyst_otsu)
res$cutoff_fxko_global_otsu = as.numeric(res$cutoff_fxko_global_otsu)
res$cutoff_pxko_cyst_otsu = as.numeric(res$cutoff_pxko_cyst_otsu)
res$cutoff_pxko_global_otsu = as.numeric(res$cutoff_pxko_global_otsu)

plot((res$cutoff_fxko_cyst_otsu/res$cutoff_fxko_global_otsu),
      (res$cutoff_pxko_cyst_otsu/res$cutoff_pxko_global_otsu))
#jj1 = which(res$treatment == 'noRA')
#points(res$cutoff_fxko_cyst_otsu[jj1]/res$cutoff_fxko_global_otsu[jj1],
#        res$cutoff_pxko_cyst_otsu[jj1]/res$cutoff_pxko_global_otsu[jj1], col = 'red', pch = 2)

abline(v = 2.1)
abline(h = 1.7)

res$outliers = FALSE
res$outliers[which(res$cutoff_fxko_cyst_otsu/res$cutoff_fxko_global_otsu > 2.1 |
                      res$cutoff_pxko_cyst_otsu/res$cutoff_pxko_global_otsu > 1.7)] = TRUE

res = res[!res$outliers, ]

# 
# jj = which(res$treatment == 'RA')
# length(which(res$nb_shh_global_multiotsu[jj]>0))/length(jj)
# length(which(res$pct_shh_genotype[jj]>0))/length(jj)
# 
# plot(res$cutoff_fxko_cyst_otsu, res$cutoff_pxko_cyst_otsu)
# jj1 = which(res$treatment == 'dox')
# points(res$cutoff_fxko_cyst_otsu[jj1], res$cutoff_pxko_cyst_otsu[jj1], col = 'red', pch = 2)
# 
# jj1 = which(res$treatment == 'RA')
# points(res$cutoff_fxko_cyst_otsu[jj1], res$cutoff_pxko_cyst_otsu[jj1], col = 'blue', pch = 2)
# 
# 
# res$outliers[which(res$cutoff_fxko_cyst_otsu > 1.8 |
#                       res$cutoff_pxko_cyst_otsu > 2.5)] = TRUE
#  
# res = res[!res$outliers, ]
#  
# jj = which(res$treatment == 'RA')
# length(which(res$nb_shh_global_multiotsu[jj]>0))/length(jj)
# length(which(res$nb_shh_cyst_otsu[jj]>0))/length(jj)

#saveRDS(res, file = paste0(outDir, 'd6_TetOn_TetOn_Shh_Nkx22_Filtering_size_genotype.rds'))

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


##### Shh 
res$genotype_pct_cyst =res$genotype_pct_otsu

df = c()
cutoffs = seq(0, 1, by = 0.2)

for(n in 1:(length(cutoffs)-1))
{
  # n = 1
  index_group = which(res$treatment == 'dox'& 
                        res$genotype_pct_cyst >=cutoffs[n] & res$genotype_pct_cyst <= cutoffs[n+1])
  
  df = c(df, length(which(res$manual_score_shh[index_group] ==1))/length(index_group))
  
}

index_wt = which(res$treatment == 'RA') 
pct_wt = length(which(res$manual_score_shh[index_wt] ==1))/length(index_wt)

df = data.frame(group = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'), 
                pct_shh = c(pct_wt, df))

df$group = factor(df$group, levels = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'))

df


ggplot(data=df, aes(x=group, y=pct_shh, fill=group)) +
  geom_bar(stat="identity")+
  #theme_minimal() +
  ylab("% cyst with Shh") + 
  xlab("% genotype TetOn-FoxA2") +
  ylim(0, 0.6) + 
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") + 
  theme_classic() +  
  scale_fill_manual(values=c("darkgreen", rep("deepskyblue", 5))) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'TetOn_TetOn_day6_DVpatterning_Shh_vsWT_manual_v2.pdf'), height = 6, width = 8)



######### Nkx22
res$genotype_pct_cyst =res$genotype_pct_otsu

df = c()
cutoffs = seq(0, 1, by = 0.2)

for(n in 1:(length(cutoffs)-1))
{
  # n = 1
  index_group = which(res$treatment == 'dox'& 
                        res$genotype_pct_cyst >=cutoffs[n] & res$genotype_pct_cyst <= cutoffs[n+1])
  
  df = c(df, length(which(res$manual_score_nkx22[index_group] ==1))/length(index_group))
  
}

index_wt = which(res$treatment == 'RA') 
pct_wt = length(which(res$manual_score_nkx22[index_wt] ==1))/length(index_wt)

df = data.frame(group = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'), 
                pct_shh = c(pct_wt, df))

df$group = factor(df$group, levels = c('WT_RA',  '0-20%', '20-40%', '40-60%', '60-80%', '80-100%'))

df


ggplot(data=df, aes(x=group, y=pct_shh, fill=group)) +
  geom_bar(stat="identity")+
  #theme_minimal() +
  ylab("% cyst with Nkx2.2") + 
  xlab("% genotype TetOn-FoxA2") +
  #ylim(0, 1.0) + 
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") + 
  theme_classic() +  
  scale_fill_manual(values=c("darkgreen", rep("deepskyblue", 5))) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'TetOn_TetOn_day6_DVpatterning_Nkx22_vsWT_manual_v2.pdf'), height = 6, width = 8)


########################################################
########################################################
# Section IV : WT_FoxA2 chimeras and WT_Pax6 chimeras
# Fig 4B
########################################################
########################################################
reproduce_figurePlot = FALSE
if(reproduce_figurePlot){
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
  
  
}

res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "FoxAKO_WTchim_d4/",
                             "cyst_size_genotype_FoxA2_Pax6_cystThresholds_test_otsuThresholdGenotype_meanThresholdFoxA2_",
                             "originalGenetype_Normallization_globallocalThreshods_v10.csv"), 
               header = TRUE, row.names = c(1))

res$pct_foxa2_wt = res$nb_foxa2_wt/(res$nb_foxa2_wt + res$nb_pax6_wt)
res$pct_ko = res$nb_ko_all/(res$nb_wt_all + res$nb_ko_all)
res$pct_ko_global = res$nb_ko_all_global/(res$nb_wt_all_global + res$nb_ko_all_global)

plot(res$pct_ko, res$pct_foxa2_wt)
plot(res$pct_ko_global, res$pct_foxa2_wt)

res$pct_foxa2 = res$nb_foxa2_all/(res$nb_foxa2_all + res$nb_pax6_all + res$nb_double_all)
res$pct_foxa2_global = res$nb_foxa2_all_global/(res$nb_foxa2_all_global + res$nb_pax6_all_global)

plot(res$pct_ko, res$pct_foxa2)

plot(res$pct_ko_global, res$pct_foxa2)

plot(res$pct_ko_global, res$pct_foxa2_global)



# kk = which(res$pct_ko_global > 0.02 & res$pct_ko_global < 0.98)
# plot(res$pct_ko_global[kk], res$pct_foxa2[kk])
# 
# 
# plot(res$pct_ko, res$pct_foxa2_global)
# plot(res$pct_ko_global, res$pct_foxa2_global)


# res$pct_foxa2_ko = res$nb_foxa2_ko/(res$nb_foxa2_ko + res$nb_pax6_ko)
# 
# plot(res$pct_ko, res$pct_foxa2_ko)
# 
# res$pct_foxa2_genotype = res$nb_foxa2_genotype/(res$nb_foxa2_genotype + res$nb_pax6_genotype)
# 
# plot(res$pct_ko, res$pct_foxa2_genotype)
# 
# res$pct_foxa2_all = res$nb_foxa2_all/(res$nb_foxa2_all + res$nb_pax6_all)
# 
# plot(res$pct_ko, res$pct_foxa2_all)
# 
# res$pct_foxa2 = (1 - res$pct_ko) * res$pct_foxa2_wt
# plot(res$pct_ko, res$pct_foxa2)
# 
# plot(res$pct_foxa2, res$pct_foxa2_genotype)


## size filtering 
hist(log10(res$cyst_size), breaks = 50)

res = res[which(res$cyst_size > 10^3.5), ]

# plot(res$pct_ko_global, res$pct_ko)
# abline(0, 1, lwd =2.0, col = 'red')
# abline(v = c(0.01, 0.99), col = 'red')
# 
# res$pct_ko[which(res$pct_ko_global > 0.99)] = 1.0
# res$pct_ko[which(res$pct_ko_global < 0.01)] = 0

plot(res$pct_ko, res$pct_foxa2)

plot(res$pct_ko, res$pct_foxa2_global)

plot(res$pct_ko_global, res$pct_foxa2, ylim = range(c(res$pct_foxa2, res$pct_foxa2_global)))
points(res$pct_ko_global, res$pct_foxa2_global, col = 'blue')


plot(res$pct_foxa2_global, res$pct_foxa2, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lwd = 2.0, col = 'blue')
abline(v = 0.01, col = 'red')


#plot(res$pct_ko, res$pct_foxa2)

plot(res$cutoff_wt, res$cutoff_ko)
abline(v = 20)
abline(h = 5)


res = res[which(res$cutoff_wt < 20 & res$cutoff_ko <5), ]

plot(res$pct_ko, res$pct_foxa2_wt)

plot(res$cutoff_foxa2, res$cutoff_pax6)
abline(v = 17)
abline(h = 3.5)

plot(res$pct_foxa2_all, res$pct_foxa2_genotype)
plot(res$pct_foxa2, res$pct_foxa2_genotype)

plot(res$pct_foxa2, res$pct_foxa2_all)
plot(res$pct_ko, res$pct_foxa2_ko);
abline(h = 0.07)

plot(res$pct_foxa2, res$pct_foxa2_all)

plot(res$pct_foxa2, res$pct_foxa2_genotype)
jj0 = which(res$pct_foxa2_ko > 0.1)
points(res$pct_foxa2[jj0], res$pct_foxa2_genotype[jj0], col = 'red')


plot(res$pct_ko, res$pct_foxa2)
jj1 = which(res$cutoff_foxa2 > 6)
points(res$pct_ko[jj1], res$pct_foxa2[jj1], pch = 16, col = 'red')

res = res[which(res$pct_foxa2_ko < 0.1 & res$cutoff_foxa2 < 6), ]

#plot(res$pct_ko, (res$nb_foxa2 + res$nb_pax6))

#res = res[which(res$nb_pax6 > 50), ]

ggplot(res, aes(x=pct_ko, y=pct_foxa2_wt)) +
  geom_point() + 
  geom_smooth(method=loess) +
  ylab("% FoxA2+ in WT ") + 
  xlab("genotype % FoxA2-/-") + 
  theme_classic() +  
  xlim(0, 1) + 
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'FoxA2KO_WTchimeras_day4_v2.pdf'), height = 6, width = 8)


kk = which(res$pct_ko_global < 0.99)
plot(res$pct_ko_global[kk], res$pct_foxa2[kk])

ggplot(res[kk, ], aes(x=pct_ko_global, y=pct_foxa2)) +
  geom_point() + 
  geom_smooth(method = loess) +
  #geom_smooth(method = lm, formula = y ~ splines::ns(x, 2)) +
  ylab("% FoxA2+ in WT ") + 
  xlab("genotype % FoxA2-/-") + 
  theme_classic() +  
  xlim(0, 1) + 
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))

ggsave(filename = paste0(outDir, 'FoxA2_cyst_KOWTchimeras_day4_v3.pdf'), height = 6, width = 8)


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



########################################################
########################################################
# Section : 2xTetOn distribution of states and Nkx22
# 
########################################################
########################################################
##########################################
# only day6
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "d6_2xTetOn/",
                             #"cyst_size_NKX22_FoxA2Pax6_states_maualGlobalThreshods_v1.csv"), 
                             "cyst_size_NKX22_FoxA2Pax6_states_maualGlobalThreshods_DAPIcounts_v2.csv"), 
               header = TRUE, row.names = c(1))

res$genotype = sapply(res$image, function(x){unlist(strsplit(x, '_'))[7]})
res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[8]})
res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[4]})
res$treatment = gsub('-wash-d2', '', res$treatment)

res$condition = paste0(res$genotype, '_', res$treatment)

res$total = res$nb_foxa2_pax6_pp + res$nb_foxa2_pax6_np + res$nb_foxa2_pax6_pn + res$nb_foxa2_pax6_nn
res$pct_nkx = res$nb_global_nkx/res$cyst_size

res$pct_pp = res$nb_foxa2_pax6_pp/res$nb_dapi_cyst_otsu
res$pct_pn = res$nb_foxa2_pax6_pn/res$nb_dapi_cyst_otsu
res$pct_np = res$nb_foxa2_pax6_np/res$nb_dapi_cyst_otsu
res$pct_nn = res$nb_foxa2_pax6_nn/res$nb_dapi_cyst_otsu


## size filtering 
hist(log10(res$cyst_size), breaks = 50)

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5.5), ]

res = data.frame(res, stringsAsFactors = FALSE)

levels_cc = c("wt_noRA", 'wt_RA', 'wt_dox', 
              "n26_noRA", 'n26_RA', 'n26_dox')

res$condition = factor(res$condition, levels = levels_cc)

ggplot(res, aes(x=condition, y=pct_nkx, fill = genotype)) +
  geom_boxplot()+
  ylab("% NTOs with Nkx2.2 ") + 
  xlab("") + 
  #ylim(0, 1.0) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14))


## count the cyst numbers
hist(log10(res$pct_nkx[which(res$condition == 'wt_noRA')]), breaks = 50)

df = data.frame(condition = levels_cc, 
                genotype = res$genotype[match(levels_cc, res$condition)], 
                treatment = res$treatment[match(levels_cc, res$condition)], stringsAsFactors = FALSE)

df$pct_nkx = NA
df$sd = NA
cutoff = 0.01

for(n in 1:nrow(df))
{
  # n = 1
  kk = which(res$condition == df$condition[n])
  df$pct_nkx[n] = length(which(res$pct_nkx[kk] > cutoff)) / length(kk)
  
  
  images = unique(res$image[kk])
  pct_nkx_image = c()
  for(i in 1:length(images))
  {
    # i = 1
    ii_kk = kk[which(res$image[kk] == images[i])]
    pct_nkx_image = c(pct_nkx_image, length(which(res$pct_nkx[ii_kk] > cutoff)) / length(ii_kk))
  }
  
  df$sd[n] = var(pct_nkx_image)
}

library(viridis)
df$condition = factor(df$condition, levels = levels_cc)
df$genotype = factor(df$genotype, levels = c('wt', 'n26'))
#df$sd = 0.2
#df$sd = sqrt(df$sd)

ggplot(df, aes(x=condition, y=pct_nkx, fill = genotype)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=pct_nkx-sd, ymax=pct_nkx+sd), width=.2,
                position=position_dodge(.9)) +
  ylab("% NTOs with Nkx2.2 ") + 
  xlab("") + 
  ylim(0, 1.0) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 0.6),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14)) +
  scale_fill_manual(values=c("darkgreen", "#31688EFF")) 
  #scale_fill_manual(values=c("darkgreen", rep("deepskyblue", 5))) +

ggsave(filename = paste0(outDir, '2xTetOn_day6_percentages_NTOs_postiveNKX22_v2.pdf'), 
       width = 6, height = 4)



kk_sels = match(c('condition', 'image', 'pct_pp', 'pct_pn', 'pct_np', 'pct_nn'), colnames(res))
xx = res[, kk_sels]
xx = xx %>%
  tidyr::pivot_longer(starts_with("pct"), names_to = "states", values_to = "pct")

xx = data.frame(xx)
xx$cc = paste0(xx$condition, "_", xx$states)

df = xx[match(unique(xx$cc), xx$cc), ]
df$sd_pct = 0
df$sd_pct2 = 0

for(n in 1:nrow(df))
{
  # n = 1
  df$pct[n] = mean(xx$pct[which(xx$cc == df$cc[n])])
  df$sd_pct2[n] = var(xx$pct[which(xx$cc == df$cc[n])])
  
  jj = which(xx$cc == df$cc[n])
  images = unique(xx$image[jj])
  pct_image = c()
  for(i in 1:length(images))
  {
    ii_jj = jj[which(xx$image[jj] == images[i])]
    pct_image = c(pct_image, mean(xx$pct[ii_jj]))
  }
  
  df$sd_pct[n] = var(pct_image)
}

df$states = factor(df$states, levels = c('pct_pn', 'pct_pp', 'pct_np', 'pct_nn'))

df$sd_pct = sqrt(df$sd_pct2)

error_bars = df %>%
  arrange(condition, desc(states)) %>%
  # for each cyl group, calculate new value by cumulative sum
  group_by(condition) %>%
  mutate(mean_hp_new = cumsum(pct)) %>%
  ungroup()

ggplot(df, aes(x = condition, y = pct)) +
  geom_bar(stat = 'identity', aes(fill = states)) +
  geom_errorbar(data = error_bars,
                aes(x = condition, ymax = mean_hp_new + sd_pct, ymin = mean_hp_new - sd_pct), 
                width = 0.2, position=position_dodge(.9)) +
# ggplot(data=df, aes(x=condition, y=pct, fill=states)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=pct-sd, ymax=pct+sd), width=.2,
#                 position=position_dodge(.9)) +
  ylab("% of states ") + 
  xlab("") + 
  ylim(0, 1.2) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14)) +
  scale_fill_manual(values=c("darkgreen", "darkorange", "red", 'gray')) 

ggsave(filename = paste0(outDir, '2xTetOn_day6_percentages_cystStates_v3.pdf'), 
       width = 8, height = 5)

# ## distribution of pn states  
# cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
# 
# jj = which(res$condition == 'wt_RA' | res$condition == "n26_RA" | res$condition == "n26_dox")
# ggplot(res[jj, ], aes(x=pct_pn, fill = condition)) +
#   geom_density(alpha = 0.7) +
#   scale_fill_manual(values = c("darkgreen", "green2", "#72D8FF")) +
#   ylab("Density of % FoxA2+ Pax6- ") + 
#   xlab("") + 
#   #ylim(0, 1.0) + 
#   theme_classic() +  
#   theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
#         axis.text.y = element_text(angle = 0, size = 14)) +
#   theme(legend.key = element_blank()) + 
#   theme(plot.margin=unit(c(1,3,1,1),"cm"))+
#   #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
#   theme(legend.title = element_blank(), 
#         legend.text = element_text(size = 14))
# 
# ggsave(filename = paste0(outDir, '2xTetOn_day6_densityPlot_percentages_statePN.pdf'), 
#        width = 8, height = 5)
# 
# 
# jj = which(res$condition == 'wt_RA' | res$condition == "n26_RA" | res$condition == "n26_dox")
# ggplot(res[jj, ], aes(x=pct_np, fill = condition)) +
#   geom_density(alpha = 0.7) +
#   scale_fill_manual(values = c("darkgreen", "green2", "#72D8FF")) +
#   ylab("Density of % FoxA2+ Pax6- ") + 
#   xlab("") + 
#   #ylim(0, 1.0) + 
#   theme_classic() +  
#   theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
#         axis.text.y = element_text(angle = 0, size = 14)) +
#   theme(legend.key = element_blank()) + 
#   theme(plot.margin=unit(c(1,3,1,1),"cm"))+
#   #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
#   theme(legend.title = element_blank(), 
#         legend.text = element_text(size = 14))
# 
# ggsave(filename = paste0(outDir, '2xTetOn_day6_densityPlot_percentages_stateNP.pdf'), 
#        width = 8, height = 5)
# 


########################################################
########################################################
# Section X: analyze the embryo data
# 
########################################################
########################################################
inputDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/', 
                  'embryo_features/test_v2/')

list_files = list.files(path = inputDir, pattern = '*.csv', full.names = TRUE)

res = c()
for(n in 1:length(list_files))
{
  # n = 1
  cat(n, basename(list_files[n]), '\n')
  xx = read.csv(file = list_files[n])
  xx = xx[, c(2:9, grep('^area_foxa2|intensity_mean', colnames(xx)))]
  cc = gsub('.csv', '', basename(list_files[n]))
  cc = gsub('featuresCollection_250712_30xsil-041umZ_|featuresCollection_250713_30xsil-041umZ_', '', cc)
  xx = data.frame(condition = rep(cc, nrow(xx)), xx)
  res = rbind(res, xx)
  rm(xx)
  
}

colnames(res)[(ncol(res)-4):ncol(res)] = c('area', 'foxa2', 'pax6', 'sox2', 'dapi')

res$condition = gsub('E85-','', res$condition)
res$condition = gsub('crop-','', res$condition)
res$embryo = sapply(res$condition, function(x){unlist(strsplit(as.character(x), '_'))[1]})
res$position = sapply(res$condition, function(x){unlist(strsplit(as.character(x), '_'))[2]})

plot(res$volume, res$area, log = 'xy')
abline(0, 1, lwd = 2.0, col = 'red')

res$area = log10(res$volume)

hist(res$area, breaks = 100)
abline(v = c(4, 4.9))

plot(res$sphericity_legland, res$solidity)
plot(res$sphericity_wadell, res$solidity)

plot(res$sphericity_legland, res$sphericity_wadell)
abline(h = 0.8)

plot(res$area, res$sphericity_wadell, cex = 0.1)
abline(v = c(4., 4.9), col = 'red', lwd = 1.5)
abline(h = c(0.6), col = 'red', lwd = 1.5)

res = res[which(res$area > 4. & res$area < 4.9 & res$sphericity_wadell > 0.6), ]


res$foxa2 = log10(res$foxa2)
res$pax6 = log10(res$pax6)

hist(res$foxa2, breaks = 100)
abline(v = c(3.2))

hist(res$pax6, breaks = 100)
abline(v = c(3.6, 3.2))

saveRDS(res, file = paste0(outDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
#saveRDS(res, file = paste0(outDir, 'embryo_features_filtered.size_v1.rds'))

embs = unique(res$embryo)
cc = unique(res$position)

res = data.frame(res)

cc = c("CLE2",  "CLE1", "somite7", "somite6", "somite5", "somite4", "somite3", "somite2", "somite1")



library(ggpubr)

for(e in embs)
{
  # e = embs[3]
  cc_emb = cc[which(!is.na(match(cc, unique(res$position[which(res$embryo == e)]))))]
  
  for(n in 1:length(cc_emb))
  {
    # n =5
    c = cc_emb[n]
    kk = which(res$position == cc_emb[n] & res$embryo == e)
    
    eval(parse(text = paste0("p", n, " = ggplot(res[kk, ], aes(x=foxa2, y=pax6)) +
    geom_point(size = 1) +
    geom_density_2d() + 
    xlim(2.5, 4.0) +
    ylim(2.5, 4.0) +
    geom_hline(yintercept = c(3.2)) +
    geom_vline(xintercept = 3.0) + 
    theme_classic() + 
    ggtitle(c)            
                           ")))
    
  }
  
  if(length(cc_emb) == 8){
    ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, 
              #labels = c("A", "B", "C"),
              ncol = 4, nrow = 2)
    ggsave(filename = paste0(outDir, 'firstTeset_embryo_', e, '.pdf'),  
           width = 16, height = 6)
    
  }
  
  
  if(length(cc_emb) == 9){
    ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,
              #labels = c("A", "B", "C"),
              ncol = 3, nrow = 3)
    
    ggsave(filename = paste0(outDir, 'firstTeset_embryo_', e, '.pdf'),  
           width = 16, height = 12)
    
  }
  
}


##########################################
# plot test for embyo H and L 
##########################################
res  = readRDS(file = paste0(outDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
res = data.frame(res)

res = res[which(res$embryo != 'embB'), ]

embs = unique(res$embryo)
cc = unique(res$position)

cc = c("CLE2",  "CLE1", "somite7", "somite6", "somite5", "somite4", "somite3", "somite2", "somite1")

for(n in 1:length(cc))
{
  # n = 1
  c = cc[n]
  kk = which(res$position == c)
  
  plot = ggplot(res[kk, ], aes(x=foxa2, y=pax6, color = embryo)) +
    geom_point(size = 1) +
    #geom_density_2d() + 
    xlim(2.5, 4.0) +
    ylim(2.5, 4.0) +
    geom_hline(yintercept = c(3.2)) +
    geom_vline(xintercept = 3.0) + 
    theme_classic() + 
    ggtitle(c) +
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    theme(legend.position = c(0.9,1.0), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14)) +
    xlab("FoxA2") + 
    ylab("Pax6") 
  
  ggsave(filename = paste0(outDir, 'firstTest_embryoHI_', c, '.pdf'),  
         width = 8, height = 6)
  
}


## plot the Embryo H and I with position somite4-1
res  = readRDS(file = paste0(outDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
res = data.frame(res)

res = res[which(res$embryo != 'embB'), ]

cc = c("somite4", "somite3", "somite2", "somite1")
res = res[which(!is.na(match(res$position, cc))), ]

ggplot(res, aes(x=foxa2, y=pax6, color = embryo, shape=position)) +
  geom_point(size = 1.5) +
  #geom_density_2d() + 
  xlim(2.7, 3.8) +
  ylim(2.8, 3.7) +
  geom_hline(yintercept = c(3.2)) +
  geom_vline(xintercept = 3.0) + 
  theme_classic() + 
  #ggtitle(c) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position="top") + 
  #theme(legend.text = element_text(colour="blue", size=10, face="bold")) +
  guides(size = 2) +
  #theme(legend.position = c(0.9,1.0), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 16)) +
  xlab("FoxA2") + 
  ylab("Pax6") +
  scale_color_manual(values=c('#E69F00', '#56B4E9'))

ggsave(filename = paste0(outDir, 'scatterPlot_embryoHI.pdf'),  
       width = 8, height = 6)


##########################################
# plot FoxA2 and Pax6 intensity for only somite 1-4
##########################################
res  = readRDS(file = paste0(outDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
res = data.frame(res)

embs = unique(res$embryo)
cc = unique(res$position)

cc = c("somite4", "somite3", "somite2", "somite1")

res = res[which(!is.na(match(res$position, cc))), ]


for(n in 1:length(embs))
{
  # n = 1
  e = embs[n]
  kk = which(res$embryo == e)
  
  plot = ggplot(res[kk, ], aes(x=foxa2, y=pax6, color = position)) +
    geom_point(size = 1) +
    #geom_density_2d() + 
    xlim(2.7, 4.0) +
    ylim(2.7, 4.0) +
    #geom_hline(yintercept = c(3.2)) +
    #geom_vline(xintercept = 3.0) + 
    theme_classic() + 
    ggtitle(e) +
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(0.9,1.0), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 16)) +
    xlab("FoxA2") + 
    ylab("Pax6") +
    scale_color_brewer(palette = "Set1")
  plot
  #ggsave(filename = paste0(outDir, 'firstTest_embryoHI_', c, '.pdf'),  
  #       width = 8, height = 6)
  ggsave(filename = paste0(outDir, 'ScatterPlot_somite1_4_embryo_', e, '.pdf'),  
         width = 8, height = 6)
  
}


##########################################
# multi-gaussian clustering of FoxA2, Pax6 and Sox2 
##########################################
resDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/results/",
                "figures_tables_R13547_10x_mNT_20240522/")
res  = readRDS(file = paste0(resDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
res = data.frame(res)

outDir = paste0(resDir, '/embryo_mclusters/')
if(!dir.exists(outDir)) dir.create(outDir)


embs = unique(res$embryo)
cc = unique(res$position)

cc = c("somite4", "somite3", "somite2", "somite1")

res = res[which(!is.na(match(res$position, cc))), ]

res$pax6 = 10^res$pax6
res$foxa2 = 10^res$foxa2

library(mclust)
source(paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/scripts/",
              "functions_plotMclust.R"))

for(n in 1:length(embs))
{
  # n = 1
  e = embs[n]
  kk = which(res$embryo == e)
  
  #hist(res$foxa2[kk], breaks = 100)
  #hist(res$pax6[kk], breaks = 100)
  #hist(res$sox2[kk], breaks = 100)
  
  plot = ggplot(res[kk, ], aes(x=pax6, y=foxa2, color = position)) +
    geom_point(size = 1) +
    #geom_density_2d() + 
    #xlim(2.7, 4.0) +
    #ylim(2.7, 4.0) +
    #geom_hline(yintercept = c(3.2)) +
    #geom_vline(xintercept = 3.0) + 
    theme_classic() + 
    ggtitle(e) +
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(0.9,1.0), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 16)) +
    xlab("Pax6") + 
    ylab("Foxa2") +
    scale_color_brewer(palette = "Set1")
  plot
  
  ggsave(filename = paste0(outDir, 'ScatterPlot_somite1_4_embryo_', e, '.pdf'),  
         width = 10, height = 6)
  
  ##########################################
  # pooling all treatment and time
  ##########################################
  mat = res[kk, c(11:13)]
  #mat[, 1] = 10^mat[,1]
  #mat[, 2] = 10^mat[,2]
  metadata = res[kk, c(1, 15:16) ]
  
  #print(dim(mat))
  #print(dim(metadata))
  
  for(nb_clusters in c(3:5))
  {
    # nb_clusters = 4
    cat(e, ': nb of cluster -- ', nb_clusters, '\n')
    
    set.seed(1000)
    mb = Mclust(mat, G = nb_clusters, control = emControl(itmax=500, tol = 1.e-6))
    
    cat('loglike --', mb$loglik, "\n")
    
    # probality for an observation to be in a given cluster
    #head(mb$z)
    
    # get probabilities, means, variances
    #summary(mb, parameters = TRUE)
    
    #saveRDS(mb, file = paste0(outDir, '/res_mclust_nbClusters.', nb_clusters, '.rds'))
    
    clusters = mb$classification
    clusters = clusters[match(rownames(mat), names(clusters))]  
    
    # keep = table(metadata$condition, mb$classification)
    # 
    # #Compare amount of the data within each cluster
    # write.csv2(keep, file = paste0(outDir, '/cellNumbers_perCluster_perCondition_nbClusters_', 
    #                                nb_clusters, '.csv'))
    # 
    # for(n in 1:nrow(keep)){
    #   keep[n, ] = keep[n, ]/sum(keep[n,])
    # }
    # 
    # write.csv2(keep, file = paste0(outDir, '/cellProportions_perCluster_perCondition_nbCluste_', 
    #                                nb_clusters, '.csv'))
    
    keep = data.frame(mat, metadata[match(rownames(mat), rownames(metadata)), ], stringsAsFactors = FALSE)
    keep = data.frame(keep, clusters, stringsAsFactors = FALSE)
    
    table(keep$condition, keep$clusters)
    write.csv2(keep, file = paste0(outDir, 'data_metadata_clusterIDs_', e, '_nbClusters_',
                                  nb_clusters, '.csv'))
    
    #metadata$cluster = mb$classification
    cc = unique(clusters)
    cc = cc[order(cc)]
    
    pdf(paste0(outDir, "markerIntensity_", e, "_nbClusters_", nb_clusters, ".pdf"), 
        height = 3*nb_clusters, width =16, useDingbats = FALSE)
    
    attach(mtcars)
    par(mfrow=c(length(cc), ncol(mat))) 
    for(n in 1:length(cc))
    {
      c = cc[n];
      for(m in 1:ncol(mat))
      {
        # c = 1; m = 1;
        hist(mat[which(clusters == c), m], breaks = 50, xlim = range(mat),
             xlab = '', ylab = paste0('cluster_', c), main = colnames(mat)[m],
             col = n);
        
      }
    }
    
    dev.off()
    
    
    pdf(paste0(outDir, "clusterProjection_", e, "_nbClusters_",  nb_clusters, ".pdf"), 
        height = 12, width =16, useDingbats = FALSE)
    
    #After the data is fit into the model, we plot the model based on clustering results.
    # plot(mb, "density")
    
    
    # S3 method for Mclust
    # plot(mb, what = c("classification", "density"), 
    #      dimens = NULL, xlab = NULL, ylab = NULL,
    #      addEllipses = TRUE, main = FALSE)
    
    plot.Mclust_cutomized(mb, what=c("classification"), cex = 0.5, 
                          addEllipses = TRUE, cex_clusterlabels = 2.0)
    
    
    #plot.surface_customized(mb)
    dev.off()
    
    
  }
  
  
}




##########################################
# percentages of FoxA2+, Pax6+ and ++ cells based on the image thresholding
##########################################
inputDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/', 
                  'embryo_features/test_foxa2_pax6')

list_files = list.files(path = inputDir, pattern = '*.csv', full.names = TRUE)

res = c()
for(n in 1:length(list_files))
{
  # n = 1
  cat(n, basename(list_files[n]), '\n')
  xx = read.csv(file = list_files[n])
  xx = xx[, -1]
  cc = gsub('.csv', '', basename(list_files[n]))
  cc = gsub('250712_30xsil-041umZ_E85-|embryo_voxel_counts_global_otsu_li_mean|250713_30xsil-041umZ_E85-', '', cc)
  xx = data.frame(condition = rep(cc, nrow(xx)), xx)
  res = rbind(res, xx)
  rm(xx)
  
}

colnames(res)[3:4] = c('cell_index', 'cell_size')

#res$condition = gsub('E85-','', res$condition)
res$condition = gsub('crop-','', res$condition)
res$embryo = sapply(res$condition, function(x){unlist(strsplit(as.character(x), '_'))[1]})
res$position = sapply(res$condition, function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$id = paste0(res$condition, '_', res$cell_index)

xx = readRDS(file = paste0(outDir, 'embryo_features_filtered.size.sphericity_v2.rds'))
xx$id = paste0(xx$condition, '_', xx$label)

mm = match(xx$id, res$id)

res = res[mm, ]

#saveRDS(res, file = paste0(outDir, 'embryo_FoxA2_Pax6_detection.rds'))

## define FoxA2+, Pax6+ or FoxA2+Pax6+
res = readRDS(file = paste0(outDir, 'embryo_FoxA2_Pax6_detection.rds'))
res$pct_foxa2 = res$nb_foxa2_otsu_global/res$cell_size
res$pct_pax6 = res$nb_pax6_otsu_global/res$cell_size

embs = unique(res$embryo)
cc = unique(res$position)
cc = c("CLE2",  "CLE1", "somite7", "somite6", "somite5", "somite4", "somite3", "somite2", "somite1")

for(n in 1:length(embs))
{
  # n = 1
  e = embs[n]
  kk = which(res$embryo == e)
  xx = res[kk, ]
  
  c = unique(xx$position)
  
  pcts = c()
  cutoff = 0.02
  for(i in 1:length(c))
  {
    jj = which(xx$position == c[i])
    
    index_foxa2 = which(xx$pct_foxa2[jj] > cutoff & xx$pct_pax6[jj] < cutoff)
    index_pax6 = which(xx$pct_foxa2[jj] < cutoff & xx$pct_pax6[jj] > cutoff)
    index_double = which(xx$pct_foxa2[jj] > cutoff & xx$pct_pax6[jj] > cutoff)
    
    pcts = rbind(pcts, c(e, c[i], 'FoxA2+', length(index_foxa2)/length(jj)))
    pcts = rbind(pcts, c(e, c[i], 'Pax6+', length(index_pax6)/length(jj)))
    pcts = rbind(pcts, c(e, c[i], 'doublePos', length(index_double)/length(jj)))
    pcts = rbind(pcts, c(e, c[i], 'doubleNeg', 
                         (length(jj) - length(index_double) - length(index_foxa2) - length(index_pax6))/length(jj)
                         )
                 )
  }
  
  pcts = data.frame(pcts)
  colnames(pcts) = c('embryo', 'position', 'group', 'pct')
  pcts$position = factor(pcts$position, levels = cc)
  pcts$position = droplevels(pcts$position)
  pcts$pct = as.numeric(pcts$pct)
  pcts$group = factor(pcts$group, levels = c('FoxA2+', 'Pax6+', 'doublePos', 'doubleNeg'))
  
  
  ggplot(data=pcts, aes(x=position, y=pct, fill=group)) +
    #geom_bar(stat="identity", color="black", position=position_dodge())+
    geom_bar(stat="identity") + 
    theme_minimal() + 
    scale_fill_manual(values=c('darkgreen',"red", '#E69F00', 'gray')) + # Use custom colors
    ggtitle(e) +
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 14)) +
    theme(legend.key = element_blank()) + 
    #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
    #theme(legend.position = c(0.9,1.0), legend.direction = "vertical") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 14)) +
    xlab("") + 
    ylab("%") 
  
  ggsave(filename = paste0(outDir, 'firstTest_otsuThresholding_embryo_', e, 'doublePos_pct.pdf'),  
         width = 12, height = 6)
  
}
