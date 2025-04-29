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

########################################################
########################################################
# Section I : distribution of FoxA2+ in WT (figure 1K and 1L) 
# 
########################################################
########################################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                             "images_data/results/test_WT_FoxA2_pct/wt_FoxA2_Pax6_voxel_counts.csv"), 
               header = TRUE, row.names = c(1))

res$pct_foxa2 = res$nb_foxa2_postive/res$total_postives

plot(res$cyst_size, res$pct_foxa2)
abline(v = c(5000, 10000), lwd = 2.0, col = 'red')

hist(res$cutoff_foxa2, breaks = 40)
abline(v = c(0.6, 2.0), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_foxa2 > 0.6 & res$cutoff_foxa2 < 2.0), ]

plot(res$cyst_size, res$pct_foxa2)
abline(v = c(5000, 10000), lwd = 2.0, col = 'red')

hist(res$cutoff_pax6, breaks = 40)
abline(v = c(0.6, 2.0), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_pax6 > 0.6 & res$cutoff_foxa2 < 2.0), ]

plot(res$cyst_size, res$pct_foxa2)
abline(v = c(5000, 10000), lwd = 2.0, col = 'red')


plot(res$total, res$pct_foxa2)
plot(res$cyst_size, res$total) ## here total is the number of FoxA2+, Pax6+ and FoxA2+&Pax6+ 
abline(a = 0, b = 0.5, lwd = 2.0, col = 'red')

res = res[which(res$cyst_size > 10000), ]

hist(res$pct_foxa2, breaks = 20)
plot(density(res$pct_foxa2, adjust = 1.5))

xx = log10(res$cyst_size[which(res$cyst_size>10^4)])
hist(xx,  probability = TRUE)
lines(density(xx), lwd = 2, col = "chocolate3")

rd = rnorm(n = 500, mean = mean(xx), sd = sd(xx))
lines(density(rd), lwd = 2, col = "red")

set.seed(2025)
a = mean(xx)/log10(200)
rd2 = rnorm(n = 10000, mean = mean(xx)/a, sd = sd(xx)/a)
rd2 = 10^rd2
hist(rd2, breaks = 20)


rd2 = floor(rd2)

pct = c()
pct2 = c()
set.seed(2025)
for(size in rd2)
{
  prob = runif(n = 1, min = 0, max = 1);
  pct = c(pct, rbinom(n = 1, size = size, prob = prob)/size)
  pct2 = c(pct2, rbinom(n = 1, size = size, prob = 0.3)/size)
  
}



plot(density(res$pct_foxa2, adjust = 1.5), col = 'darkgreen', lwd = 3.0,  xlim = c(0, 1), ylim=c(0, 4))
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

#res = res[which(res$time == 'd4' & res$treatment == 'RA'), ]
res = res[which(res$treatment == 'RA'), ]

res$genotype_pct_pxko = res$nb_pxko/res$genotype_total
res$pct_foxa2 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/res$genotype_total*2

## size filtering 
hist(log10(res$cyst_size), breaks = 50)


res = res[which(res$cyst_size > 10^4), ]

plot(res$cutoff_foxa2, res$cutoff_pax6)
abline(v = 0.25, col = 'red')
abline(h = 0.25, col = 'red')

hist(res$cutoff_foxa2, breaks = 40)
abline(v = c(0.25, 1.7), lwd = 2.0, col = 'red')
res = res[which(res$cutoff_foxa2 > 0.2 & res$cutoff_foxa2 < 1.7), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_pax6, breaks = 40)
abline(v = c(0.25, 1.2), lwd = 2.0, col = 'red')
res = res[which(res$cutoff_foxa2 < 1.2 & res$cutoff_pax6 > 0.25), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')


hist(res$cutoff_wt, breaks = 40)
abline(v = c(0.1, 1.5), lwd = 2.0, col = 'red')
res = res[which(res$cutoff_wt > 0.1 & res$cutoff_foxa2 < 1.5), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_ko, breaks = 20)
abline(v = c(0.2, 0.9), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_ko < 0.9), ]


plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

# Add regression lines
ggplot(res, aes(x=genotype_pct_pxko, y=pct_foxa2, color=time, shape=time)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=time))+
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

ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day3_4_6_FoxA2pct.pdf'), height = 6, width = 10)

res = res[which(res$time == 'd4'), ]
ggplot(res, aes(x=genotype_pct_pxko, y=pct_foxa2)) +
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

ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day4_FoxA2pct.pdf'), height = 6, width = 10)


##########################################
# Pax6KO and FoxA2KO quantifying the genotypes and marker genes  
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "PKO_FKO_d6_stain_DVpatterning/",
                             "cyst_size_genotype_FoxA2_markers_cystThresholds_noMarkerNormalization.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})
res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-2umZ','', res$time)
res$marker = sapply(res$image, function(x){unlist(strsplit(x, '_'))[6]})
res$marker = gsub('FA2-', '', res$marker)

res$genotype_pct_pxko = res$nb_pxko/res$genotype_total
res$pct_foxa2 = (res$nb_foxa2_fxko + res$nb_foxa2_pxko)/res$genotype_total
res$pct_marker = (res$nb_marker_fxko + res$nb_marker_pxko)/res$genotype_total

res$pct_foxa2_cyst = res$nb_foxa2_cyst/res$cyst_size
res$pct_marker_cyst = res$nb_marker_cyst/res$cyst_size

res$pct_foxa2_global = res$nb_foxa2_cyst_cutoffImage/res$cyst_size
res$pct_marker_global = res$nb_marker_cyst_cutoffImage/res$cyst_size

plot(res$pct_foxa2, res$pct_foxa2_cyst)

plot(res$pct_foxa2_global, res$pct_foxa2_cyst)
abline(0, 1, lwd =2.0, col = 'red')

plot(res$pct_foxa2, res$pct_foxa2_global)
abline(0, 1, lwd =2.0, col = 'red')

saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_all_nomarkerNormalization.rds'))

kk = which(res$treatment == 'RA')
#res = res[which(res$time == 'd4' & res$treatment == 'RA'), ]
res = res[kk, ]


## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(4, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^4 & res$cyst_size < 10^5.5), ]

saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$cutoff_foxa2, res$cutoff_marker)
abline(v = 0.25, col = 'red')
abline(v = 1.5, col = 'red')

hist(res$cutoff_foxa2, breaks = 40)
abline(v = c(1.5, 8), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_foxa2 > 1.5 & res$cutoff_foxa2 < 8), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_marker, breaks = 40)
abline(v = c(0.25, 1.2), lwd = 2.0, col = 'red')
res = res[which(res$cutoff_marker < 2.0), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')


hist(res$cutoff_fxko, breaks = 40)
abline(v = c(0.2, 2), lwd = 2.0, col = 'red')

res = res[which(res$cutoff_fxko > 0.2 & res$cutoff_fxko < 2.0), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

hist(res$cutoff_pxko, breaks = 40)
abline(v = c(0.2, 4.2), lwd = 2.0, col = 'red')

plot(res$cutoff_fxko, res$cutoff_pxko)
res = res[which(res$cutoff_pxko < 4.2), ]

res$pct_foxa2_false = res$nb_foxa2_fxko/res$nb_fxko

hist(res$pct_foxa2_false, breaks = 50)
abline(v = c(0.1, 0.05), lwd = 2.0, col = 'red')
#hist(res$cutoff_fxko, breaks = 20)

res = res[which(res$pct_foxa2_false < 0.05), ]

plot(res$genotype_pct_pxko, res$pct_foxa2)
abline(0, 1, lwd = 2.0, col = 'red')

## plot the pct FoxA2 in function of genotype pct
ggplot(res, aes(x=genotype_pct_pxko, y=pct_foxa2)) +
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
                             "cyst_size_genotype_FoxA2_Pax6_cystThresholds_globalThreshold_quantile.csv"), 
               header = TRUE, row.names = c(1))

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
res$pct_foxa2_cyst = res$nb_foxa2_cyst/(res$nb_foxa2_cyst + res$nb_pax6_cyst)

#res$genotype_pct_pxko = res$nb_fxko/res$genotype_total
#res$pct_foxa2 = (res$fxko_nb_foxa2+res$pxko_nb_foxa2)/res$genotype_total*2


## select only the d4
res = res[which(res$time == 'd4'), ]

#res = res[which(res$time == 'd4' & res$treatment == 'RA'), ]
#res = res[which(res$treatment == 'RA'), ]

plot(res$genotype_pct_cyst, res$genotype_pct_global);abline(0, 1, col = 'red', lwd = 2.0)


plot(res$pct_foxa2_cyst, res$pct_foxa2_global);abline(0, 1, col = 'red', lwd = 2.0)
#plot(res$pct_foxa2_cyst, res$pct_foxa2_cyst2);abline(0, 1, col = 'red', lwd = 2.0)

plot(res$genotype_pct_global, res$pct_foxa2_global)

########################
## filtering steps 
## size filtering 
hist(log10(res$cyst_size), breaks = 50)

hist(log10(res$cyst_size[which(res$treatment == 'dox')]), breaks = 100)

res = res[which(res$cyst_size > 10^3.5 & res$cyst_size < 10^5), ]

table(res$treatment)

plot(res$cutoff_fxko, res$cutoff_pxko, cex = 0.7)
abline(v = c(150, 1400), col = 'red')
abline(h = 1000, col = 'red')
kk = which(res$treatment == 'dox')
points(res$cutoff_fxko[kk], res$cutoff_pxko[kk], col = 'blue', cex = 0.7)
kk = which(res$treatment == 'noRA')
points(res$cutoff_fxko[kk], res$cutoff_pxko[kk], col = 'red', cex = 0.7, pch =4)

hist(res$cutoff_fxko, breaks = 100);abline(v = c(150, 1400))
hist(res$cutoff_pxko, breaks = 100);abline(v = c(1200))

jj1 = which(res$cutoff_fxko > 150  & res$cutoff_fxko < 1400 & res$cutoff_pxko > 1000)
res = res[jj1, ]
table(res$treatment)

plot(res$cutoff_foxa2, res$cutoff_pax6, xlim = c(0, 7000), ylim = c(0, 1500), type = 'n')
abline(v = 1500, col = 'red')
abline(h = 400, col = 'red')
kk = which(res$treatment == 'noRA')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'red', pch = 2)

kk = which(res$treatment == 'dox')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'blue', pch = 16)

kk = which(res$treatment == 'RA')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'orange', pch = 4)

kk = which(res$treatment == 'RAdox')
points(res$cutoff_foxa2[kk], res$cutoff_pax6[kk], col = 'green', pch = 1)


hist(res$cutoff_foxa2[which(res$treatment == 'RA')], breaks = 40)
abline(v = c(1000, 4000))

hist(res$cutoff_foxa2[which(res$treatment == 'dox')], breaks = 40)
abline(v = c(1000, 6000))

hist(res$cutoff_pax6[which(res$treatment == 'dox')], breaks = 40)
abline(v = c(1000, 6000))

hist(res$cutoff_pax6[which(res$treatment == 'RAdox')], breaks = 40)
abline(v = c(450, 1500))


j1 = which(res$cutoff_foxa2 > 1000 & res$cutoff_foxa2 <4000 & res$treatment == 'RA')
j2 = which(res$cutoff_pax6 > 450 & res$cutoff_pax6 < 1500 & res$treatment == 'RAdox')
j3 = which(res$cutoff_foxa2 > 1000 & res$cutoff_foxa2 <6000 & res$cutoff_pax6 > 250 & res$treatment == 'dox')

jj2 = c(j1, j2, j3)

res =res[jj2, ]

aa = res[which(res$treatment != 'noRA'), ]

ggplot(aa, aes(x=genotype_pct_cyst, y=pct_foxa2_cyst, color=treatment)) +
  geom_point() + 
  geom_smooth(method=loess, aes(fill=treatment))+
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

ggsave(filename = paste0(outDir, 'TetOnFx_TetOnPx_day4_FoxA2pct_genotypepct_v2.pdf'), height = 6, width = 10)



ggplot(aa, aes(x=genotype_pct_global, y=pct_foxa2_global, color=treatment)) +
  geom_point() + 
  #geom_smooth(method=lm, aes(fill=treatment))+
  geom_smooth(method=loess, aes(fill=treatment))+
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

ggplot(aa, aes(x=genotype_pct_global2, y=pct_foxa2_global2, color=treatment)) +
  geom_point() + 
  #geom_smooth(method=lm, aes(fill=treatment))+
  geom_smooth(method=loess, aes(fill=treatment))+
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



ggsave(filename = paste0(outDir, 'TetOnFx_TetOnPx_day4_FoxA2pct_genotypepct.pdf'), height = 6, width = 10)





ggsave(filename = paste0(outDir, 'TetOnFx_TetOnPx_day4_FoxA2pct_genotypepct_v2.pdf'), height = 6, width = 10)



p3 = ggplot(aa, aes(x=genotype_pct_pxko, y=pct_foxa2, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=treatment))+
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


p1 + p2 + p3

##########################################
# d6 stain patterning (markers SHH, NKX22) 
##########################################
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/",
                             "TetOnF_TetOnP_chim_d6_genotype_stainDVpatterning/",
                             "cyst_size_genotype_Nkx22_Shh_cystThresholds_globalThreshold_quantile95.csv"), 
               header = TRUE, row.names = c(1))

res$treatment = sapply(res$image, function(x){unlist(strsplit(x, '_'))[5]})

res$time = sapply(res$image, function(x){unlist(strsplit(x, '_'))[3]})
res$time = gsub('-Nk22-Shh','', res$time)

res$genotype_pct_global = res$nb_fxko_global/(res$nb_fxko_global + res$nb_pxko_global)
res$genotype_pct_global2 = res$nb_fxko_global/res$cyst_size

res$genotype_pct_cyst = res$nb_fxko_cyst/(res$nb_fxko_cyst + res$nb_pxko_cyst)
res$genotype_pct_cyst2 = res$nb_fxko_cyst/res$cyst_size


res$genotype_pct_pxko = res$nb_fxko/res$genotype_total

res$pct_nkx22_global = res$nb_foxa2_global/res$cyst_size
res$pct_nkx22_cyst = res$nb_foxa2_cyst/res$cyst_size

res$pct_shh_global = res$nb_pax6_global/res$cyst_size
res$pct_shh_cyst = res$nb_pax6_cyst/res$cyst_size

## size filtering 
hist(log10(res$cyst_size), breaks = 50)
abline(v = c(3.7, 5.5), col = 'red')

res = res[which(res$cyst_size > 10^4. & res$cyst_size < 10^5.5), ]

## filtering 
plot((res$nb_fxko_global + res$nb_pxko_global), res$cyst_size, log = 'xy')
abline(0, 1, lwd = 2.0, col = 'red')


plot((res$nb_fxko_cyst + res$nb_pxko_cyst), res$cyst_size ,log = 'xy')
abline(0, 1, lwd = 2.0, col = 'red')

plot(res$pct_nkx22_global, res$pct_nkx22_cyst)
abline(0, 1, lwd =2.0, col = 'red')

plot(res$pct_shh_global, res$pct_shh_cyst)
abline(0, 1, lwd =2.0, col = 'red')


plot(res$genotype_pct_global, res$genotype_pct_cyst)
abline(0, 1, lwd =2.0, col = 'red')

ggplot(res, aes(x=quantile95_pax6, fill=treatment)) +
  geom_histogram(alpha = 0.5, position = 'identity') + 
  #geom_histogram(alpha=0.7, adjust = 1.5) +
  #scale_fill_manual(values=c("darkgreen", "#999999")) + 
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

hist(res$quantile95_pax6[which(res$treatment == 'noRA')], breaks = 40)
hist(res$quantile95_pax6[which(res$treatment != 'noRA')], breaks = 40)
abline(v = 400, lwd = 2.0, col = 'red')

res$pct_shh_cyst[which(res$quantile95_pax6 <400)] = 0

res = res[which(res$treatment != 'noRA'), ]



plot(res$genotype_pct_global, res$genotype_pct_cyst)
abline(0, 1, lwd =2.0, col = 'red')

plot(res$pct_shh_global, res$pct_shh_cyst)
abline(0, 1, lwd =2.0, col = 'red')

plot(res$pct_nkx22_global, res$pct_nkx22_cyst)
abline(0, 1, lwd =2.0, col = 'red')

saveRDS(res, file = paste0(outDir, 'res_TetOnFx_TetOnPx_day6_DVpatterning_sizeFiltering.rds'))


#saveRDS(res, file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))
#res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_sizeFiltering.rds'))
#res = readRDS(file = paste0(outDir, 'res_FXKO_PXKO_day6_DVpatterning_all_global.rds'))
res = readRDS(file = paste0(outDir, 'res_TetOnFx_TetOnPx_day6_DVpatterning_sizeFiltering.rds'))


plot(res$cutoff_pxko - res$cutoff_pxko_image, res$cutoff_fxko - res$cutoff_fxko_image)

diff1 = res$cutoff_fxko - res$cutoff_fxko_image
diff2 = res$cutoff_pxko - res$cutoff_pxko_image

plot(diff1, diff2)

hist(diff1, breaks = 50)
abline(v =c(0.7, -0.5), col = 'red')

jj = which(diff1 > (-0.5) & diff1 < 0.7)
res = res[jj, ]


hist(diff2, breaks = 50)
abline(v =c(0.3, -0.2), col = 'red')

jj = which(diff2 > (-0.2) & diff2 < 0.3)
res = res[jj, ]

saveRDS(res, file = paste0(outDir, 'res_TetOnFx_TetOnPx_day6_DVpatterning_sizeFiltering_genotypeFiltering.rds'))

res = readRDS(file = paste0(outDir, 'res_TetOnFx_TetOnPx_day6_DVpatterning_sizeFiltering_genotypeFiltering.rds'))

diff3 = res$cutoff_foxa2 - res$cutoff_foxa2_image
res = res[which(diff3 > (-10) & diff3 < 10), ]

diff3 = res$cutoff_foxa2 - res$cutoff_foxa2_image
hist(diff3, breaks = 40)

hist(res$cutoff_pax6, breaks = 40)
abline(v = 0.8)

aa = res[which(res$cutoff_pax6 < 0.8), ]

ggplot(aa, aes(x=genotype_pct_cyst, y=pct_nkx22_cyst, color=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=loess, aes(fill=treatment))+
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

res$pct_shh_cyst[which(res$quantile95_pax6 <300)] = 0
aa = res[which(res$treatment != 'noRA' & res$quantile95_pax6 >400), ]
ggplot(aa, aes(x=genotype_pct_cyst2, y=pct_shh_cyst, color=treatment)) +
  geom_point(size = 2.5) + 
  geom_smooth(method=loess, aes(fill=treatment))+
  ylab("% marker+ ") + 
  xlab("genotype % TetOn-FoxA2") + 
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


# res = res[which(res$cutoff_foxa2 > 1.5 & res$cutoff_foxa2 < 8), ]
# 
# plot(res$genotype_pct_pxko, res$pct_foxa2)
# abline(0, 1, lwd = 2.0, col = 'red')
# 
# hist(res$cutoff_marker, breaks = 40)
# abline(v = c(0.25, 1.2), lwd = 2.0, col = 'red')
# res = res[which(res$cutoff_marker < 2.0), ]
# 
# plot(res$genotype_pct_pxko, res$pct_foxa2)
# abline(0, 1, lwd = 2.0, col = 'red')
# 
# 
# hist(res$cutoff_fxko, breaks = 40)
# abline(v = c(0.2, 2), lwd = 2.0, col = 'red')
# 
# res = res[which(res$cutoff_fxko > 0.2 & res$cutoff_fxko < 2.0), ]
# 
# plot(res$genotype_pct_pxko, res$pct_foxa2)
# abline(0, 1, lwd = 2.0, col = 'red')
# 
# hist(res$cutoff_pxko, breaks = 40)
# abline(v = c(0.2, 4.2), lwd = 2.0, col = 'red')
# 
# plot(res$cutoff_fxko, res$cutoff_pxko)
# res = res[which(res$cutoff_pxko < 4.2), ]
# 
# res$pct_foxa2_false = res$nb_foxa2_fxko/res$nb_fxko
# 
# hist(res$pct_foxa2_false, breaks = 50)
# abline(v = c(0.1, 0.05), lwd = 2.0, col = 'red')
# #hist(res$cutoff_fxko, breaks = 20)
# 
# res = res[which(res$pct_foxa2_false < 0.05), ]
# 
# plot(res$genotype_pct_pxko, res$pct_foxa2)
# abline(0, 1, lwd = 2.0, col = 'red')
# 
# ## plot the pct FoxA2 in function of genotype pct
# ggplot(res, aes(x=genotype_pct_cyst, y=pct_shh_cyst, color = )) +
#   geom_point() + 
#   geom_smooth(method=lm) +
#   ylab("% FoxA2+ ") + 
#   xlab("genotype % TetOn-FoxA2") + 
#   theme_bw() +  
#   theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
#         axis.text.y = element_text(angle = 0, size = 12)) +
#   theme(legend.key = element_blank()) + 
#   theme(plot.margin=unit(c(1,3,1,1),"cm"))+
#   #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
#   theme(legend.title = element_blank(), 
#         legend.text = element_text(size = 14))
# 
# ggsave(filename = paste0(outDir, 'RA_WT_KO_KO_day6_DVpatterning_FoxA2pct.pdf'), height = 6, width = 10)

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



aa = res

# Add regression lines
ggplot(aa, aes(x=genotype_pct_global, y=pct_shh_global, color=treatment)) +
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

ggplot(aa, aes(x=genotype_pct_cyst, y=pct_shh_cyst, color=treatment)) +
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





