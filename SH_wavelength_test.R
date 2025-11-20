rm(list =ls())

require(ggplot2)
library(plotly)
library(ggrepel)
library(tidyverse)
library(patchwork)

dataDir = "res/test_organoid_KOKO_TetOnTetOn_timepoints_v1"

figureDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/",
                   "projects/RA_competence/results/figures_tables_R13547_10x_mNT_20240522/SHanalysis_test")

if(!dir.exists(figureDir)){
  dir.create(figureDir)
}

########################################################
########################################################
# Section I : # Test wavelength parameter to regress out the organoid size
# 
########################################################
########################################################

ExploraryAnalaysis_testWavelength = FALSE
if(ExploraryAnalaysis_testWavelength)
{
  res = readRDS(file = paste0(dataDir, 
                              '/saved_res_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered_FoxA2SDfiltered.rds'))
  hist(log10(res$cyst_size))
  
  #res = res[which(res$time != "d3.5" & res$time != "d5"), ]
  res$r2 = res$cyst_r^2
  
  cc.levels = levels = c("WT", "KO_KO", "TetOn_TetON_RA", 'TetOn_TetON_dox')
  res$condition = factor(res$condition, levels = cc.levels)
  
  jj1 = which(res$condition == 'WT')
  jj2 = which(res$condition == 'KO_KO')
  jj3 = which(res$condition == 'TetOn_TetON_RA')
  jj4 = which(res$condition == 'TetOn_TetON_dox')
  
  jj = which(res$time == 'd3' & res$condition != 'TetOn_TetON_dox')
  
  res$intensity_mean_foxa2_log = log10(res$intensity_mean_foxa2)
  res$cv_log = log10(res$intensity_std_foxa2/res$intensity_mean_foxa2)
  res$cv = (res$intensity_std_foxa2/res$intensity_mean_foxa2)
  res$cc = paste0(res$condition, '_', res$time)
  
  ##########################################
  # search for the size-invariance parameter based on FoxA2 distribution
  ##########################################
  ggplot(res[c(jj1, jj3), ],  aes(x=cyst_size_log, y=intensity_std_foxa2_log, color = time, shape = condition)) + 
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
  
  
  
  ggplot(res[c(jj1, jj3), ],  aes(x=intensity_mean_foxa2_log, 
                                  y=intensity_std_foxa2_log, color = time, shape = condition)) + 
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
  

  ggplot(res[c(jj1), ],  aes(x=cyst_size_log, 
                                  y=cv, color = time, shape = condition)) + 
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
  
  ggsave(filename = paste0(figureDir, '/searching_sizeInvariance_metric_CV_WT.pdf'), height = 6, width = 10)
  
  
  
  ggplot(res[c(jj1, jj3), ],  aes(x=cyst_size_log, 
                             y=cv, color = time, shape = condition)) + 
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
  
  ggsave(filename = paste0(figureDir, '/searching_sizeInvariance_metric_CV_WT_TetOnTetONRA.pdf'), height = 6, width = 10)
  
  
  
  ##########################################
  # test the wavelength metric 
  ##########################################
  jj1 = which(res$condition == 'WT')
  jj2 = which(res$condition == 'KO_KO')
  jj3 = which(res$condition == 'TetOn_TetON_RA')
  jj4 = which(res$condition == 'TetOn_TetON_dox')
  
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
  
  #aa = aa[, -1]
  xx = aa
  for(n in 1:nrow(xx))
  {
    xx[n, ] = 2*pi*res$cyst_r[n]*3 /(c(0:100) + 0.5)
  }
  
  plot(c(1,1), type = 'n', 
       xlim = c(50, 1000),
       ylim = c(10^-5, 1), log = 'xy', xlab = 'wavelength', ylab = 'variance %')
  
  nb_sample = 50
  
  index = c(which(res$cc == 'WT_d3')[1:nb_sample], 
            which(res$cc == 'WT_d4')[1:nb_sample],
            which(res$cc == 'WT_d5')[1:nb_sample]
            )
  cols = rep(c('darkred', 'orange', 'darkblue'), each = nb_sample)
  for(n in 1:length(index))
  {
    points(xx[index[n], ], aa[index[n],], col = cols[n], type = 'b')
  }
  
  #matplot(log(t(xx[c(1:100), ])), type = "l")
  lseq <- function(from=1, to=100000, length.out=6) {
    # logarithmic spaced sequence
    # blatantly stolen from library("emdbook"), because need only this
    exp(seq(log(from), log(to), length.out = length.out))
  }
  
  lambda = lseq(20, 500, length.out = 20)
  vv = matrix(0, nrow = nrow(aa), ncol = length(lambda))
  colnames(vv) = paste0('lambda_', c(1:length(lambda)))
  rownames(vv) = rownames(aa)
  
  library('gam')
  library(KernSmooth)
  library(locfit)
  
  for(n in 1:nrow(vv))
  {
    # n = 1
    cat(n, '\n')
    y0 = log10(as.numeric(aa[n,]))
    tt = log10(as.numeric(xx[n,]))
    fit <- gam(y0 ~ s(tt, df=6), family = gaussian) 
    
    prediction = predict.Gam(fit, newdata = data.frame(tt = log10(lambda)))
    
    #plot(tt, y0, log = '')
    #points(log10(lambda), prediction, type='b', col = 'red')
    vv[n, ] = prediction
    
    #fit = locfit(y0~lp(tt, nn = 0.3, h=0.05, deg=2))
    #plot(tt, y0, log = 'xy')
    #lines(fit)
  }
  
  save(lambda, vv, file = paste0(figureDir, '/lambda_wavelength_2saved.Rdata'))
  
  load(file = paste0(figureDir, '/lambda_wavelength_2saved.Rdata'))
  #colnames(vv) = lambda
  vv = data.frame(genotype = res$condition, time = res$time, condition = res$cc, vv)
  
  #saveRDS(vv, file = paste0(dataDir, '/saved_wavelength_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered.rds'))
  
  
  library(reshape2)
  xx = reshape2::melt(vv, id.vars = c(1:3), variable.name = "lambda", value.name = "variance")
  
  xx$lambda = as.numeric(gsub('lambda_', '', as.character(xx$lambda)))
  xx$lambda = log10(lambda)[xx$lambda]
  
  
  library(ggplot2)
  # Basic scatter plot
  kk = which(xx$time != 'd3.5' & xx$time != 'd5')
  ggplot(xx[kk, ], aes(x=lambda, y=variance, color = time, shape = genotype)) + 
    #geom_point(size = 0.0) +
    geom_smooth(data = xx[kk, ], method=loess, aes(x=lambda, y=variance, color = time, linetype = genotype),
                se  = FALSE)+
    #geom_line(aes(y = mean, color = time), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = time), alpha = .2) +
    #xlab("cyst surface (r2 log2)") + 
    #ylab("lmax") + 
    #geom_hline(yintercept=c(2.5, 2.0), linetype="dashed", color = "red") + 
    #geom_vline(xintercept=c(3.75), linetype="dashed", color = "red") + 
    theme_bw()
  
  ggsave(filename = paste0(figureDir, '/searching_sizeInvariance_metric_WavelengthTest.pdf'), height = 6, width = 10)
  
  
  load(file = paste0(figureDir, '/lambda_wavelength_2saved.Rdata'))
  #colnames(vv) = lambda
  
  metadata = res[, c(1:31, (ncol(res)-7):ncol(res))]
  colnames(metadata)[which(colnames(metadata) == 'condition')] = 'genotype'
  colnames(metadata)[which(colnames(metadata) == 'cc')] = 'condition'
  vv = data.frame(vv, metadata)
  
  saveRDS(vv, file = paste0(dataDir, '/saved_wavelength_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered.rds'))
  
  
}


##########################################
# ## prepare the data for singlecellExperiment
##########################################
Test_umap_summary = FALSE
if(Test_umap_summary){
  library(SingleCellExperiment)
  
  res = readRDS(file = paste0(dataDir, '/saved_wavelength_WT_KOKO_TetOnTetOn_metadataFeatureCollection_filtered.rds'))
  
  counts = as.matrix(t(res[, 1:20]))
  mat = res[, -c(1:20)]
  #mat$treatment = 'RA_wt'
  #mat$condition = mat$time
  metadata = mat
  metadata$marker_class = 'type'
  metadata = data.frame(metadata, stringsAsFactors = FALSE)
  
  sce <- SingleCellExperiment(assays=list(counts=counts),
                              colData=metadata, 
                              metadata = metadata)
  
  y <- assay(sce, "counts")
  #y <- asinh(sweep(y, 1, cf, "/"))
  assay(sce, "exprs", FALSE) <- y
  
  sce$sample_id = as.character(sce$genotype)
  #sce$condition = gsub('noRA_d2', "beforeRA_d2", sce$condition)
  
  cc.levels = levels = c("WT", "KO_KO", "TetOn_TetON_RA", 'TetOn_TetON_dox')
  
  sce$condition = factor(sce$condition, levels = cc.levels)
  sce@metadata$condition = sce$condition
  
  rowData(sce)$marker_name = rownames(sce)
  rowData(sce)$channel_name = NULL
  rowData(sce)$marker_class = 'type'
  
  saveRDS(sce, file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/data/", 
                             "image_SHout_dlogl_sce_wt_conditions_metadata_wavelength_v4.rds"))
  
  table(sce$genotype, sce$time)
  
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


