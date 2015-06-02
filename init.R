if(exists('parent_dir'))
  setwd(parent_dir)
library(ggplot2)
library(RColorBrewer)
source('jrb_R_scripts/read_countDir.R')
source('jrb_R_scripts/norm_counts.R')
source('jrb_R_scripts/jrb_matrix.R')
source('jrb_R_scripts/heatmap.3-split.R')
source('functions_misc.R')
source('functions_ngsplot.R')
source('functions_movingAverage.R')
load('ref//dictionaries.save')
load('ref//enst_dicts.save')
load('ref//gsea_dataset.save')
load('ref//gsea_membership.save')

in_parent_dir = function(file){
  return(paste(parent_dir, file, sep = '/'))
}


if(!exists('parent_dir'))
  parent_dir = getwd()
out_dir = paste0('figures_',date2fname())
dir.create(out_dir, showWarnings = F)
setwd(paste(parent_dir, out_dir, sep = '/'))

if(!file.exists('my_fe_corrected.save')){
  source(in_parent_dir('process_mycounts.R'))
}else
  load('my_fe_corrected.save')


my_fe = my_fe[,c(4:6,1:3)]
detect_thresh = 2 #this threshold is used repeatedly for FE/FC thresholding
gsea_sig_thresh = 12 #-log10 pvalue threshold for gsea testing
primeColors = rep(RColorBrewer::brewer.pal(n = 3, 'Dark2'),2)

