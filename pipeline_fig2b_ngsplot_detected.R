#sorts data for each mark and plots alongside corresponding mark

if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else
  source('init.R')

markData_4me3_4ac = my_fe
keep = apply(markData_4me3_4ac, 1, max) > detect_thresh
markData_4me3_4ac = markData_4me3_4ac[keep,]

lines = c('MCF10A', 'MCF7', 'MDA231')
mods = c('H3K4AC', 'H3K4ME3')
other_mod = function(m){
  if(m == mods[1])
    return(mods[2])
  else if(m == mods[2])
    return(mods[1])
  return('unrecognized mod')
}

applyWindow = function(dat, win = 10){
  if(win < 2)
    return(dat)
  out = matrix(0, nrow = nrow(dat), ncol = ncol(dat)/win)
  for(i in 1:ncol(out)){
    start = (i-1)*win+1
    end = i*win
    #out[,i] = apply(dat[,start:end],1,median)
    out[,i] = rowMeans(dat[,start:end])
  }
  return(out)
}
#pdf('ngs_custom.pdf')
if(!exists('all_profiles')){
  all_profiles = list()
  for(n in sub(' ', '_', colnames(markData_4me3_4ac))){
    
    in_dir = dir(in_parent_dir('ngsplot_data'), full.names = T, pattern = n)
    
    fname = paste(in_dir, '/hm1.txt', sep = "")
    
    tmp = read.table(fname, stringsAsFactors = F)
    ensgs = tmp[2:nrow(tmp),1]
    strand = tmp[2:nrow(tmp),4]
    dat = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
    rownames(dat) = ensgs
    dat = dat[intersect(rownames(dat), ensg2cut[rownames(markData_4me3_4ac)]),]
    all_profiles[[n]] = dat
  }
}

me3_i = 1:3
ac_i = 4:6
xs = 0:100 * 20 - 1000   

pdf(paste('Figure_2b_ngsplot_JB-',date2fname(),'.pdf', sep = ''), width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
ngsplot_3lines = function(indexes, name){
  plot(1, xlim = c(-1000, 1000), ylim = c(0,3.5), type = 'n', xlab = 'bp from TSS', ylab = 'FE')
  title(name)
  for(i in indexes){
    dat = all_profiles[[i]]
    lines(xs, colMeans(dat), col = primeColors[i], lwd = 3)
  }
}
ngsplot_3lines(me3_i, 'H3K4me3')
ngsplot_3lines(ac_i, 'H3K4ac')
legend(x = 'topright', legend = c('MCF10A', 'MCF-7', 'MDA-MB-231'), fill = primeColors, bty = 'n')
dev.off()
