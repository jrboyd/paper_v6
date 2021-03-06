#sorts data for each mark and plots alongside corresponding mark

if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else
  source('init.R')

keep = apply(my_fe, 1, max) > detect_thresh
my_fe = my_fe[keep,]

lines = c('MCF10A', 'MCF7', 'MDA231')
mods = c('H3K4ME3', 'H3K4AC')
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

#do some spiffy line plots
pdf(paste('Figure_3b_ranked_lines_JB-',date2fname(),'.pdf', sep = ''), width = 10, height = 5)
mat_lay = matrix(1:6, nrow = 3, ncol = 2, byrow = T)
mat_lay = cbind((1:3), mat_lay+3)#line labels
mat_lay = rbind((c(-1,1,2)), mat_lay+2)#mark labels
mat_lay = rbind(c(0,1,1), mat_lay+1)#mark labels

nf = layout(mat_lay, widths = c(.2,1,1), heights = c(.5,.2,1,1,1))
par(mai = c(.4,.5,.1,.1))
textPlot = function(txt){
  plot(c(0,1), c(0,1), type = 'n', axes = F)
  text(.5,.5, txt)
}
markColors = c('#31a354', '#3182bd')
names(markColors) = mods
forePlot = function(x, m){
  #colors = markColors[m]
  lines(x, col = 'black', lwd = 3)
}
toptops = list()
tops = list()

bgPlot = function(x, m){
  x = movingAverage(x, n = 3, centered = T)
  colors = markColors[m]
  plot(x, ylim = c(0,8), type = 'l', lwd = 1, col = colors, xlab = '', ylab = '')
  perc_95 = length(x)/25*5
  perc_97 = length(x)/25*3
  perc_99 = length(x)/25*1
  perc_98 = length(x)/25*2
  perc_kept = perc_98
  lines(c(perc_kept, perc_kept),c(0,6), lty = 3)
  toptop = names(x)[1:perc_kept]
  toptops[[length(toptops) + 1]] <<- ensg2sym[toptop]
  top = names(x)
  tops[[length(tops) + 1]] <<- ensg2sym[top]
  #print(ensg2sym[top])
}
par(mai = rep(0,4))
textPlot(paste('The most highly methylated promoters are also highly acetylated.\nThe most highly acetylated promoters are not necessarily as highly methylated.'))
for(m in mods){
  textPlot(paste(m, 'over 75th percentile sorted by', m,'\n',other_mod(m),'in background'))
}
for(l in lines){
  textPlot(l)
}
par(mai = c(.4,.5,.1,.1))
for(l in lines){
  ac_dat = my_fe[, paste(l,'H3K4AC')]
  me_dat = my_fe[, paste(l,'H3K4ME3')]
  q_thresh = .75
  me_thresh = quantile(me_dat, q_thresh)
  keep = me_dat > me_thresh
  ac_dat = ac_dat[keep]
  me_dat = me_dat[keep]
  
  o = order(me_dat, decreasing = T)
  bgPlot(ac_dat[o], mods[1])
  forePlot(me_dat[o], mods[2])
  
  ac_dat = my_fe[, paste(l,'H3K4AC')]
  me_dat = my_fe[, paste(l,'H3K4ME3')]
  
  
  ac_thresh = quantile(ac_dat, q_thresh)
  
  keep = ac_dat > ac_thresh
  ac_dat = ac_dat[keep]
  me_dat = me_dat[keep]
  
  o = order(ac_dat, decreasing = T)
  bgPlot(me_dat[o], mods[2])
  forePlot(ac_dat[o], mods[1])
}
dev.off()
thresh98 = apply(my_fe, 2, function(x)quantile(x, .98))
thresh75 = apply(my_fe, 2, function(x)quantile(x, .75))
pass98 = ifelse(my_fe, F, F)
for(i in 1:ncol(pass98)){
  pass98[,i] = my_fe[,i] >= thresh98[i]
}
pass75 = ifelse(my_fe, F, F)
for(i in 1:ncol(pass98)){
  pass75[,i] = my_fe[,i] >= thresh75[i]
}

library(limma)
library(venneuler)
#pdf('percentile venns.pdf')
pdf(paste('Figure_3b_ranked_venns_JB-',date2fname(),'.pdf', sep = ''))
layout(matrix(1:4, nrow = 2))
par(mai = c(0,0,1,0), xpd = NA)

perc_vennDiagram = function(x){
  tmp = vennCounts(x)
  tmp[1,4] = 0
  tmp[,4] = format(round(tmp[,4] / sum(tmp[,4]), digits = 2), nsmall = 2)
  tmp[1,4] = ""
  vennDiagram(tmp, names = lines)
}

perc_vennDiagram(pass98[,1:3])
title('98th percentile, H3K4me3')
perc_vennDiagram(pass75[,1:3])
title('75th percentile, H3K4me3')

perc_vennDiagram(pass98[,4:6])
title('98th percentile, H3K4ac')
perc_vennDiagram(pass75[,4:6])
title('75th percentile, H3K4ac')

vennDiagram(pass98[,1:3], names = lines, )
title('98th percentile, H3K4me3')
vennDiagram(pass75[,1:3], names = lines)
title('75th percentile, H3K4me3')

vennDiagram(pass98[,4:6], names = lines)
title('98th percentile, H3K4ac')
vennDiagram(pass75[,4:6], names = lines)
title('75th percentile, H3K4ac')

plot(venneuler(pass98[,1:3]))
plot(venneuler(pass75[,1:3]))

plot(venneuler(pass98[,4:6]))
plot(venneuler(pass75[,4:6]))
dev.off()

#output full lists of genes in top 2%
final = matrix('', nrow = max(colSums(pass98)), ncol = ncol(pass98))
colnames(final) = colnames(pass98)
for(i in 1:ncol(pass98)){
  cname = colnames(pass98)[i]
  keep = pass98[,i]
  vals = sort(my_fe[keep,i], decreasing = T)
  final[1:length(vals),i] = ensg2sym[names(vals)]
}
write.table(final, file =  paste0('Figure_3b_ranked_full_lists_JB-',date2fname(), '.csv' ), sep = ',', row.names = F, col.names = T, quote = F)

#output uniques
final_uniq = matrix('', nrow = nrow(final), ncol = 1.5 * ncol(final))

layout(1:3)
rng = 1:3
for(i in rng){#TERRIBLE CODE, write lists for mark overlaps
  tmp_me = final[,i]
  tmp_ac = final[,i+3]
  tmp_all = union(tmp_me, tmp_ac)
  membs = matrix(F, ncol = 2, nrow = length(tmp_all))
  rownames(membs) = tmp_all
  membs[tmp_me,1] = T
  membs[tmp_ac,2] = T
  vennDiagram(membs)
  me_only = rownames(membs)[membs[,1] & !membs[,2]]
  ac_only = rownames(membs)[!membs[,1] & membs[,2]]
  both = rownames(membs)[membs[,1] & membs[,2]]
  
  final_uniq[1:length(me_only),1+(i-1)*3] = me_only
  final_uniq[1:length(both),2+(i-1)*3] = both
  final_uniq[1:length(ac_only),3+(i-1)*3] = ac_only
  
#   is_uniq = rep(T, length(tmp))
#   for(j in 1:length(tmp)){
#     if(sum(tmp[j] == final[,rng]) > 1){
#       is_uniq[j] = F
#     }
#   }
#   kept = tmp[is_uniq]
#   final_uniq[1:length(kept),i] = kept
}

colnames(final_uniq) = unlist(lapply(lines, function(x)paste(x, c('H3K4me3 top 2%', 'both top 2%', 'H3K4AC top 2%'))))
write.table(final_uniq, file =  paste0('Figure_3b_ranked_unique_lists_JB-',date2fname(), '.csv' ), sep = ',', row.names = F, col.names = T, quote = F)
