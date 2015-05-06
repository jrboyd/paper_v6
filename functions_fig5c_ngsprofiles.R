#outputs all ngs profiles for gsea lists in gsea_passing.save.
#separately plots groups in uniq_passing.save, highlighting whether mark trend hold for full gsea list.

##dependencies


markData_4me3_4ac = my_fe
load('pw_uniq2gsea_results.save')
load('pw_uniquely.save')

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

# cut2ensg = names(ensg2cut)
# names(cut2ensg) = ensg2cut
# ngs_ensg = cut2ensg[rownames(all_profiles[[1]])]
# NGS_SYMBOL_SET = unique(ensg2sym[ngs_ensg])
NGS_SYMBOL_SET = rownames(all_profiles[[1]])

lines = c('MCF10A', 'MCF7', 'MDA231')
library(RColorBrewer)
l2col = RColorBrewer::brewer.pal(n = 3, 'Dark2')
names(l2col) = lines
l2col_bg = RColorBrewer::brewer.pal(n = 3, 'Set2')
names(l2col_bg) = lines
mods = c('H3K4AC', 'H3K4ME3')
other_mod = function(m){
  if(m == mods[1])
    return(mods[2])
  else if(m == mods[2])
    return(mods[1])
  return('unrecognized mod')
}

ngsplot_passing = function(res, line_a, line_b, pthresh = 9, outdir_name = 'output_pw'){
  if(!file.exists(outdir_name))
    dir.create(outdir_name)
  grp_name = paste(line_a,'vs', line_b, sep = '_')
  dir_name = paste(outdir_name, '//ngs_', grp_name, sep = '')
  if(!file.exists(dir_name))
    dir.create(dir_name)
  
  
  passing = res[[3]] > pthresh
  uniq_passing = unlist(apply(passing, 1, function(x)return(colnames(passing)[x])))
  ngsea_passing = apply(passing, 1, function(x)return(sum(x)))
  gsea_passing = unlist(sapply(names(ngsea_passing), function(x){
    npassing = ngsea_passing[x]
    out = rep(x, npassing)
    names(out) = NULL
    return(out)
  }) )
  names(gsea_passing) = NULL
  names(uniq_passing) = NULL
  uniquely_K4_membership = res[[2]]
  membership = gsea_master[,gsea_passing]
  all_fg = character()
  for(j in 1:length(gsea_passing)){
    list_name = gsea_passing[j]
    best_uniq = uniq_passing[j]
    fname = paste(dir_name, '//', best_uniq,'_', list_name, '.pdf', sep = '')
    pdf(fname, width = 10, height = 8)
    toPlot = membership[,list_name]
    toPlot = names(toPlot)[toPlot] #character vector of gene symbols in list
    
    isSigUniq = uniquely_K4_membership[,best_uniq]
    uniq_name = sub('_uniqAConly', ' uniquely marked by ac, not me3', best_uniq)
    uniq_name = sub('_uniqMEonly', ' uniquely marked by me3, not ac', uniq_name)
    uniq_name = sub('_both', ' uniquely marked by both me3 and ac', uniq_name)
    if(sum(isSigUniq) == 0){
      next
    }
    isUniq = uniquely_K4_membership[intersect(toPlot, rownames(uniquely_K4_membership[isSigUniq,])),,drop = F]
    fg_toPlot = rownames(isUniq) #those genes that are uniq in a cell line
    isBg = rep(T, length(toPlot))
    names(isBg) = toPlot
    isBg[fg_toPlot] = F
    bg_toPlot = toPlot[isBg]#those genes that are not uniq in a cell line
    all_fg = union(all_fg, fg_toPlot)
    fg_toPlot = sym2cut(fg_toPlot)
    bg_toPlot = sym2cut(bg_toPlot)
    plotNGS_wBG(fg_toPlot, bg_toPlot, list_name, uniq_name, ymax = 5, linesToPlot = c(line_a, line_b))
    dev.off()
  }
  return(all_fg)
}