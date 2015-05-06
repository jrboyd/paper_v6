##produce heatmaps of selected gsea lists based on gseasigmaps
if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else if(file.exists('init.R')){
  source('init.R')
}else{
  setwd('..')
  source('init.R')
}
load('pw_uniquely.save')
load('pw_uniq2gsea_results.save')




plot_GSEA_heatmaps = function(res){
  pthresh = 9
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
  
  
#   keep = sapply(gsea_passing, function(x)return(any(x == listNames)))
#   
#   gsea_passing = gsea_passing[keep]
#   uniq_passing = uniq_passing[keep]
  membership = gsea_master[,gsea_passing]
  
  for(i in 1:length(gsea_passing)){
    g = gsea_passing[i]
    u = uniq_passing[i]
    gsea_list = rownames(membership)[membership[,g]]
    uniq_list = rownames(uniquely_K4_membership)[uniquely_K4_membership[,u]]
    is_uniq = intersect(gsea_list, uniq_list)
    not_uniq = sapply(gsea_list, function(x)return(!any(x == is_uniq)))
    not_uniq = names(not_uniq)[not_uniq]
    
    keep = sapply(ensg2sym[rownames(my_fe)], function(x)return(any(x == is_uniq)))
    is_uniq_fe = my_fe[keep,]
    hres = heatmap.3(is_uniq_fe / rowSums(is_uniq_fe), nsplits = 2, classCount = 3, main = paste(g, u, sep = '\n'))
    tmp = hres[[3]]
    rownames(tmp) = ensg2sym[rownames(tmp)]
    hres[[3]] = tmp
    plot.HeatmapLists(hres)
    
    keep = sapply(ensg2sym[rownames(my_fe)], function(x)return(any(x == not_uniq)))
    not_uniq_fe = my_fe[keep,]
    hres = heatmap.3(not_uniq_fe / rowSums(not_uniq_fe), nsplits = 2, main = paste(g, paste('NOT', u), sep = '\n'))
    tmp = hres[[3]]
    rownames(tmp) = ensg2sym[rownames(tmp)]
    hres[[3]] = tmp
    plot.HeatmapLists(hres)
    
  }
}
out_dir = paste('Figure_5e_gsea_heatmaps_JB-',date2fname(), sep = '')
dir.create(out_dir)
setwd(out_dir)
pdf('MCF10A vs MCF7.pdf')
plot_GSEA_heatmaps(res_10a_v_7)
dev.off()

pdf('MCF10A vs MDA231.pdf')
plot_GSEA_heatmaps(res_10a_v_231)
dev.off()

pdf('MCF7 vs MDA231.pdf')
plot_GSEA_heatmaps(res_7_v_231)
dev.off()
setwd('..')
# 
# listNames = c('GOZGIT_ESR1_TARGETS_DN',
#               'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP',
#               'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN',
#               'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP',
#               'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN',
#               'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_UP',
#               'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN')
# 
# expected_good = c('CREIGHTON_ENDOCRINE_THERAPY_RESISTANCE_1',
#                   'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP',
#                   'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN',
#                   'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN',
#                   'HUANG_DASATINIB_RESISTANCE_UP'
# )
# # pdf('examples.pdf')
# my_fe = ifelse(my_fe > 0, my_fe, 0)
# my_fe[,1:3] = my_fe[,1:3] / (rowSums( my_fe[,1:3])+1)
# my_fe[,4:6] = my_fe[,4:6] / (rowSums( my_fe[,4:6])+1)
# for(ln in listNames){
#   geneList = rownames(gsea_master)[gsea_master[,ln]]
#   keep = sapply(ensg2sym[rownames(my_fe)], function(x)return(any(x == geneList)))
#   list_fe = my_fe[keep,]
#   heatmap.3(list_fe, dendrogram = 'n', trace = 'n', Colv = F, scale = 'n', classCount = 5, clusterSort = F, nsplits = 2, main = ln, cex.main = .5)
#   
# #   pw_uniq = uniq_7_v_231
# #   non_uniq = !sapply(ensg2sym[rownames(list_fe)], function(x)return(any(x == rownames(pw_uniq))))
# #   toPlot = matrix(0, nrow = 0, ncol = ncol(list_fe))
# #   rowSeps = numeric()
# #   for(i in 1:ncol(pw_uniq)){
# #     uniq = rownames(pw_uniq)[pw_uniq[,i]]
# #     keep = sapply(ensg2sym[rownames(list_fe)], function(x)return(any(x == uniq)))
# #     keep = names(keep)[keep]
# #     dat = list_fe[keep,,drop = F]
# # #     key = pw_key[[i]]
# # #     key_sub = pw_key_sub[[i]]
# #     if(nrow(dat) > 2){
# #       clust = hclust(dist(dat))
# #       o = clust$order
# #       dat = dat[o,]
# #     }
# #     rowSeps = c(rowSeps, max(rowSeps, 0)+nrow(dat))
# #     toPlot = rbind(toPlot, dat)
# #   }
# #   
# #   heatmap.2(toPlot, dendrogram = 'n', trace = 'n', cexCol = .8, Rowv = F, rowsep = rowSeps, sepcolor = 'black', Colv = F, scale = 'n', margins = c(10,10), main = ln)
# #   heatmap.2(list_fe[non_uniq,], dendrogram = 'n', trace = 'n', cexCol = .8, Rowv = T, rowsep = rowSeps, sepcolor = 'black', Colv = F, scale = 'n', margins = c(10,10), main = ln)
# }
# dev.off()
