#starting with uniquely sets, tests them for enrichment using binom.test in gsea lists.
#a heatmap of -log10 pvals is written to pdf (gsea lists are rows, uniq groups are columns)
#gsea lists passing a threshold are added to gsea_passing.save
#the corresponding best uniquely group is added to uniq_passing.save

##dependencies
if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else if(file.exists('init.R')){
  source('init.R')
}else{
  setwd('..')
  source('init.R')
}
load('pw_uniq2gsea_results.save')

##parameters
#name of heatmap pdf
out_dir = paste('Figure_5d_list_overlaps_JB-',date2fname(),sep = '')
dir.create(out_dir)
pdfName = paste0(out_dir, '/gsea_PAIR_important_genes.pdf')
savePlot = T
#starting pvalue - will be reduced until multiple gsea lists pass
pthresh = 9
#bg size, all detected?  all passing FE?  all genes?
bg_size = 20000
#gene most be present in at least this many gsea lists
min_shared = 1






get_important_genes = function(res, pdfName){
  if(!file.exists(dirname(pdfName)))
    dir.create(dirname(pdfName))
  gsea_membs = gsea_membership[,res[[1]]]
  keep = rowSums(gsea_membs) > 0
  gsea_membs = gsea_membs[keep,]
  gsea_membs = ifelse(gsea_membs, 1, 0)
  
  uniquely_K4_membership = res[[2]]
  
  pvals = res[[3]]
  
  
  if(savePlot) pdf(pdfName, width = 12, height = 8)
  for(i in 1:ncol(uniquely_K4_membership)){
    grp_pvals = pvals[,colnames(uniquely_K4_membership)[i]]
    keep = names(grp_pvals)[grp_pvals > 9]
    gsea_membership_sig = gsea_membs[,keep, drop = F]
    keep = rowSums(gsea_membership_sig) > 0
    gsea_membership_sig = gsea_membership_sig[keep,, drop = F]
    uniq_in_gsea_sig = uniquely_K4_membership[intersect(rownames(gsea_membership_sig), rownames(uniquely_K4_membership)),]
    uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
    toPlot = t(ifelse(gsea_membership_sig[uniq_members,, drop = F], 1, 0))
    if(min(dim(toPlot)) < 1)
      next
    print(dim(toPlot))
    r_order = order(rowSums(toPlot), decreasing = T)
    c_order = order(colSums(toPlot), decreasing = T)
    toPlot = toPlot[r_order, c_order, drop = F]
    
    if(nrow(toPlot) > 1){
      j = 1
      tot = 0
      end = ncol(toPlot)
      while(tot < ncol(toPlot) & j <= nrow(toPlot)){
        start = tot + 1
        o = order(toPlot[j,start:end], decreasing = T)
        n = sum(toPlot[j,start:end])
        #   if(n == 0)
        #     next
        
        o = o + start - 1
        tot = tot + n 
        toPlot[,start:end] = toPlot[,o, drop = F]
        
        j = j + 1
        if(j > nrow(toPlot))
          next
        o = order(toPlot[j,start:tot], decreasing = F)
        o = o + min(start, tot) - 1
        toPlot[,start:tot] = toPlot[,o, drop = F]
      }
      
      heatmap.2((toPlot), scale = 'n', margins = c(12,16), cexCol = min(20/ncol(toPlot),.4), cexRow = min(10/nrow(toPlot), .6),Colv = F, Rowv = F, col = c('white', 'black'), trace = 'n', dendrogram = 'n', main = colnames(uniq_in_gsea_sig)[i], key = F)
    }}
  if(savePlot) dev.off()
  
  output = character(length = ncol(uniq_in_gsea_sig))
  names(output) = colnames(uniq_in_gsea_sig)
  output_shared = output
  
  for(i in 1:ncol(uniquely_K4_membership)){
    grp_pvals = pvals[,colnames(uniquely_K4_membership)[i], drop = F]
    keep = rownames(grp_pvals)[grp_pvals > 9]
    if(is.null(keep)){
      next
    }
    gsea_membership_sig = gsea_membs[,keep, drop = F]
    keep = rowSums(gsea_membership_sig) > 0
    gsea_membership_sig = gsea_membership_sig[keep,, drop = F]
    uniq_in_gsea_sig = uniquely_K4_membership[intersect(rownames(gsea_membership_sig), rownames(uniquely_K4_membership)),]
    uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
    output[i] = paste(uniq_members, collapse = ',')
    
    keep = rownames(grp_pvals)[grp_pvals > 9]
    gsea_membership_sig = gsea_membs[,keep, drop = F]
    keep = rowSums(gsea_membership_sig) > 1
    gsea_membership_sig = gsea_membership_sig[keep,, drop = F]
    uniq_in_gsea_sig = uniquely_K4_membership[intersect(rownames(gsea_membership_sig), rownames(uniquely_K4_membership)),]
    uniq_members = rownames(uniq_in_gsea_sig)[uniq_in_gsea_sig[,i]]
    output_shared[i] = paste(uniq_members, collapse = ',')
  }
  output = paste(names(output), output, sep = ',', collapse = '\n')
  #write.table(output, file = sub('.pdf', '.csv', pdfName), col.names = F, row.names = F, quote = F)
  
  output_shared = paste(names(output_shared), output_shared, sep = ',', collapse = '\n')
  #write.table(output_shared, file = sub('.pdf', '_shared.csv', pdfName), col.names = F, row.names = F, quote = F)
  return(list(output, output_shared))
}

a = get_important_genes(res_10a_v_7, sub('PAIR', '10a_v_7', pdfName))
b = get_important_genes(res_10a_v_231, sub('PAIR', '10a_v_231', pdfName))
c = get_important_genes(res_7_v_231, sub('PAIR', '7_v_231', pdfName))
