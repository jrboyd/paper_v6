if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else
  source('init.R')
figGen = T
markData_4me3_4ac = my_fe
isMarked = markData_4me3_4ac > detect_thresh
# enst_isMarked = markData_4me3_4ac > 1
# uniq_ensg = unique(enst2ensg)
# ensg_isMarked = 

library(limma)
cellLines = c("MCF-10a", "MCF-7", "MDA-MB-231")

fname = paste0("Figure_2c_vennDiagrams_JB_",date2fname(),".pdf")
if(figGen)
  pdf(fname, width = 16, height = 8.5)
layout(matrix(1:2, ncol = 2))
vennDiagram(isMarked[, 1:3], names = cellLines)
title("H3K4me3")

vennDiagram(isMarked[, 4:6], names = cellLines)
title("H3K4ac")
if(figGen)
  dev.off()
print(paste("wrote venn diagrams to", fname)) 

q_uniques = list(
  c(T,F,F),
  c(F,T,F),
  c(F,F,T))
testRow = function(x)return(all(x == q))

# get_list = function(q, me_OR_ac){
#   key = 1:3
#   if(me_OR_ac == 'me'){
#     key = 4:6
#   }
#   keep = apply(isMarked[, key], 1, testRow)
#   ensgs = rownames(isMarked)[keep]
#   keep = !is.na(ensg2sym[enst2ensg[ensgs]])
#   if(sum(keep) != length(ensgs))
#     warning('not all ensgs mapped to symbol')
#   ensgs = ensgs[keep]
#   return(ensg2sym[enst2ensg[ensgs]])
# }
# 
# for(q in q_uniques){
#   ac_list = get_list(q, 'ac')
#   fname = paste('vennList uniquely ', colnames(isMarked)[1:3][q], '.txt', sep = '')
#   write.table(ac_list, file = fname, row.names = F, col.names = F, quote = F)
#   me_list = get_list(q, 'me')
#   fname = paste('vennList uniquely ', colnames(isMarked)[4:6][q], '.txt', sep = '')
#   write.table(me_list, file = fname, row.names = F, col.names = F, quote = F)
# }
# 
# q_pairs = lapply(q_uniques, function(x)return(!x))
# 
# for(q in q_pairs){
#   ac_list = get_list(q, 'ac')
#   fname = paste('vennList double ', 
#                 paste(colnames(isMarked)[1:3][q], collapse = '-'),
#                 '.txt', sep = '')
#   write.table(ac_list, file = fname, row.names = F, col.names = F, quote = F)
#   me_list = get_list(q, 'me')
#   fname = paste('vennList double ', 
#                 paste(colnames(isMarked)[4:6][q], collapse = '-'),
#                 '.txt', sep = '')
#   write.table(me_list, file = fname, row.names = F, col.names = F, quote = F)
# }
# 
# q_universal = c(T,T,T)
# q = q_universal
# ac_list = get_list(q, 'ac')
# fname = paste('vennList universal K4AC.txt')
# write.table(ac_list, file = fname, row.names = F, col.names = F, quote = F)
# me_list = get_list(q, 'me')
# fname = paste('vennList universal K4me3.txt')
# write.table(me_list, file = fname, row.names = F, col.names = F, quote = F)
# 
# q_never = c(F,F,F)
# q = q_never
# ac_list = get_list(q, 'ac')
# fname = paste('vennList never K4AC.txt')
# write.table(ac_list, file = fname, row.names = F, col.names = F, quote = F)
# me_list = get_list(q, 'me')
# fname = paste('vennList never K4me3.txt')
# write.table(me_list, file = fname, row.names = F, col.names = F, quote = F)
# 
# print('wrote all venn lists')
