#write all supplemental data tables

if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else if(file.exists('init.R')){
  source('init.R')
}else{
  setwd('..')
  source('init.R')
}
  
write.csv(my_fe, file = 'supp_logFE.csv', quote = F)


get_uniqSets = function(a, b){
  my_fc = matrix(F, nrow = nrow(my_fe), ncol = 0)
  my_fc = cbind(my_fc, my_fe[,b] - (my_fe[,a]))
  my_fc = cbind(my_fc, my_fe[,b+3] - (my_fe[,a+3] ))
  
  p1 = paste(colnames(my_fe)[b], '/', colnames(my_fe)[a])
  p2 = paste(colnames(my_fe)[b+3], '/', colnames(my_fe)[a+3])
  
  colnames(my_fc) = c(p1, p2)
  
  return(my_fc)
}

fc_10a_v_7 = get_uniqSets(1, 2)
fc_10a_v_231 = get_uniqSets(1, 3)
fc_7_v_231 = get_uniqSets(2, 3)
write.csv(cbind(fc_10a_v_7, fc_10a_v_231, fc_7_v_231), file = 'supp_pairwise_logFC.csv', quote = F)
