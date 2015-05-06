sym2dat = function(sym){
  ensg = names(ensg2sym)[ensg2sym == sym]
  print(my_fe[ensg,])
  print(mycounts[ensg,])
  return(my_fe[ensg,])
}

date2fname = function(){
  d = date()
  tmp = strsplit(d, ' ')[[1]]
  d = paste(tmp[c(2:3,5)],collapse = '-')
  return(d)
}
