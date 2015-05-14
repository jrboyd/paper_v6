sym2dat = function(sym){
  ensg = sapply(sym, function(x){
    return(names(ensg2sym)[ensg2sym == x])
  })
  print(my_fe[ensg,])
  #print(mycounts[ensg,])
  return(my_fe[ensg,])
}

date2fname = function(){
  d = date()
  tmp = strsplit(d, ' ')[[1]]
  d = paste(tmp[c(2:3,5)],collapse = '-')
  return(d)
}
