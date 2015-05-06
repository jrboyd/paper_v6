#executes all figure generation code

if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else if(file.exists('init.R')){
  source('init.R')
}else{
  setwd('..')
  source('init.R')
}
setwd('..')

for(f in dir(pattern = 'pipeline_fig.+\\.R')){
  print(f)
  source(in_parent_dir(f))
}
