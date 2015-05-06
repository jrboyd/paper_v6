

mycounts = countDir.read(in_parent_dir('counts_mine/'))
tmp = read.table(in_parent_dir('million_mapped_reads.txt'), stringsAsFactors = F)
mm_reads = tmp[,2]
names(mm_reads) = sapply(strsplit(tmp[,1], '_'), function(x)return(paste(x[1:3], collapse = ' ')))

tmp = jrb.split_colnames(mycounts)
cell_lines = tmp[1,]
histone_mods = tmp[2,]
reps = tmp[3,]

colnames(mycounts) = paste(cell_lines, histone_mods, reps)

mycpm = mycounts

for(cn in colnames(mycpm)){
  mycpm[,cn] = mycpm[,cn] / mm_reads[cn]
}

myrpkm = mycpm / 2 # since all promoters are 2kb

myrpkm_agg = matrix(0, nrow = nrow(myrpkm), ncol = 0)
rownames(myrpkm_agg) = rownames(myrpkm)
for(cl in unique(cell_lines)){
  for(hm in unique(histone_mods)){
    match = cell_lines == cl & histone_mods == hm
    myrpkm_agg = cbind(myrpkm_agg, rowMeans(myrpkm[,match]))
    colnames(myrpkm_agg)[ncol(myrpkm_agg)] = paste(cl, hm)                   
  }
}

floor = quantile(myrpkm_agg, .05)
myrpkm_aggfloor = ifelse(myrpkm_agg < floor, floor, myrpkm_agg)

myfe = norm.fe(myrpkm_aggfloor)
keep = apply(myfe, 1, max) > 1
myfe = myfe[keep,]
layout(matrix(1:3, nrow = 1))
plot(log2(myfe[,2:1]), xlim = c(0,6), ylim = c(0,6))
plot(log2(myfe[,4:3]), xlim = c(0,6), ylim = c(0,6))
plot(log2(myfe[,6:5]), xlim = c(0,6), ylim = c(0,6))
myfe = myfe[,c(1,3,5,2,4,6)]
myfe = log2(myfe)
library(limma)
layout(1)
boxplot(myfe)
vennDiagram(myfe[,1:3] > 2)
vennDiagram(myfe[,4:6] > 2)
my_fe = myfe

save(my_fe, file = 'my_fe_corrected.save')
