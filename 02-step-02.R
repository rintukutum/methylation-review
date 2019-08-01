rm(list=ls())
load('./data/gene.train.data.RData')
gene.mm.na <- lapply(
  gene.train.data,
  function(x){
    naCount <- apply(x[-1],2,function(x){length(which(is.na(x)))})
    return((naCount/nrow(x))*100)
  }
)
filtered.gene.train <- list()
for(i in 1:length(gene.train.data)){
  idx <- gene.mm.na[[i]] <= 10
  filtered <- gene.train.data[[i]][,-1][,idx]
  rownames(filtered) <- gene.train.data[[i]]$sampleID
  filtered.gene.train[[i]] <- filtered
}
names(filtered.gene.train) <- names(gene.train.data)

gtf <- lapply(
  filtered.gene.train,
  function(x){
   naCount <- apply(x,1,function(x){
     length(which(is.na(x)))
   })
   
   return(x[naCount == 0,])
 })
save(gtf,
     file = './data/gtf.RData')