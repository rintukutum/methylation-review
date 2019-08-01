rm(list=ls())
source('./func-room.R')
load('./data/rM.data.RData')
idxMgenes <- sapply(
  rM.data,
  function(x){!any(is.na(x$imp.features))})
mGenes <- names(rM.data)[idxMgenes]
mData <- rM.data[idxMgenes]
rMod(mData)
rm(list=ls())
load('./data/rMOD.RData')
rmod <- rMOD$mod
mdperf <- list()
for(i in 1:length(rmod)){
  mat <- data.frame(rmod[[i]]$confusion)
  mat$gene <- names(rmod)[i]
  mat$original <- rownames(mat)
  mdperf[[i]] <- mat
}
mdperf.df <- plyr::ldply(mdperf)
tmp <- mdperf.df[,colnames(mdperf.df)[c(4:5,1:3)]]
tmp$gene[(seq(1,nrow(tmp),by = 2)-1)[-1]] <- '-'
write.csv(
  tmp,
  './data/Table-01.csv',
  row.names = FALSE
)
###---------
rm(list=ls())
load('./data/rM.data.RData')
source('./func-room.R')
rModx <- runMOD(rM.data)

png('./figures/MDS-plot-final-300dpi.png',
    width = 1100,
    height = 1100,
    res = 300)
vMapP1(rModx)
dev.off()
png('./figures/Importance-300dpi.png',
    width = 1100,
    height = 1100,
    res = 300)
vMapP2(rModx)
dev.off()