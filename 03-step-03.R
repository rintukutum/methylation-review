rm(list=ls())
source('./func-room.R')
load('./data/gtf.RData')

metadata <- read.csv(
  './data/meta-information.csv',
  stringsAsFactors = FALSE
)
br.all <- list()
for(i in 1:length(gtf)){
  cat(paste0(i, '. ',names(gtf)[i],'\n'))
  br.all[[i]] <- fS(
    x = gtf[[i]],
    seed = i
  )
}
names(br.all) <- names(gtf)
rM(br.all)


