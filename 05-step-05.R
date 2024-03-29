rm(list=ls())
library(plyr)
library(ggplot2)
library(reshape)
source('./func-room.R')
load('./data/gtf.RData')
t.out.gene <- list()
for(i in 1:length(gtf)){
  
  cat(paste0(i,'. ',names(gtf)[i],'\n'))
  gene.x <- gtf[[i]]
  Class <- sapply(
    rownames(gene.x),
    function(x){
      strsplit(x,split='')[[1]][1]
    }
  )
  x.out <- apply(gene.x,
                 2,
                 function(x){
                   x <- as.numeric(x)
                   
                   df_ <- data.frame(x = x,y=Class)
                   tout <- t.test(x~y,data = df_)
                   x.out <- data.frame(
                     t.stat = tout$statistic,
                     p.value = tout$p.value
                   )
                 })
  x.df <- ldply(
    x.out
  )
  rownames(x.df) <- colnames(gene.x)
  t.out.gene[[i]] <- x.df
}
names(t.out.gene) <- names(gtf)
t.out.gene <- lapply(
  t.out.gene,
  function(x){
    colnames(x)[1] <- 'position'
    rownames(x) <- 1:nrow(x)
    x
  }
)
gene.stats <- ldply(
  t.out.gene
)
colnames(gene.stats)[1] <- 'gene'

sig.val <- c(0.05,0.01,0.001,0.0001)
out.l <- list()
perc.sig <- c()
for(i in 1:length(sig.val)){
  idx.sig <- gene.stats$p.value <= sig.val[i]
  out <- table(idx.sig)
  out.l[[i]] <- out
  perc.sig[i] <- round(
    (out['TRUE']/sum(out))*100
    ,
    digits = 4)
  
}
names(perc.sig) <- sig.val
sig.cri <- data.frame(
  perc.sig,
  p.val.threshold = names(perc.sig),
  count = sapply(out.l,function(x){x['TRUE']})
)
rownames(sig.cri) <- 1:nrow(sig.cri)



coV <- list()
for(i in 1:length(gtf)){
  coV[[i]] <- extractCoV(
    gtf[[i]],
    names(gtf)[i]
  )
}

names(coV) <- names(gtf)
coV.df <- ldply(coV)
idx.nan.n <- is.na(coV.df$normal)
idx.nan.c <- is.na(coV.df$case)
coV.df.f <- coV.df[!(idx.nan.c | idx.nan.n),]
coV.df.f.m <- melt(coV.df.f)
coV.order <- unlist(dlply(
  coV.df.f.m,
  '.id',
  function(x){
    median(x[x$variable == 'normal','value'])
  }
))

coV.df.f.m$.id <- factor(
  coV.df.f.m$.id,
  levels = rev(names(coV.order[order(coV.order)]))
)
png('./figures/variability-300dpi.png',
    width = 2600,
    height = 3000,
    res = 300)
ggplot(coV.df.f.m,aes(variable,value)) +
  geom_boxplot() +
  facet_wrap(facets = '.id',scales = 'free_y') +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  ylab('Coefficient of variation(SD/mean)') +
  xlab('Class')
dev.off()
unlist(dlply(coV.df.f.m,
      '.id',
      function(x){
        tout <- wilcox.test(value~variable,data = x)
        return(tout$p.value)
      }))