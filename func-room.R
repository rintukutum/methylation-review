library(randomForest)
library(Boruta)
library(plyr)
library(reshape)
library(ggplot2)
load.data <- function(){
  metadata <- read.csv(
    './data/meta-information.csv',
    stringsAsFactors = FALSE
  )
  raw.data <- read.csv(
    './data/raw-filtered-data.csv',
    stringsAsFactors = FALSE
  )
  return(
    list(
      metadata = metadata,
      raw.data = raw.data
      )
  )
}
getGeneData <- function(){
  gene.level.data <- dlply(
    raw.data,
    'Gene_name',
    function(x){
      chr.pos <- x[,2:4]
      methyl <- x[,-c(1:4)]
      methyl.l <- apply(methyl,1,na.omit)
      return(list(
        chr.pos = chr.pos,
        methyl = methyl.l
      ))
    }
  )
  return(gene.level.data)
}
selectPosition <- function(x){
  sampleID <- names(x)
  nCount <- length(grep('N',sampleID))
  dCount <- length(grep('C',sampleID))
  if(nCount >= 2 & dCount >= 2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

filterGenePosition <- function(x){
  xx.pos <- x$chr.pos
  xx.methyl <- x$methyl
  idx <- sapply(xx.methyl,selectPosition)
  out <- list(chr.pos = xx.pos[idx,],
              methyl = xx.methyl[idx]
  )
  return(out)
}
extractTrainGene <- function(x){
  feature.name <- apply(
    x$chr.pos,
    1,
    function(x){
      x <- as.character(x)
      return(
        paste0(
          x[2],
          '_',
          x[3],
          '_',
          x[1]
        )
      )
    }
  )
  #  values
  methyl.mat <- lapply(
    x$methyl,
    function(x){
      val <- data.frame(x)
      val$sampleID <- rownames(val)
      colnames(val)[1] <- 'methyl.perc'
      rownames(val) <- 1:nrow(val)
      return(val)
    }
  )
  names(methyl.mat) <- feature.name
  methyl.df <- ldply(methyl.mat)
  colnames(methyl.df)[1] <- 'feature'
  methyl.train <- cast(
    methyl.df,
    sampleID~feature,
    value='methyl.perc'
  )
  return(methyl.train)
}
fS <- function(x,seed){
  y <- sapply(rownames(x),function(x){strsplit(x,split='')[[1]][1]})
  org.var <- colnames(x)
  colnames(x) <- paste0('F',1:ncol(x))
  names(org.var) <- colnames(x)
  tr <- data.frame(class=y,x)
  set.seed(seed)
  br <- Boruta(class~., data = tr)
  return(list(org.var = org.var,br=br,tr=tr))
}
rM <- function(br.all){
  rfD <- list()
  for(i in 1:length(br.all)){
    idx <- as.character(br.all[[i]]$br$finalDecision) %in% c('Tentative','Confirmed')
    if(any(idx)){
      imp.features <- names(br.all[[i]]$br$finalDecision)[idx]
      tr <- data.frame(br.all[[i]]$tr[,c('class',imp.features)])
      imp.features <- br.all[[i]]$org.var[imp.features]
      colnames(tr) <- c('class',imp.features)
    }else{
      imp.features <- NA
      tr <- br.all[[i]]$tr
      colnames(tr) <- c('class',br.all[[i]]$org.var)
    }
    rfD[[i]] <- list(
      imp.features = imp.features,
      tr = tr
    )
  }
  names(rfD) <- names(br.all)
  rM.data <- rfD
  save(rM.data,file = './data/rM.data.RData')
  return(cat('Done\n'))
}
rMod <- function(mData){
  nFpG <- sapply(mData,function(x){ncol(x$tr)})
  mod <- list()
  for(i in 1:length(mData)){
    tr <- mData[[i]]$tr
    colnames(tr) <- gsub(' ','',colnames(tr))
    cat(paste0(i,'. ',names(mData)[i],'\n'))
    set.seed(i+34)
    mod[[i]] <- randomForest(
      class~., data=tr, ntree=20000,
      proximity=TRUE,
      importance=TRUE)
  }
  names(mod) <- names(mData)
  rMOD <- list(
    mod=mod,
    nFpG=nFpG
  )
  save(rMOD,
       file = './data/rMOD.RData')
  return(cat('Done\n'))
}
runMOD <- function(rM.data){
  sampleIDs <- lapply(
    rM.data,
    function(x){
      rownames(x$tr)
    }
  )
  all.s <- unique(unlist(sampleIDs))
  for(i in 1:length(sampleIDs)){
    if(i == 1){
      common.s <- intersect(all.s,sampleIDs[[i]])  
    }else{
      common.s <- intersect(common.s,sampleIDs[[i]])
    }
  }
  abcg1 <- sampleIDs$ABCG1
  abcg1.g <- list()
  for(i in 1:length(sampleIDs)){
    if(i == 1){
      abcg1.c <- intersect(abcg1,sampleIDs[[i]])  
      abc.tmp <- abcg1.c
    }else{
      abcg1.c <- intersect(abcg1.c,sampleIDs[[i]])
      abc.tmp <- intersect(abcg1.c,sampleIDs[[i]])
    }
    abcg1.g[[i]] <- abc.tmp
  }
  names(abcg1.g) <- names(sampleIDs)
  rdat <- lapply(
    rM.data,
    function(x){
      tr <- data.frame(x$tr[common.s,-1])
      colnames(tr) <- colnames(x$tr)[-1]
      tr
    }
  )
  for(i in 1:length(rdat)){
    tmp <- rdat[[i]]
    colnames(tmp) <- paste0(
      names(rdat)[i],
      '.',colnames(tmp))
    if(i == 1){
      rdat.mini <- tmp
    }else{
      rdat.mini <- cbind(
        rdat.mini, tmp
      )
    }
  }
  class <- sapply(
    rownames(rdat.mini),
    function(x){strsplit(x,split='')[[1]][1]})
  set.seed(562)
  bt <- Boruta(
    x = rdat.mini,
    y = as.factor(class)
  )
  imp.vars <- names(bt$finalDecision[bt$finalDecision %in% c('Tentative','Confirmed')])
  set.seed(561)
  rMod <- randomForest(
    x=rdat.mini[,imp.vars],
    y = as.factor(class),
    ntree=30000,
    proximity = TRUE
  )
  return(rMod)
}
vMapP1 <- function(rModx){
  md <- MDSplot(rModx,rModx$y)
  mdp <- data.frame(md$points,class=rModx$y)
  p1 <- ggplot(mdp,aes(x=Dim.1,y=Dim.2)) +
    geom_point(aes(fill=class),size=3.5,shape=21,alpha=0.65)+
    scale_fill_manual(
      name = 'Class',
      values=c(
        'C' = '#aa00d4ff',
        'N' = '#2c89a0ff'
      ),
      label = c('Case','Control')
    ) +
    ggtitle(paste0('Accuracy=',100-rModx$err.rate[30000]*100,'%'))+
    theme(legend.position = c(0.15,0.8415))
  xx <- rownames(rModx$importance)[order(rModx$importance[,1])]
  df_ <- data.frame(
    varName = rownames(rModx$importance),
    Importance = rModx$importance[,1]
  )
  p2 <- ggplot(df_,aes(x=varName,y=Importance)) +
    geom_bar(stat='identity') +
    coord_flip() +
    xlab('Gene(methylation position)') +
    ggtitle(paste0('Accuracy=',100-rModx$err.rate[30000]*100,'%')) +
    scale_x_discrete(limit=xx)
  return(p1)
}
vMapP2 <- function(rModx){
  xx <- rownames(rModx$importance)[order(rModx$importance[,1])]
  df_ <- data.frame(
    varName = rownames(rModx$importance),
    Importance = rModx$importance[,1]
  )
  p2 <- ggplot(df_,aes(x=varName,y=Importance)) +
    geom_bar(stat='identity') +
    coord_flip() +
    xlab('Gene(methylation position)') +
    ggtitle(paste0('Accuracy=',100-rModx$err.rate[30000]*100,'%')) +
    scale_x_discrete(limit=xx)
  return(p2)
}
extractCoV <- function(x,gene.name){
  idx.normal <- sapply(
    rownames(x),
    function(x){strsplit(x,split = '')[[1]][1]}) == 'N'
  coV <- function(x){
    (sd(x)/mean(x))
  }
  normal.x <- apply(x[idx.normal,],2,coV)
  case.x <- apply(x[!idx.normal,],2,coV)
  return(
    data.frame(
      normal = normal.x,
      case=case.x
    )
  )
}
