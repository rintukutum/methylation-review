rm(list=ls())
dir.create('./figures',showWarnings = FALSE)
dir.create('./data',showWarnings = FALSE)
source('./func-room.R')
xx <- load.data()
metadata <- xx$metadata
raw.data <- xx$raw.data
rm(xx)
gene.level.data <- getGeneData()
filter.gene.level <- lapply(
  gene.level.data,
  filterGenePosition
)
gene.train.data <- lapply(
  filter.gene.level,
  extractTrainGene
)
save(gene.train.data,
     file = './data/gene.train.data.RData'
)
