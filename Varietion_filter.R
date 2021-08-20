setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")

library(tidyverse)
library(edgeR)
library(biomaRt)
#library(goSTAG)
library(org.Dm.eg.db)
library(readr)
library(readxl)
library(plotly)
library(matrixStats)


########## Import Ensembl (bioMart) information
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
RNA_count <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RNA_count.csv",header = T)

rpkm <- cbind(RNA_count[1:4],rpkm(RNA_count[,5:22],gene.length=RNA_count$transcript_length))
RNA_rpkm <- rpkm %>% tidyr::gather(key="Experiment",value = "rpkm",-flybase_gene_id,-external_gene_name,-entrezgene,-transcript_length) %>% dplyr::mutate(group = gsub('.{3}$', '', Experiment))

# rpkm_annotate <- unique(merge(rpkm,attribute_papers_iab[,-2],by="flybase_gene_id",all.x=T) %>% dplyr::select(1:2,4,23:30,dplyr::everything()))
# rpkm_annotate[is.na(rpkm_annotate)] <- ""

RNA_group <- factor(RNA_rpkm$group)

rpkm_mean <- rpkm[,1:2]

for (exp in levels(RNA_group)){
  x <- rpkm[,grep(exp,names(rpkm))] %>% dplyr::mutate(!!paste0(exp,"_mean"):=rowMeans(.[,1:3]),!!paste0(exp,"_SD"):=rowSds(as.matrix(.[,1:3]))) %>% 
    dplyr::mutate(!!paste0(exp,"_CV"):=.[,5]/.[,4])
  rpkm_mean = cbind(rpkm_mean,x[,4:6])
}

rpkm_mean[is.na(rpkm_mean)] <- 0

# ggplot(rpkm_mean,aes(x=log2(CS_iab4_CV),y=log2(CS_iab8_CV))) + geom_point() + geom_abline(slope = 1, intercept = 0)



RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]

ggplot(RPKM1,aes(x=log2(miR_iab8_SD),y=log2(CS_iab8_SD))) + geom_point() + geom_abline(slope = 1, intercept = 0) + xlim(-10,15) + ylim(-10,15)
ggplot(RPKM5,aes(x=miR_iab8_SD/miR_iab8_mean,y=CS_iab8_SD/CS_iab8_mean)) + geom_point() + geom_abline(slope = 1, intercept = 0) + xlim(0,2) + ylim(0,2)
