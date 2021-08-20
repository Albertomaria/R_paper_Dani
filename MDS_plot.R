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
# mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
# #View(listAttributes(mart))
# tx4fly <- getBM(attributes = c("external_gene_name","flybase_gene_id","transcript_length","entrezgene"),mart=mart)

RNA_count <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RNA_count.csv",header = T)
RNA_iab8_count <- RNA_count %>% dplyr::select(-dplyr::contains("iab4"))

rpkm_iab8 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv", header = T)[,-c(3:11)] %>% dplyr::select(-dplyr::contains("iab4"))

RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]


for (list in c("RPKM","RPKM1","RPKM5")){
  gene_filter <- get(list)$flybase_gene_id
  
  RNA_filtered <- RNA_count %>% dplyr::filter(flybase_gene_id %in% gene_filter)
  
  RNA_group <- factor(gsub('.{3}$', '', names(RNA_count[,5:ncol(RNA_count)])))
  
  RNA_edger <- DGEList(counts = RNA_filtered[,5:ncol(RNA_filtered)],group=RNA_group,genes=RNA_filtered[,1:4])
  RNA_edger <- calcNormFactors(RNA_edger)
  
  MDSRNA <- plotMDS(RNA_edger,col=c(rep(1:3, each=3)))
  MDSRNA <- data.frame(PCA1=MDSRNA$x,PCA2=MDSRNA$y) %>% rownames_to_column("ID") %>% tidyr::separate(ID, into=c("Condition","Location","Experiment"), sep="_",remove=F)
  tot <- ggplot(MDSRNA,aes(x=PCA1,y=PCA2,color=Condition,shape=Location,label = Experiment)) +   geom_point(size=7,alpha=0.5) + 
    geom_text(color="black",size=4) + xlim(-2,2) + ylim(-1.5,1.5) +
    scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
    scale_shape_manual(values = c(16,15)) +
    labs(x="dimension1", y="dimension2", color="Genotype",shape=NULL,title=paste0("MDS plot (",nrow(RNA_filtered)," genes)")) +
    geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=14), axis.title = element_text(size=16))

  if(list == "RPKM"){filter_num = 0} else if(list == "RPKM1"){filter_num = 1} else{filter_num = 5}
  
  RNA_iab8_group <- factor(gsub('.{3}$', '', names(RNA_iab8_count[,5:ncol(RNA_iab8_count)])))
  keep_name <-data.frame("flybase_gene_id"=NULL)
  
  for (experiment in levels(RNA_iab8_group)){
    x <- (rpkm_iab8 %>% dplyr::select(dplyr::matches(paste0("flybase_gene_id|",experiment))) %>% dplyr::mutate(sum=rowSums(.[-1]>=filter_num)) %>% dplyr::filter(sum>=2))[1]
    keep_name <- rbind(keep_name,x)
  }
  
  keep_name <- unique(keep_name)
  
  RNA_iab8_filtered <- RNA_iab8_count %>% dplyr::filter(flybase_gene_id %in% unlist(keep_name))
  
  RNA_iab8_edger <- DGEList(counts = RNA_iab8_filtered[,5:ncol(RNA_iab8_filtered)],group=RNA_iab8_group,genes=RNA_iab8_filtered[,1:4])
  RNA_iab8_edger <- calcNormFactors(RNA_iab8_edger)
  
  MDSRNA_8 <- plotMDS(RNA_iab8_edger,col=c(rep(1:3, each=3)))
  MDSRNA_8 <- data.frame(PCA1=MDSRNA_8$x,PCA2=MDSRNA_8$y) %>% rownames_to_column("ID") %>% tidyr::separate(ID, into=c("Condition","Location","Experiment"), sep="_",remove=F)
  iab8 <- ggplot(MDSRNA_8,aes(x=PCA1,y=PCA2,color=Condition,shape=Location,label = Experiment)) +   geom_point(size=7,alpha=0.5) + 
    geom_text(color="black",size=4) + xlim(-2,2) + ylim(-1.5,1.5) +
    scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
    scale_shape_manual(values = c(16,15)) +
    labs(x="dimension1", y="dimension2", color="Genotype",shape=NULL,title=paste0("MDS plot (",nrow(RNA_iab8_filtered)," genes)")) +
    geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=14), axis.title = element_text(size=16))
  
  if(list == "RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/"
  } else if(list == "RPKM1"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/"
  } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/"}
  
  pdf(paste0(location,"compare_MDS_w-o_iab4.pdf"),width = 17,height = 10)
  print(cowplot::plot_grid(tot,iab8,align = "h",ncol = 2,axis = "t",labels = c("include iab4", "exclude iab4"),scale=c(1,1)))
  dev.off()
}


###PCA

for (list in c("RPKM","RPKM1","RPKM5")){
  gene_filter <- get(list)$flybase_gene_id
  
  RNA_filtered <- RNA_count %>% dplyr::filter(flybase_gene_id %in% gene_filter)
  
  RNA_group <- factor(gsub('.{3}$', '', names(RNA_count[,5:ncol(RNA_count)])))
  
  RNA_edger <- DGEList(counts = RNA_filtered[,5:ncol(RNA_filtered)],group=RNA_group,genes=RNA_filtered[,1:4])
  RNA_edger <- calcNormFactors(RNA_edger)
  
  if(list == "RPKM"){pca <- prcomp(t(RNA_edger$counts))} else{pca <- prcomp(t(RNA_edger$counts), scale. = TRUE)}
  pca_stat <- summary(pca)
  
  scores = as.data.frame(pca$x) %>% rownames_to_column("ID") %>% tidyr::separate(ID,into=c("Condition","Location","Experiment"),remove=F,sep="_")
  
  tot <- ggplot(data = scores, aes(x = PC1, y = PC2,color=Condition,shape=Location,label=Experiment)) + geom_point(size=7,alpha=0.5) +
    geom_text(color="black",size=4) +
    scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
    scale_shape_manual(values = c(16,15)) +
    labs(x=paste0("PC1 = ",round(pca_stat$importance[2,1]*100,2)), y=paste0("PC2 = ",round(pca_stat$importance[2,2]*100,2)), color="Genotype",shape=NULL,
         title=paste0("PCA plot (",nrow(RNA_filtered)," genes)")) +
    geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=14), axis.title = element_text(size=16))
  
  if(list == "RPKM"){filter_num = 0} else if(list == "RPKM1"){filter_num = 1} else{filter_num = 5}
  
  RNA_iab8_group <- factor(gsub('.{3}$', '', names(RNA_iab8_count[,5:ncol(RNA_iab8_count)])))
  keep_name <-data.frame("flybase_gene_id"=NULL)
  
  for (experiment in levels(RNA_iab8_group)){
    x <- (rpkm_iab8 %>% dplyr::select(dplyr::matches(paste0("flybase_gene_id|",experiment))) %>% dplyr::mutate(sum=rowSums(.[-1]>=filter_num)) %>% dplyr::filter(sum>=2))[1]
    keep_name <- rbind(keep_name,x)
  }
  
  keep_name <- unique(keep_name)
  
  RNA_iab8_filtered <- RNA_iab8_count %>% dplyr::filter(flybase_gene_id %in% unlist(keep_name))
  
  RNA_iab8_edger <- DGEList(counts = RNA_iab8_filtered[,5:ncol(RNA_iab8_filtered)],group=RNA_iab8_group,genes=RNA_iab8_filtered[,1:4])
  RNA_iab8_edger <- calcNormFactors(RNA_iab8_edger)
  
  if(list == "RPKM"){pca <- prcomp(t(RNA_iab8_edger$counts))} else{pca <- prcomp(t(RNA_iab8_edger$counts), scale. = TRUE)}
  pca_stat <- summary(pca)
  
  scores = as.data.frame(pca$x) %>% rownames_to_column("ID") %>% tidyr::separate(ID,into=c("Condition","Location","Experiment"),remove=F,sep="_")
  
  iab8 <- ggplot(data = scores, aes(x = PC1, y = PC2,color=Condition,shape=Location,label=Experiment)) + geom_point(size=7,alpha=0.5) +
    geom_text(color="black",size=4) +
    scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
    scale_shape_manual(values = c(16,15)) +
    labs(x=paste0("PC1 = ",round(pca_stat$importance[2,1]*100,2)), y=paste0("PC2 = ",round(pca_stat$importance[2,2]*100,2)), color="Genotype",shape=NULL,
         title=paste0("PCA plot (",nrow(RNA_iab8_filtered)," genes)")) +
    geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
    theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=14), axis.title = element_text(size=16))
  
  if(list == "RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/"
  } else if(list == "RPKM1"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/"
  } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/"}
  
  pdf(paste0(location,"compare_PCA_w-o_iab4.pdf"),width = 17,height = 10)
  print(cowplot::plot_grid(tot,iab8,align = "h",ncol = 2,axis = "t",labels = c("include iab4", "exclude iab4"),scale=c(1,1)))
  dev.off()
}

# ###PCA RPKM
# 
# #RNA_count <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RNA_count.csv",header = T)
# #RNA_iab8_count <- RNA_count %>% dplyr::select(-dplyr::contains("iab4"))
# 
# rpkm_iab8 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv", header = T)[,-c(3:11)] %>% dplyr::select(-dplyr::contains("iab4"))
# 
# RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv", header = T)[,-c(3:10)] 
# RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
# RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]
# 
# for (list in c("RPKM","RPKM1","RPKM5")){
#   gene_filter <- get(list)$flybase_gene_id
#   rpkm <- get(list)
#   
#   # RNA_filtered <- RNA_count %>% dplyr::filter(flybase_gene_id %in% gene_filter)
#   
#   RNA_group <- factor(gsub('.{3}$', '', names(rpkm[,5:ncol(rpkm)])))
#   
#   # RNA_edger <- DGEList(counts = RNA_filtered[,5:ncol(RNA_filtered)],group=RNA_group,genes=RNA_filtered[,1:4])
#   # RNA_edger <- calcNormFactors(RNA_edger)
#   
#   # if(list == "RPKM"){pca <- prcomp(t(RNA_edger$counts))} else{
#   pca <- prcomp(t(rpkm[,-c(1:4)]))
#   pca_stat <- summary(pca)
#   
#   scores = as.data.frame(pca$x) %>% rownames_to_column("ID") %>% tidyr::separate(ID,into=c("Condition","Location","Experiment"),remove=F,sep="_")
#   
#   tot <- ggplot(data = scores, aes(x = PC1, y = PC2,color=Condition,shape=Location,label=Experiment)) + geom_point(size=7,alpha=0.5) +
#     geom_text(color="black",size=4) +
#     scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
#     scale_shape_manual(values = c(16,15)) +
#     labs(x=paste0("PC1 = ",round(pca_stat$importance[2,1]*100,2)), y=paste0("PC2 = ",round(pca_stat$importance[2,2]*100,2)), color="Genotype",shape=NULL,
#          title=paste0("PCA plot (",nrow(rpkm)," genes)")) +
#     geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
#     geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
#     theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
#           axis.text = element_text(size=14), axis.title = element_text(size=16))
#   
#   if(list == "RPKM"){filter_num = 0} else if(list == "RPKM1"){filter_num = 1} else{filter_num = 5}
#   
#   RNA_iab8_group <- factor(gsub('.{3}$', '', names(rpkm_iab8[,4:ncol(rpkm_iab8)])))
#   keep_name <-data.frame("flybase_gene_id"=NULL)
#   
#   for (experiment in levels(RNA_iab8_group)){
#     x <- (rpkm_iab8 %>% dplyr::select(dplyr::matches(paste0("flybase_gene_id|",experiment))) %>% dplyr::mutate(sum=rowSums(.[-1]>=filter_num)) %>% dplyr::filter(sum>=2))[1]
#     keep_name <- rbind(keep_name,x)
#   }
#   
#   keep_name <- unique(keep_name)
#   
#   RNA_iab8_filtered <- rpkm %>% dplyr::filter(flybase_gene_id %in% unlist(keep_name))
#   
#   RNA_iab8_edger <- DGEList(counts = RNA_iab8_filtered[,5:ncol(RNA_iab8_filtered)],group=RNA_iab8_group,genes=RNA_iab8_filtered[,1:4])
#   RNA_iab8_edger <- calcNormFactors(RNA_iab8_edger)
#   
#   if(list == "RPKM"){pca <- prcomp(t(RNA_iab8_filtered[,-c(1:4)]))} else{pca <- prcomp(t(RNA_iab8_edger$counts), scale. = TRUE)}
#   pca_stat <- summary(pca)
#   
#   scores = as.data.frame(pca$x) %>% rownames_to_column("ID") %>% tidyr::separate(ID,into=c("Condition","Location","Experiment"),remove=F,sep="_")
#   
#   iab8 <- ggplot(data = scores, aes(x = PC1, y = PC2,color=Condition,shape=Location,label=Experiment)) + geom_point(size=7,alpha=0.5) +
#     geom_text(color="black",size=4) +
#     scale_color_manual(values = c("#3488bd","#d53d4f","#7fbd6f")) + 
#     scale_shape_manual(values = c(16,15)) +
#     labs(x=paste0("PC1 = ",round(pca_stat$importance[2,1]*100,2)), y=paste0("PC2 = ",round(pca_stat$importance[2,2]*100,2)), color="Genotype",shape=NULL,
#          title=paste0("PCA plot (",nrow(RNA_iab8_filtered)," genes)")) +
#     geom_hline(yintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
#     geom_vline(xintercept = 0, colour = "gray65", size = 0.5, linetype = "dotted") +
#     theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
#           axis.text = element_text(size=14), axis.title = element_text(size=16))
#   
#   if(list == "RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/"
#   } else if(list == "RPKM1"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/"
#   } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/"}
#   
#   pdf(paste0(location,"compare_PCA_w-o_iab4.pdf"),width = 17,height = 10)
#   print(cowplot::plot_grid(tot,iab8,align = "h",ncol = 2,axis = "t",labels = c("include iab4", "exclude iab4"),scale=c(1,1)))
#   dev.off()
# }
