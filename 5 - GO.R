setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")
#ann <- "E:/Drosophila/ucsc_dme_r6.21/Drosophila_melanogaster.BDGP6.94.chr.gff3"

library(GO.db)
library(GOstats)
library(org.Dm.eg.db)
library(tidyverse)
library(biomaRt)
library(readr)
library(readxl)
library(plotly)
library(eulerr)
library(matrixStats)
library(cowplot)


########## Import Ensembl (bioMart) information
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 

mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
#View(listAttributes(mart))
tx4fly <- getBM(attributes = c("external_gene_name","flybase_gene_id","transcript_length","entrezgene"),mart=mart)

DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T) %>% 
  merge(unique(tx4fly[,c(2,4)]),.,by="flybase_gene_id",all.y=T) %>% dplyr::filter(!is.na(entrezgene))
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T) %>%
  merge(unique(tx4fly[,c(2,4)]),.,by="flybase_gene_id",all.y=T) %>% dplyr::filter(!is.na(entrezgene))
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T) %>%
  merge(unique(tx4fly[,c(2,4)]),.,by="flybase_gene_id",all.y=T) %>% dplyr::filter(!is.na(entrezgene))

colnames(DEG_RPKM) <- c("flybase_gene_id","Entrez_ID","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:3)])))))
colnames(DEG_1RPKM) <- c("flybase_gene_id","Entrez_ID","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:3)])))))
colnames(DEG_5RPKM) <- c("flybase_gene_id","Entrez_ID","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:3)])))))

for(name_DEG in c("DEG_RPKM","DEG_1RPKM","DEG_5RPKM")){
  DEG_df <- get(name_DEG)
  
  if(name_DEG == "DEG_RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/"
  } else if(name_DEG == "DEG_1RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/"
  } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/"}
  
  UP_miR.CS <- unique(DEG_df %>% dplyr::filter(`logFC[miR_iab8 - CS_iab8]`> 0.5849625 & `PValue[miR_iab8 - CS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                        dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  UP_miR.CS <- structure(as.integer(UP_miR.CS$Entrez_ID), names = as.character(UP_miR.CS$gene_ID))
  DOWN_miR.CS <- unique(DEG_df %>% dplyr::filter(`logFC[miR_iab8 - CS_iab8]`< -0.5849625 & `PValue[miR_iab8 - CS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                          dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  DOWN_miR.CS <- structure(as.integer(DOWN_miR.CS$Entrez_ID), names = as.character(DOWN_miR.CS$gene_ID))
  
  UP_hth.CS <- unique(DEG_df %>% dplyr::filter(`logFC[hthBS_iab8 - CS_iab8]`> 0.5849625 & `PValue[hthBS_iab8 - CS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                        dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  UP_hth.CS <- structure(as.integer(UP_hth.CS$Entrez_ID), names = as.character(UP_hth.CS$gene_ID))
  
  DOWN_hth.CS <- unique(DEG_df %>% dplyr::filter(`logFC[hthBS_iab8 - CS_iab8]`< -0.5849625 & `PValue[hthBS_iab8 - CS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                          dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  DOWN_hth.CS <- structure(as.integer(DOWN_hth.CS$Entrez_ID), names = as.character(DOWN_hth.CS$gene_ID))
  
  UP_miR.hth <- unique(DEG_df %>% dplyr::filter(`logFC[miR_iab8 - hthBS_iab8]`> 0.5849625 & `PValue[miR_iab8 - hthBS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                         dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  UP_miR.hth <- structure(as.integer(UP_miR.hth$Entrez_ID), names = as.character(UP_miR.hth$gene_ID))
  
  DOWN_miR.hth <- unique(DEG_df %>% dplyr::filter(`logFC[miR_iab8 - hthBS_iab8]`< -0.5849625 & `PValue[miR_iab8 - hthBS_iab8]`<0.05) %>% dplyr::select(gene_ID=external_gene_name,Entrez_ID) %>% 
                           dplyr::group_by(gene_ID) %>% dplyr::summarise(Entrez_ID = dplyr::first(Entrez_ID))) 
  DOWN_miR.hth <- structure(as.integer(DOWN_miR.hth$Entrez_ID), names = as.character(DOWN_miR.hth$gene_ID))
  
  UP_miR.CS_vs_hth.CS_iab8 <- unique(read.csv(paste0(location,"DEG_comparison/excel/iab8[miR - CS]_VS_iab8[hthBS - CS].csv"),header = T) %>%
                                       dplyr::filter(UpUp == 1) %>% dplyr::select(Gene_name) %>% merge(.,unique(tx4fly[,c(1,4)]),by.x="Gene_name",by.y="external_gene_name",all.x=T) %>% 
                                       dplyr::group_by(Gene_name) %>% dplyr::summarise(Entrez_ID = dplyr::first(entrezgene)))
  UP_miR.CS_vs_hth.CS_iab8 <- structure(as.integer(UP_miR.CS_vs_hth.CS_iab8$Entrez_ID), names = as.character(UP_miR.CS_vs_hth.CS_iab8$Gene_name))
  
  DOWN_miR.CS_vs_hth.CS_iab8 <- unique(read.csv(paste0(location,"DEG_comparison/excel/iab8[miR - CS]_VS_iab8[hthBS - CS].csv"),header = T) %>%
                                         dplyr::filter(DownDown == 1) %>% dplyr::select(Gene_name) %>% merge(.,unique(tx4fly[,c(1,4)]),by.x="Gene_name",by.y="external_gene_name",all.x=T) %>% 
                                         dplyr::group_by(Gene_name) %>% dplyr::summarise(Entrez_ID = dplyr::first(entrezgene)))
  DOWN_miR.CS_vs_hth.CS_iab8 <- structure(as.integer(DOWN_miR.CS_vs_hth.CS_iab8$Entrez_ID), names = as.character(DOWN_miR.CS_vs_hth.CS_iab8$Gene_name))
  
  # entrezUniverse <- unlist(unique(DEG_df %>% dplyr::select(Entrez_ID,starts_with("PValue")) %>% 
  #                                   dplyr::mutate(min_PValue = rowMins(as.matrix(.[,2:10]))) %>% dplyr::filter(min_PValue  < 0.05) %>% dplyr::select(Entrez_ID)))
  entrezUniverse <- structure(as.integer(DEG_df$Entrez_ID), names = as.character(DEG_df$external_gene_name))
    
  hgCutoff <- 0.000000001
  
  for (condition in c("UP_miR.CS","DOWN_miR.CS","UP_hth.CS","DOWN_hth.CS","UP_miR.hth","DOWN_miR.hth","UP_miR.CS_vs_hth.CS_iab8","DOWN_miR.CS_vs_hth.CS_iab8")){
    select <- get(condition)
    
    paramsBP <- new("GOHyperGParams", geneIds=select, universeGeneIds=entrezUniverse, annotation="org.Dm.eg.db", ontology="BP", pvalueCutoff=hgCutoff,conditional=FALSE, testDirection="over")
    paramsMF <- new("GOHyperGParams", geneIds=select, universeGeneIds=entrezUniverse, annotation="org.Dm.eg.db", ontology="MF", pvalueCutoff=hgCutoff,conditional=FALSE, testDirection="over")
    paramsCC <- new("GOHyperGParams", geneIds=select, universeGeneIds=entrezUniverse, annotation="org.Dm.eg.db", ontology="CC", pvalueCutoff=hgCutoff,conditional=FALSE, testDirection="over")
    
    BPOver <- hyperGTest(paramsBP)
    MFOver <- hyperGTest(paramsMF)
    CCOver <- hyperGTest(paramsCC)
    
    # BP <- summary(BPOver)
    # MF <- summary(MFOver)
    # CC <- summary(CCOver)
    
    report <- rbind( summary(BPOver, pvalue=0.01) %>% dplyr::mutate(GO="BP",genes_name=NA) %>% dplyr::rename("GO_ID"=1),
                    summary(MFOver, pvalue=0.01) %>% dplyr::mutate(GO="MF",genes_name=NA) %>% dplyr::rename("GO_ID"=1),
                    summary(CCOver, pvalue=0.01) %>% dplyr::mutate(GO="CC",genes_name=NA) %>% dplyr::rename("GO_ID"=1)) %>%
      dplyr::arrange(Pvalue) %>% dplyr::mutate(Term = gsub(",","-",Term))

    GO_ext <- as.list(org.Dm.egGO2ALLEGS)
    
    for (n in 1:nrow(report)){
      if (sum(names(GO_ext) == report[n,1]) > 0 ){
        report[n,9] <- paste(names(select[select %in% GO_ext[[report[n,1]]]]),collapse=' / ')
      }
    }
                      
    write.csv(report, file=paste0(location,"GO/",condition,".csv"),row.names = F,quote = F)
  }
}

# list.files("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/excel/")

UP_miR.CS_vs_hth.CS_iab8 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/excel/iab8[miR - CS]_VS_iab8[hthBS - CS].csv",header = T) %>%
  dplyr::filter(UpUp == 1) %>% dplyr::select(Gene_name) %>% merge(.,unique(tx4fly[,c(1,4)]),by.x="Gene_name",by.y="external_gene_name",all.x=T) %>% dplyr::group_by(Gene_name) %>%
  dplyr::summarise(Entrez_ID = dplyr::first(entrezgene))

DOWN_miR.CS_vs_hth.CS_iab8 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/excel/iab8[miR - CS]_VS_iab8[hthBS - CS].csv",header = T) %>%
  dplyr::filter(DownDown == 1) %>% dplyr::select(Gene_name) %>% merge(.,unique(tx4fly[,c(1,4)]),by.x="Gene_name",by.y="external_gene_name",all.x=T) %>% dplyr::group_by(Gene_name) %>%
  dplyr::summarise(Entrez_ID = dplyr::first(entrezgene))

