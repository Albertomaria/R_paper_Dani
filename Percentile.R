setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")

library(tidyverse)
# library(edgeR)
library(biomaRt)
library(org.Dm.eg.db)
library(readr)
library(readxl)
library(plotly)
library(matrixStats)

########## Import Ensembl (bioMart) information
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
ecdf_fun <- function(x,perc) ecdf(x)(perc)

Predicted_Targets_score <- read_delim("E:/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")

iab4 <- unlist(unique((Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::filter(miRNA == "dme-miR-iab-4-5p"))$gene_ID))
iab8 <- unlist(unique((Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::filter(miRNA == "dme-miR-iab-8-5p"))$gene_ID))

DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T)
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T)

colnames(DEG_RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:2,39:42)])))))
colnames(DEG_1RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:2,39:42)])))))
colnames(DEG_5RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:2,39:42)])))))

for(name_DEG in c("DEG_RPKM","DEG_1RPKM","DEG_5RPKM")){
  DEG_df <- get(name_DEG)
  
  DEG_percentile <- DEG_df %>% tidyr::gather(key="Comparison","Value",-flybase_gene_id,-external_gene_name) %>% tidyr::separate(Comparison,c("Stat","Comparison"),sep="\\[") %>% 
  dplyr::mutate(Comparison = gsub("]","",Comparison)) %>% tidyr::spread(key=Stat,value=Value) #%>% dplyr::filter(PValue < 0.05,logFC > 0)
  
  DEG_iab4_percentile <- DEG_percentile %>% dplyr::filter(flybase_gene_id %in% iab4) 
  DEG_iab8_percentile <- DEG_percentile %>% dplyr::filter(flybase_gene_id %in% iab8)
  
  for(percentile in c("DEG_percentile","DEG_iab4_percentile","DEG_iab8_percentile")){
    DEG_per_percentile = get(percentile)
    for (comp in unique(DEG_per_percentile$Comparison)){
      xd <- data.frame(density(subset(DEG_per_percentile,Comparison == comp)$logFC)[c("x", "y")])
      Hth_point <- (subset(DEG_per_percentile,Comparison == comp) %>% dplyr::filter(external_gene_name == "hth"))$logFC
      Hth_y <- xd[min(which(xd$x>Hth_point)),2]
      if(is.na(Hth_y)){Hth_y = min(xd$y)}
      Hth_x <- xd[min(which(xd$x>Hth_point)),1]
      if(is.na(Hth_x)){Hth_x = min(xd$x)}
      lx = 1*max(xd$x)/10
      ly = 0.25*max(xd$y)/10
      ax = 1.5*max(xd$x)/10
      ay = 0.275*max(xd$y)/10
    
      ngene <- data.frame(table(subset(DEG_per_percentile,Comparison == comp)$logFC>Hth_point))
      
      percent <- round(ecdf_fun(subset(DEG_per_percentile,Comparison == comp)$logFC,Hth_point) * 100,2)
      
      graph <- if(length(percent) != 0 ) {
        ggplot(xd, aes(x, y)) + 
          geom_area(data = subset(xd, x >= Hth_x), fill = "#baffc9") + geom_area(data = subset(xd, x < Hth_x), fill = "#c0c5ce") +
          geom_line(size = 0.5 ) + labs(x="logFC",y="density\n",title="") +
          annotate("segment", x = Hth_x + lx, xend = Hth_x, y =Hth_y + ly , yend = Hth_y, colour = "#4363d8", size=1.5, alpha=1, arrow=arrow(angle=30,length = unit(0.1, "inches"), ends = "last", type = "closed")) +
          annotate("text", x = Hth_x + ax, y =Hth_y + ay,label = paste0("Hth (percentile =",percent,"%)"),size=5,hjust=0) +
          annotate("text", x = min(xd$x) / 2, y =max(xd$y)/2 + ay,label = paste0(ngene[which(ngene$Var1=="FALSE"),2]-1,"\ngenes"),size=5,hjust=0.5) +
          annotate("text", x = max(xd$x) / 2, y =max(xd$y)/2 + ay,label = paste0(ngene[which(ngene$Var1=="TRUE"),2],"\ngenes"),size=5,hjust=0.5) +
          theme(panel.background = element_rect(fill = "white", colour = "black"),
                axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
                axis.text.x = element_text(hjust = 1),
                plot.title = element_text(size=20, face="bold",hjust=0.5))
        } else {
          ggplot(xd, aes(x, y)) + 
            geom_area(data = subset(xd, x >= Hth_x), fill = "#baffc9") + geom_area(data = subset(xd, x < Hth_x), fill = "#c0c5ce") +
            geom_line(size = 0.5 ) + labs(x="logFC",y="density\n",title="") +
            annotate("text", x = max(xd$x)/2, y = max(xd$y),label = paste0("Hth (percentile = NC)"),size=5,hjust=0) +
            theme(panel.background = element_rect(fill = "white", colour = "black"),
                  axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
                  axis.text.x = element_text(hjust = 1),
                  plot.title = element_text(size=20, face="bold",hjust=0.5))
        }
      
      if(name_DEG == "DEG_RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Percentile/Total_FC/"
      } else if(name_DEG == "DEG_1RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/Percentile/Total_FC/"
      } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/Percentile/Total_FC/"}
      
      if(percentile=="DEG_percentile"){
        pdf(paste0(location,"Hth_total/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(graph)
        dev.off()
      } else if(percentile=="DEG_iab4_percentile") {
        pdf(paste0(location,"Hth_iab4/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(graph)
        dev.off()          
      } else {
        pdf(paste0(location,"Hth_iab8/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(graph)
        dev.off()
      }
    }
  }
}


rm(list = c("ecdf_fun","DEG_df","DEG_per_percentile","xd","Hth_point","Hth_x","Hth_y","percent","comp"))
