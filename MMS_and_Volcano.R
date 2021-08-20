setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")
#ann <- "E:/Drosophila/ucsc_dme_r6.21/Drosophila_melanogaster.BDGP6.94.chr.gff3"

library(tidyverse)
library(edgeR)
library(biomaRt)
# library(mixOmics)
# library(RColorBrewer)
# library(goSTAG)
# library(org.Dm.eg.db)
library(readr)
library(readxl)
library(plotly)
library(eulerr)
library(matrixStats)
library(cowplot)

##### dataset
attribute_papers_iab <- read.csv("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/attribute_paper_iab.csv",header = T)

Slattery_2013 <- attribute_papers_iab[,c(1:2,4)]  %>% dplyr::filter(!is.na(Slattery_2013)) 

Slattery_2011 <- attribute_papers_iab[,c(1:2,3)] %>% dplyr::filter(!is.na(Slattery_2011)) 


DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T)
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T)

colnames(DEG_RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:2,39:42)])))))
colnames(DEG_1RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:2,39:42)])))))
colnames(DEG_5RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:2,39:42)])))))

DEG_RPKM <- merge(merge(DEG_RPKM,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)
DEG_1RPKM <- merge(merge(DEG_1RPKM,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)
DEG_5RPKM <- merge(merge(DEG_5RPKM,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)

RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv", header = T)[,-c(3:12)] 
#RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
#RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]

##### comparison table
RNA_group <- c(rep("CS_iab4",3),rep("CS_iab8",3),rep("hthBS_iab4",3),rep("hthBS_iab8",3),rep("miR_iab4",3),rep("miR_iab8",3))

design <- model.matrix(~0+RNA_group)
colnames(design) <- gsub("RNA_group","",colnames(design))

my.contrasts <- makeContrasts(`CS[iab8-iab4]`= CS_iab8 - CS_iab4,
                              `hthBS[iab8-iab4]`= hthBS_iab8 - hthBS_iab4,
                              `miR[iab8-iab4]`= miR_iab8 - miR_iab4,
                              
                              `iab4[miR-CS]`= miR_iab4 - CS_iab4,
                              `iab4[hthBS-CS]`= hthBS_iab4 - CS_iab4,
                              `iab4[miR-hthBS]`= miR_iab4 - hthBS_iab4,
                              
                              `iab8[miR-CS]`= miR_iab8 - CS_iab8,
                              `iab8[hthBS-CS]`= hthBS_iab8 - CS_iab8,
                              `iab8[miR - hthBS]`= miR_iab8 - hthBS_iab8,
                              levels=design)

##### analysis
# name_DEG <- "DEG_1RPKM"

i=1

for(name_DEG in c("DEG_RPKM","DEG_1RPKM","DEG_5RPKM")){
  DEG_df <- get(name_DEG)
  while(i<=dim(my.contrasts)[2]){
    condition <- my.contrasts[which(my.contrasts[,i]!=0),i]
    name <- paste0("[",names(which(condition>0))," - ",names(which(condition<0)),"]")
    table <- DEG_df %>% select(flybase_gene_id,external_gene_name,dplyr::contains(name),Slattery_2011,Slattery_2013) 
    colnames(table) <- c("flybase_gene_id","external_gene_name","logFC","logCPM","F","PValue","Slattery_2011","Slattery_2013")
    
    match <- paste("flybase_gene_id",names(condition)[1],names(condition)[2],sep="|")
    mean_rpkm <- RPKM %>% dplyr::select(dplyr::matches(match)) %>% dplyr::mutate(mean = log2(rowMeans(.[,2:7])+1)) %>% dplyr::select(flybase_gene_id=1 ,`log2[RPKM+1]`=8) %>%
      dplyr::mutate(exp=paste(names(which(condition>0)),names(which(condition<0)),sep=" - "))
    table = merge(table,mean_rpkm,by="flybase_gene_id",all.x=T)
    table_filtered = table %>% dplyr::filter(PValue<0.05)
    
    markD11 = table_filtered %>% dplyr::filter(is.na(Slattery_2011))
    mark11 = table_filtered %>% dplyr::filter(!is.na(Slattery_2011))
    markD13 = table_filtered %>% dplyr::filter(is.na(Slattery_2013))
    mark13 = table_filtered %>% dplyr::filter(!is.na(Slattery_2013))
    
    pt11_MC.up <- format(prop.test(x=c(table(mark11$logFC>0)[2],table(markD11$logFC>0)[2]),n=c(nrow(mark11),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt11_MC.down <- format(prop.test(x=c(table(mark11$logFC<0)[2],table(markD11$logFC<0)[2]),n=c(nrow(mark11),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt11_HC.up <- format(prop.test(x=c(table((mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence"))$logFC>0)[2],table(markD11$logFC>0)[2]),
                                n=c(nrow(mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence")),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt11_HC.down <- format(prop.test(x=c(table((mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence"))$logFC<0)[2],table(markD11$logFC<0)[2]),
                                n=c(nrow(mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence")),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    
    pt13.up <- format(prop.test(x=c(table(mark13$logFC>0)[2],table(markD13$logFC>0)[2]),n=c(nrow(mark13),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13.down <- format(prop.test(x=c(table(mark13$logFC<0)[2],table(markD13$logFC<0)[2]),n=c(nrow(mark13),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_EA.up <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W"))$logFC>0)[2],table(markD13$logFC>0)[2]),
                                n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_EA.down <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W"))$logFC<0)[2],table(markD13$logFC<0)[2]),
                                n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_EW.up <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W"))$logFC>0)[2],table(markD13$logFC>0)[2]),
                                n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_EW.down <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W"))$logFC<0)[2],table(markD13$logFC<0)[2]),
                                n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_W.up <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth W>EA"))$logFC>0)[2],table(markD13$logFC>0)[2]),
                               n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth W>EA")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    pt13_W.down <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth W>EA"))$logFC<0)[2],table(markD13$logFC<0)[2]),
                               n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth W>EA")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
    
    sc_11 <- ggplot(table,aes(y=logFC,x=`log2[RPKM+1]`,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
      geom_point(data=mark11 ,aes(y=logFC,x=`log2[RPKM+1]`, fill=Slattery_2011), size=2,shape=21) +
      scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
      scale_alpha(guide="none") +
      labs(y = paste0("logFC[",names(which(condition>0))," - ",names(which(condition<0)),"]"),x="log2[RPKM+1]") +
      geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5),
            legend.position="bottom")
    
    sc_13 <- ggplot(table,aes(y=logFC,x=`log2[RPKM+1]`,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
      geom_point(data=mark13 ,aes(y=logFC,x=`log2[RPKM+1]`, fill=Slattery_2013), size=2,shape=21) +
      scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
      scale_alpha(guide="none") +
      labs(y = paste0("logFC[",names(which(condition>0))," - ",names(which(condition<0)),"]"),x="log2[RPKM+1]") +
      geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5),
            legend.position="bottom")
    
    sc111 <- sc_11 + geom_vline(xintercept = 1, size = 0.5, color="#990000" )
    sc511 <- sc_11 + geom_vline(xintercept = 2.585, size = 0.5, color="#009900" )
    sc113 <- sc_13 + geom_vline(xintercept = 1, size = 0.5, color="#990000" )
    sc513 <- sc_13 + geom_vline(xintercept = 2.585, size = 0.5, color="#009900" )
    
    vp <- ggplot(table,aes(x= exp, y=`log2[RPKM+1]`)) + geom_violin(fill="#dedede") +
      coord_flip() +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_blank(),axis.title = element_blank(),
            plot.title = element_blank())
  
    vp1 <- vp + geom_hline(yintercept = 1, size = 0.5, color="#990000" ) +
      annotate("text",x=1.48,y=0.3,label = paste0("< 1 RPKM\n",table(table$`log2[RPKM+1]`>1)[1]), color="#990000",size=3) +
      annotate("text",x=1.48,y=1.7,label = paste0("> 1 RPKM\n",table(table$`log2[RPKM+1]`>1)[2]), color="#990000",size=3)
    vp5 <- vp + geom_hline(yintercept = 2.585, size = 0.5, color="#009900" ) +
      annotate("text",x=1.48,y=1.7,label = paste0("< 5 RPKM\n",table(table$`log2[RPKM+1]`>2.585)[1]), color="#009900",size=3) +
      annotate("text",x=1.48,y=3.3,label = paste0("> 5 RPKM\n",table(table$`log2[RPKM+1]`>2.585)[2]), color="#009900",size=3)
    
    Volp11 <- ggplot(table,aes(y=-log10(PValue),x=logFC,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
      geom_point(data=mark11,aes(y=-log10(PValue),x=logFC, fill=Slattery_2011), size=3,shape=21) +
      scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
      scale_alpha(guide="none") +
      labs(y = "log10[pValue]",title = paste0(names(which(condition>0))," - ",names(which(condition<0))),x="logFC") +
      annotate("text",x=-7,y=17,label=paste0("Down-regulate pValue = ",pt11_MC.down),color="#3cd44b",size=3,hjust=0) +
      annotate("text",x=7,y=17,label=paste0("Up-regulate pValue = ",pt11_MC.up),color="#3cd44b",size=3,hjust=1) +
      annotate("text",x=0,y=17,label="Medioum confidence Hth target",color="#3cd44b",size=4,hjust=0.5) +
      annotate("text",x=-7,y=16.5,label=paste0("Down-regulate pValue = ",pt11_HC.down),color="#e6194b",size=3,hjust=0) +
      annotate("text",x=7,y=16.5,label=paste0("Up-regulate pValue = ",pt11_HC.up),color="#e6194b",size=3,hjust=1) +
      annotate("text",x=0,y=16.5,label="High confidence Hth target",color="#e6194b",size=4,hjust=0.5) +
      geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
      geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5),
            legend.position="bottom")
    
    Volp13 <- ggplot(table,aes(y=-log10(PValue),x=logFC,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
      geom_point(data=mark13,aes(y=-log10(PValue),x=logFC, fill=Slattery_2013), size=3,shape=21) +
      scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
      scale_alpha(guide="none") +
      labs(y = "log10[pValue]",title = paste0(names(which(condition>0))," - ",names(which(condition<0))),x="logFC") +
      annotate("text",x=-7,y=17,label=paste0("Down-regulate pValue = ",pt13.down),size=3,hjust=0) +
      annotate("text",x=7,y=17,label=paste0("Up-regulate pValue = ",pt13.up),size=3,hjust=1) +
      annotate("text",x=0,y=17,label="Hth target Total [Slattery 2013]",size=4,hjust=0.5) +
      
      annotate("text",x=-7,y=16.5,label=paste0("Down-regulate pValue = ",pt13_EW.down),color="#e6194b",size=3,hjust=0) +
      annotate("text",x=7,y=16.5,label=paste0("Up-regulate pValue = ",pt13_EW.up),color="#e6194b",size=3,hjust=1) +
      annotate("text",x=0,y=16.5,label="Hth target EA=W [Slattery 2013]",color="#e6194b",size=4,hjust=0.5) +
      
      annotate("text",x=-7,y=16,label=paste0("Down-regulate pValue = ",pt13_EA.down),color="#3cd44b",size=3,hjust=0) +
      annotate("text",x=7,y=16,label=paste0("Up-regulate pValue = ",pt13_EA.up),color="#3cd44b",size=3,hjust=1) +
      annotate("text",x=0,y=16,label="Hth target EA>W [Slattery 2013]",color="#3cd44b",size=4,hjust=0.5) +
      
      annotate("text",x=-7,y=15.5,label=paste0("Down-regulate pValue = ",pt13_W.down),color="#ffe119",size=3,hjust=0) +
      annotate("text",x=7,y=15.5,label=paste0("Up-regulate pValue = ",pt13_W.up),color="#ffe119",size=3,hjust=1) +
      annotate("text",x=0,y=15.5,label="Hth target W>EA [Slattery 2013]",color="#ffe119",size=4,hjust=0.5) + 
      
      geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
      geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5),
            legend.position="bottom")
  
    ps111 <- ggplotly(sc111,tooltip = "label")
    ps113 <- ggplotly(sc113,tooltip = "label")
    ps511 <- ggplotly(sc511,tooltip = "label")
    ps513 <- ggplotly(sc513,tooltip = "label")
    pv11 <- ggplotly(Volp11,tooltip = "label")
    pv13 <- ggplotly(Volp13,tooltip = "label")
    
    
    if(name_DEG == "DEG_RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/"
    } else if(name_DEG == "DEG_1RPKM"){location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/"
    } else{location = "L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/"}
    
    
    pdf(paste0(location,"MMD/1RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(cowplot::plot_grid(vp1,sc111,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
    dev.off()
    htmlwidgets::saveWidget(ps111, paste0(location,"MMD/1RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))

    pdf(paste0(location,"MMD/1RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(cowplot::plot_grid(vp1,sc113,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
    dev.off()
    htmlwidgets::saveWidget(ps113, paste0(location,"MMD/1RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))

    pdf(paste0(location,"MMD/5RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(cowplot::plot_grid(vp5,sc511,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
    dev.off()
    htmlwidgets::saveWidget(ps511, paste0(location,"MMD/5RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))

    pdf(paste0(location,"MMD/5RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(cowplot::plot_grid(vp5,sc513,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
    dev.off()
    htmlwidgets::saveWidget(ps513, paste0(location,"MMD/5RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))

    pdf(paste0(location,"Volcano/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(Volp11)
    dev.off()
    htmlwidgets::saveWidget(pv11, paste0(location,"Volcano/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))
    
    pdf(paste0(location,"Volcano/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(Volp13)
    dev.off()
    htmlwidgets::saveWidget(pv13, paste0(location,"Volcano/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))
    
    i=i+1
  }
}

rm(list = c("i","Slattery_2011","Slattery_2013","DEG_df_miR","condition","name","table","DEG_df","match","mean_rpkm","mark11","mark13","markD11","markD13",ls()[grep("sc",ls())],ls()[grep("vp",ls())],
            ls()[grep("Volp",ls())],ls()[grep("ps",ls())],ls()[grep("pv",ls())],ls()[grep("pt",ls())]))
