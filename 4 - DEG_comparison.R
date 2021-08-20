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

########## Import Ensembl (bioMart) information
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
# View(listAttributes(mart))
tx5fly <- getBM(attributes = c("external_gene_name","flybase_gene_id","transcript_length","entrezgene","name_1006","namespace_1003"),mart=mart) %>% 
  dplyr::filter(namespace_1003 == "biological_process") %>% dplyr::group_by(flybase_gene_id) %>% 
  dplyr::summarise(external_gene_name=dplyr::first(external_gene_name), transcript_length=max(transcript_length), 
                   entrezgene=dplyr::first(entrezgene), N_BP_term = n(), BP_term = paste(name_1006,collapse = " / "))

rpkm <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv")[,c(1:2,13:30)]
annotation <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv")[,c(1,4:5)] %>% 
  dplyr::mutate(hth=ifelse(Slattery_2011!="" & Slattery_2013 =="","Slattery_2011",
                           ifelse(Slattery_2011=="" & Slattery_2013 !="","Slattery_2013",
                                  ifelse(Slattery_2011!="" & Slattery_2013 !="","both",""))))

DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T)
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T)

# DEG_df_miR <- merge(merge(DEG_df,iab4[,-2],by = "flybase_gene_id",all.x=T),iab8[,-2],by = "flybase_gene_id",all.x=T)
colnames(DEG_RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:2,39:42)])))))
colnames(DEG_1RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:2,39:42)])))))
colnames(DEG_5RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:2,39:42)])))))

########## compare_scatter

significance.threshold <- 0.05
FC.threshold <- log2(1.5)

condition = c(gsub("logFC","", names(DEG_RPKM)[grep("logFC",names(DEG_RPKM))]))

for(name_DEG in c("DEG_RPKM","DEG_1RPKM","DEG_5RPKM")){
  DEG_df <- get(name_DEG)
  
  for ( j in (1:length(condition))){
    k <- j + 1
    while (k <= length(condition)){
      a <- c(condition[j],condition[k])
      
      xlabel=a[1]
      ylabel=a[2]
      
      part1 <- names(which((table(gsub("\\[|\\]","",unlist(strsplit(strsplit(xlabel," - ")[[1]],"_")))))>1))
      xrename <- paste0(part1,gsub("_","",gsub(part1,"",xlabel)))
      part2 <- names(which((table(gsub("\\[|\\]","",unlist(strsplit(strsplit(ylabel," - ")[[1]],"_")))))>1))
      yrename <- paste0(part2,gsub("_","",gsub(part2,"",ylabel)))
      
      xdf <- DEG_df[,c(1,2,grep(xlabel,names(DEG_df),fixed=T))] %>% dplyr::rename("logFC"=3, "logCPM"=4, "F"=5, "PValue"=6)
      ydf <- DEG_df[,c(1,2,grep(ylabel,names(DEG_df),fixed=T))] %>% dplyr::rename("logFC"=3, "logCPM"=4, "F"=5, "PValue"=6)
      
      sig.matrix <- data.frame(cbind(xdf$logFC, 2^xdf$logFC, xdf$PValue,p.adjust(xdf$PValue, method="BH",n=length(xdf$PValue)),
                                     ydf$logFC, 2^ydf$logFC, ydf$PValue,p.adjust(ydf$PValue, method="BH",n=length(ydf$PValue))))
      
      colnames(sig.matrix) <- c("XlogFC","originalFCvalX","XpValue","XadjPv","YlogFC","originalFCvalY","YpValue","YadjPv")
      sig.matrix["Gene_name"] <- xdf$external_gene_name
      sig.matrix["Gene_ID"] <- ydf$flybase_gene_id
      
      
      resultMatrix <- sig.matrix %>% dplyr::filter(XpValue<significance.threshold|YpValue<significance.threshold) %>%
        merge(.,annotation[,c(1,4)],by.x="Gene_ID",by.y="flybase_gene_id") %>%
        dplyr::mutate(level1Up = ifelse(XlogFC > FC.threshold,1,""),
                      level1Down = ifelse(XlogFC < -FC.threshold,1,""),
                      level2Up = ifelse(YlogFC > FC.threshold,1,""),
                      level2Down = ifelse(YlogFC < -FC.threshold,1,"")) %>%  
        dplyr::mutate( UpUp = ifelse(level1Up == 1 & level2Up == 1,1,""),
                       DownUp = ifelse(level1Down == 1 & level2Up == 1,1,""),
                       UpDown = ifelse(level1Up == 1 & level2Down == 1,1,""),
                       DownDown = ifelse(level1Down == 1 & level2Down == 1,1,"")) %>%
        dplyr::mutate(level1Up = ifelse(UpUp==1|UpDown==1,"",level1Up),
                      level1Down = ifelse(DownUp==1|DownDown==1,"",level1Down),
                      level2Up = ifelse(UpUp==1|DownUp==1,"",level2Up),
                      level2Down = ifelse(UpDown==1|DownDown==1,"",level2Down)) %>%
        dplyr::mutate(Direction = ifelse(level1Up==1 | level1Down==1,paste(xrename, " only"),"-"),
                      Direction = ifelse(level2Up==1 | level2Down==1,paste(yrename, " only"),Direction),
                      Direction = ifelse(UpUp==1 | DownDown==1,"homodirectional",Direction),
                      Direction = ifelse(UpDown==1 | DownUp==1,"opposite change",Direction)) %>%
        dplyr::mutate(similar=1/abs(YlogFC-XlogFC),bigger = abs(YlogFC+XlogFC)) %>%
        dplyr::mutate(similar=similar/max(similar),bigger = bigger/max(bigger)) %>% rowwise() %>% dplyr::mutate(weight = weighted.mean(c(similar,bigger),c(1,2)))
        
      
      minx = round(min(xdf$logFC),0) - 1
      miny = round(min(ydf$logFC),0) - 1
      maxx = round(max(xdf$logFC),0) + 1
      maxy = round(max(ydf$logFC),0) + 1
      
      resultMatrix$Direction <- factor(resultMatrix$Direction,levels = c("-",paste(xrename, " only"),paste(yrename, " only"),"homodirectional","opposite change"))
      
      DEG_only <- subset(resultMatrix,level1Up==1|level1Down==1|level2Up==1|level2Down==1|UpUp==1|DownUp==1|UpDown==1|DownDown==1)
      
      Spearman_all = round(cor.test(resultMatrix$XlogFC, resultMatrix$YlogFC, method="spearman")$estimate, 2)
      Spearman_DEG = round(cor.test(DEG_only$XlogFC, DEG_only$YlogFC, method="spearman")$estimate, 2)
      Spearman_hth = round(cor.test(subset(resultMatrix,hth!="")$XlogFC, subset(resultMatrix,hth!="")$YlogFC, method="spearman")$estimate, 2)
      
      DEG_list <- data.frame(table((DEG_only[,c(10:20)] %>% tidyr::gather(key=level,value=value,-Gene_name,-Direction) %>% dplyr::filter(value!=""))$level))
      DEG_list <- merge(DEG_list,data.frame("Var1" = c("DownUp","UpDown","DownDown","UpUp","level2Down","level2Up","level1Down","level1Up")),by="Var1",all=T)
      DEG_list$Var1 <- factor(DEG_list$Var1,levels = c("DownUp","UpDown","DownDown","UpUp","level2Down","level2Up","level1Down","level1Up"))
      
      
      sp <- ggplot(resultMatrix, aes(x=XlogFC, y=YlogFC, label=Gene_name, color=Direction)) + geom_point(size=1,alpha=0.4)  +
        scale_color_manual(values = c("grey","steelblue3","gold1","chartreuse3","red2")) +
        geom_vline(xintercept = c(FC.threshold,-FC.threshold), colour = "gray65",size = 0.5, linetype = "dotted") +
        geom_hline(yintercept = c(FC.threshold,-FC.threshold), colour = "gray65",size = 0.5, linetype = "dotted") +
        annotate("text",x=maxx,y=maxy,label = "+/+") + annotate("text",x=maxx,y=maxy-0.5,label = DEG_list[which(DEG_list$Var1=="UpUp"),2],size=2) +
        annotate("text",x=maxx,y=miny,label = "+/-") + annotate("text",x=maxx,y=miny+0.5,label = DEG_list[which(DEG_list$Var1=="UpDown"),2],size=2) +
        annotate("text",x=minx,y=miny,label = "-/-") + annotate("text",x=minx,y=miny+0.5,label = DEG_list[which(DEG_list$Var1=="DownDown"),2],size=2) +
        annotate("text",x=minx,y=maxy,label = "-/+") + annotate("text",x=minx,y=maxy-0.5,label = DEG_list[which(DEG_list$Var1=="DownUp"),2],size=2) +
        annotate("text",x=0,y=maxy,label = DEG_list[which(DEG_list$Var1=="level2Up"),2],size=2) +
        annotate("text",x=0,y=miny,label = DEG_list[which(DEG_list$Var1=="level2Down"),2],size=2) +
        annotate("text",y=0,x=maxx,label = DEG_list[which(DEG_list$Var1=="level1Up"),2],size=2) +
        annotate("text",y=0,x=minx,label = DEG_list[which(DEG_list$Var1=="level1Down"),2],size=2) +
        xlim(minx,maxx) + ylim(miny,maxy) + 
        labs(x = xrename, y = yrename ) +
        annotate("text", x = maxx - 0.5, y = miny + 1 , label = paste("Spearman all:", Spearman_all," "),hjust=1) +
        annotate("text", x = maxx - 0.5, y = miny + 0.5, label = paste("Spearman DEGs:",Spearman_DEG," "),hjust=1) +
        theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
              axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
              axis.text.x = element_text(hjust = 1),
              plot.title = element_text(size=20, face="bold",hjust=0.5))
      
      hsp <- ggplot(resultMatrix, aes(x=XlogFC, y=YlogFC, label=Gene_name, color=hth)) + geom_point(size=1,alpha=0.4)  +
        scale_color_manual(values = c("grey","steelblue3","gold1","chartreuse3","red2")) +
        geom_vline(xintercept = c(FC.threshold,-FC.threshold), colour = "gray65",size = 0.5, linetype = "dotted") +
        geom_hline(yintercept = c(FC.threshold,-FC.threshold), colour = "gray65",size = 0.5, linetype = "dotted") +
        annotate("text",x=maxx,y=maxy,label = "+/+") + annotate("text",x=maxx,y=maxy-0.5,label = DEG_list[which(DEG_list$Var1=="UpUp"),2],size=2) +
        annotate("text",x=maxx,y=miny,label = "+/-") + annotate("text",x=maxx,y=miny+0.5,label = DEG_list[which(DEG_list$Var1=="UpDown"),2],size=2) +
        annotate("text",x=minx,y=miny,label = "-/-") + annotate("text",x=minx,y=miny+0.5,label = DEG_list[which(DEG_list$Var1=="DownDown"),2],size=2) +
        annotate("text",x=minx,y=maxy,label = "-/+") + annotate("text",x=minx,y=maxy-0.5,label = DEG_list[which(DEG_list$Var1=="DownUp"),2],size=2) +
        annotate("text",x=0,y=maxy,label = DEG_list[which(DEG_list$Var1=="level2Up"),2],size=2) +
        annotate("text",x=0,y=miny,label = DEG_list[which(DEG_list$Var1=="level2Down"),2],size=2) +
        annotate("text",y=0,x=maxx,label = DEG_list[which(DEG_list$Var1=="level1Up"),2],size=2) +
        annotate("text",y=0,x=minx,label = DEG_list[which(DEG_list$Var1=="level1Down"),2],size=2) +
        xlim(minx,maxx) + ylim(miny,maxy) + 
        labs(x = xrename, y = yrename ) +
        annotate("text", x = maxx - 0.5, y = miny + 1 , label = paste("Spearman all:", Spearman_all," "),hjust=1) +
        annotate("text", x = maxx - 0.5, y = miny + 0.5, label = paste("Spearman Hth:",Spearman_hth," "),hjust=1) +
        theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
              axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
              axis.text.x = element_text(hjust = 1),
              plot.title = element_text(size=20, face="bold",hjust=0.5))
      
      # sp
      ps <- ggplotly(sp,tooltip = "label")
      psh <- ggplotly(hsp,tooltip = "label")
      
      bp <- ggplot(DEG_list,aes(x=Var1,y=Freq,fill=Var1)) + geom_bar(stat = "identity") + coord_flip() + 
        labs(x="Level", y="# of genes", title = paste0(xrename," - ",yrename)) +
        scale_fill_manual(values = c("#f46666","#be0000","#a3e166","#51a400","#ffe766","#e5c100","#83b4dc","#4785b8"),guide="none") +
        theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
              axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
              axis.text.x = element_text(hjust = 1),
              axis.text.y = element_text(hjust = 1, angle = 45),
              plot.title = element_text(size=20, face="bold",hjust=0.5))
      
      list1 <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(XpValue<significance.threshold & XlogFC > FC.threshold) %>% dplyr::select(Gene_name)
      list2 <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(YpValue<significance.threshold & YlogFC > FC.threshold) %>% dplyr::select(Gene_name)
      list3 <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(XpValue<significance.threshold & XlogFC < -FC.threshold) %>% dplyr::select(Gene_name)
      list4 <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(YpValue<significance.threshold & YlogFC < -FC.threshold) %>% dplyr::select(Gene_name)
      
      Nlist1 <- length(unlist(list1))
      Nlist2 <- length(unlist(list2))
      Nlist3 <- length(unlist(list3))
      Nlist4 <- length(unlist(list4))
      
      Nlist1_2 <- sum(list1$Gene_name %in% list2$Gene_name)
      Nlist3_4 <- sum(list3$Gene_name %in% list4$Gene_name)
      
      fitUP <- euler(c("exp1" = Nlist1-Nlist1_2, "exp2" = Nlist2-Nlist1_2, "exp1&exp2" = Nlist1_2), shape = "ellipse")
      fitDOWN <- euler(c("exp1" = Nlist3-Nlist3_4, "exp2" = Nlist4-Nlist3_4, "exp1&exp2" = Nlist3_4), shape = "ellipse")
      
      plotUP <- plot(fitUP, quantities = list(cex = 1), fills = list(fill = c("gray67", "darkgoldenrod1"), alpha = 0.5), legend = list(labels = c(xrename, yrename)))
      plotDOWN <- plot(fitDOWN, quantities = list(cex = 1), fills = list(fill = c("gray67", "darkgoldenrod1"), alpha = 0.5), legend = list(labels = c(xrename, yrename)))
      
      list1a <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(level1Up==1 | UpUp==1 | UpDown==1) %>% dplyr::select(Gene_name)
      list2a <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(level2Up==1 | UpUp==1 | DownUp==1) %>% dplyr::select(Gene_name)
      list3a <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(level1Down==1 | DownUp==1 | DownDown==1) %>% dplyr::select(Gene_name)
      list4a <- resultMatrix %>% dplyr::filter(!is.na(Gene_name)) %>% dplyr::filter(level2Down==1 | UpDown==1 | DownDown==1) %>% dplyr::select(Gene_name)
      
      Nalist1 <- length(unlist(list1a))
      Nalist2 <- length(unlist(list2a))
      Nalist3 <- length(unlist(list3a))
      Nalist4 <- length(unlist(list4a))
      
      Nalist1_2 <- sum(list1a$Gene_name %in% list2a$Gene_name)
      Nalist3_4 <- sum(list3a$Gene_name %in% list4a$Gene_name)
      
      fitUPa <- euler(c("exp1" = Nalist1-Nalist1_2, "exp2" = Nalist2-Nalist1_2, "exp1&exp2" = Nalist1_2), shape = "ellipse")
      fitDOWNa <- euler(c("exp1" = Nalist3-Nalist3_4, "exp2" = Nalist4-Nalist3_4, "exp1&exp2" = Nalist3_4), shape = "ellipse")
      
      plotUPa <- plot(fitUPa, quantities = list(cex = 1), fills = list(fill = c("dodgerblue4", "gainsboro"), alpha = 0.5), legend = list(labels = c(xrename, yrename)))
      plotDOWNa <- plot(fitDOWNa, quantities = list(cex = 1), fills = list(fill = c("dodgerblue4", "gainsboro"), alpha = 0.5), legend = list(labels = c(xrename, yrename)))
      
      xlabel=a[1]
      ylabel=a[2]
      
      condition1 <-paste(gsub("\\[|\\]","",unlist(strsplit(xlabel," - ")[[1]])),collapse = "|")
      condition2 <-paste(gsub("\\[|\\]","",unlist(strsplit(ylabel," - ")[[1]])),collapse = "|")
      
      total <- merge(resultMatrix, rpkm[,c(1:2,grep(condition1,names(rpkm)),grep(condition2,names(rpkm)))], by.x="Gene_name" , by.y="external_gene_name",all.x=T)
      total <- merge(total, tx5fly[,c(1,5:6)], by.x="Gene_ID" , by.y="flybase_gene_id",all.x=T)
      
      if(name_DEG=="DEG_RPKM"){
        write.csv(total,paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/excel/",xrename,"_VS_",yrename,".csv"),row.names = F)
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/pdf/",xrename,"_VS_",yrename,".pdf"),width = 20,height = 10,useDingbats=FALSE)
        print(sp)
        print(bp)
        print(hsp)
        print(cowplot::plot_grid(plotUPa,plotDOWNa,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [logFC>0.58]", "Downregulated gene [logFC<-0.58]")))
        print(cowplot::plot_grid(plotUP,plotDOWN,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [pValue(single experiment) < 0.05, logFC>0.58]", "pValue(single experiment) < 0.05, Downregulated gene [logFC<-0.58]")))
        dev.off()
        htmlwidgets::saveWidget(ps, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/interactive/",xrename,"_VS_",yrename,".html"))
        htmlwidgets::saveWidget(psh, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_comparison/hth_interactive/",xrename,"_VS_",yrename,".html"))
      } else if(name_DEG=="DEG_1RPKM") {
        write.csv(total,paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_comparison/excel/",xrename,"_VS_",yrename,".csv"),row.names = F)
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_comparison/pdf/",xrename,"_VS_",yrename,".pdf"),width = 15,height = 10,useDingbats=FALSE)
        print(sp)
        print(bp)
        print(hsp)
        print(cowplot::plot_grid(plotUPa,plotDOWNa,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [logFC>0.58]", "Downregulated gene [logFC<-0.58]")))
        print(cowplot::plot_grid(plotUP,plotDOWN,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [pValue(single experiment) < 0.05, logFC>0.58]", "pValue(single experiment) < 0.05, Downregulated gene [logFC<-0.58]")))
        dev.off()
        htmlwidgets::saveWidget(ps, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_comparison/interactive/",xrename,"_VS_",yrename,".html"))
        htmlwidgets::saveWidget(psh, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_comparison/hth_interactive/",xrename,"_VS_",yrename,".html"))
      } else {
        write.csv(total,paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_comparison/excel/",xrename,"_VS_",yrename,".csv"),row.names = F)
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_comparison/pdf/",xrename,"_VS_",yrename,".pdf"),width = 15,height = 10,useDingbats=FALSE)
        print(sp)
        print(bp)
        print(hsp)
        print(cowplot::plot_grid(plotUPa,plotDOWNa,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [logFC>0.58]", "Downregulated gene [logFC<-0.58]")))
        print(cowplot::plot_grid(plotUP,plotDOWN,align = "v",axis = "l",nrow = 2,labels = c("Upregulated gene [pValue(single experiment) < 0.05, logFC>0.58]", "pValue(single experiment) < 0.05, Downregulated gene [logFC<-0.58]")))
        dev.off()
        htmlwidgets::saveWidget(ps, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_comparison/interactive/",xrename,"_VS_",yrename,".html"))
        htmlwidgets::saveWidget(psh, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_comparison/hth_interactive/",xrename,"_VS_",yrename,".html"))
      }
      
      k = k + 1
    }
  }
}



