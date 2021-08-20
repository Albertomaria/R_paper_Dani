setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")

library(tidyverse)
# library(edgeR)
library(biomaRt)
#library(goSTAG)
library(org.Dm.eg.db)
library(readr)
library(readxl)
library(plotly)
library(matrixStats)


########## Import Ensembl (bioMart) information
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
#View(listAttributes(mart))
tx4fly <- getBM(attributes = c("external_gene_name","flybase_gene_id","transcript_length","entrezgene"),mart=mart)

# to_check=c("Ubx","abd-A","Abd-B","hth","exd","Tdc2","Ilp7")

Predicted_Targets_score <- read_delim("E:/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")
Predicted_Targets_score <- Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(1:3,5:6,11,7:8) %>% 
  dplyr::filter(miRNA %in% c("dme-miR-279-3p", "dme-miR-315-5p","dme-miR-iab-4-3p","dme-miR-iab-4-5p","dme-miR-iab-8-5p")) %>%
  dplyr::mutate(miRNA = gsub("iab-4","iab4",miRNA),miRNA = gsub("iab-8","iab8",miRNA),miRNA = gsub("dme-miR-","",miRNA),Site_Type=ifelse(Site_Type==1,"7mer-a1",ifelse(Site_Type==2,"7mer-m8","8mer"))) %>% 
  tidyr::unite("miRNA",c("miRNA","Site_Type")) #%>% tidyr::spread(key=miRNA,value=weighted_context_score)
Predicted_Targets_score

Slattery_2013 <- read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pgen.1003753.s001.XLSX",
                            sheet = "Hth ChIP Peaks", col_types = c("blank", "blank", "blank", "text", "blank", "blank", "blank", "blank", "blank", 
                                                                    "text", "blank", "blank", "blank", "blank", "blank", "text"), skip = 1)

colnames(Slattery_2013) <- c("EA=W","W>EA","EA>W")
Hth <- data.frame("gene_ID"=NULL,"Seed"=NULL)
for (col in (1:3)){
  x <- Slattery_2013[,col] %>% dplyr::mutate(Seed = names(Slattery_2013)[col]) %>% dplyr::select(gene_ID=1,2)
  Hth <- rbind(Hth,x)
}
Hth <- unique(Hth %>% tidyr::separate_rows(gene_ID,sep=";") %>% dplyr::group_by(gene_ID,Seed) %>% dplyr::summarise(Count=n()))

Hth <- rbind(Hth, 
             read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s007.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
               dplyr::select(gene_ID = 1) %>% tidyr::separate_rows(gene_ID,sep=";") %>% 
               dplyr::mutate(Seed="High-Confidence") %>% group_by(gene_ID, Seed) %>% dplyr::summarise(Count=n()),
             read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s009.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
               dplyr::select(gene_ID = 1) %>% tidyr::separate_rows(gene_ID,sep=";") %>% 
               dplyr::mutate(Seed="Medium-Confidence") %>% group_by(gene_ID, Seed) %>% dplyr::summarise(Count=n()))

conversion <- read.delim("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/unknown_conversion.txt", header=T)

Hth <- merge(conversion,Hth,by.x="Submitted_id",by.y="gene_ID",all.y=T) %>% select(gene_ID=2,gene_name=4,5,6)
Hth$gene_ID <- as.character(Hth$gene_ID)
Hth$gene_name <- as.character(Hth$gene_name)


rm(list = c("x","mart","col","Slattery_2013","conversion"))


########## Comulative plot (ECDF) FoldChange
seed <- list()

for (mir in levels(factor(Predicted_Targets_score$miRNA))){
  # assign(mir, Predicted_Targets_score %>% dplyr::filter(miRNA %in% mir))
  seed[[mir]] <- Predicted_Targets_score %>% dplyr::filter(miRNA %in% mir)
}

for (hth in levels(factor(Hth$Seed))){
  # assign(mir, Predicted_Targets_score %>% dplyr::filter(miRNA %in% mir))
  seed[[hth]] <- Hth %>% dplyr::filter(Seed %in% hth)
}

DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T)
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T)

colnames(DEG_RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:2,39:42)])))))
colnames(DEG_1RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:2,39:42)])))))
colnames(DEG_5RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:2,39:42)])))))

for(name_DEG in c("DEG_RPKM","DEG_1RPKM","DEG_5RPKM")){
  DEG_df <- get(name_DEG)
  
  for (name in names(DEG_df)[grep("logFC",names(DEG_df))]){
  comparison <- gsub("logFC","",name)
  ddf <- cbind(DEG_df[,c(1:2)],DEG_df[,grep(comparison,names(DEG_df),fixed=T)]) %>% 
    dplyr::rename("gene_ID"=1,"logFC"=3,"logCPM"=4,"F"=5,"PValue"=6) # %>% dplyr::filter(PValue<0.05) %>% dplyr::filter(logFC<2,logFC>-2)
  
  for(condition in names(seed)){
    if(class(unique(seed[[condition]][,1])) == "character"){
      gene <- data.frame("gene_ID"=unique(seed[[condition]][,1])) %>% dplyr::mutate(!!condition :="Yes")
    } else {
      gene <- unique(seed[[condition]][,1]) %>% dplyr::mutate(!!condition :="Yes")
    }
    ddf <- merge(ddf,gene,by="gene_ID",all.x=T)
  }
  
  ddf <- ddf %>% dplyr::mutate(No_iab4_site = ifelse(is.na(`iab4-5p_7mer-a1`) & is.na(`iab4-5p_7mer-m8`) & is.na(`iab4-5p_8mer`),"Yes",NA))
  ddf <- ddf %>% dplyr::mutate(No_iab8_site = ifelse(is.na(`iab8-5p_7mer-a1`) & is.na(`iab8-5p_7mer-m8`) & is.na(`iab8-5p_8mer`),"Yes",NA))
  
  df <- data.frame(matrix(NA, nrow = 22, ncol = 22))
  colnames(df) <- c(names(ddf)[7:28])
  rownames(df) <- c(names(ddf)[7:28])
  
  for (x in 1:length(df)){
    for (y in 1:length(df)){
      xn <- colnames(df)[x]
      yn <- colnames(df)[y]
      pv <- format(ks.test((ddf %>% dplyr::filter(!is.na(get(yn))))$logFC, (ddf %>% dplyr::filter(!is.na(get(xn))))$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
      if (pv == "0e+00"){pv = "< 2.2e-16"}
      df[y,x] <- pv
    }
  }
  
  ECDF4 <- ggplot(ddf , aes(logFC)) + stat_ecdf(geom = "step",size=1) +  
    stat_ecdf(data = subset(ddf,!is.na(No_iab4_site)),aes(logFC), geom = "step",size=1) + 
      annotate("text",x=2.5,y=0.03,label=paste0("No iab4 target [",table(ddf$No_iab4_site),"]"),hjust = 1) +
      stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_7mer-a1`)),aes(logFC), geom = "step",color="#008744",size=1) + 
      annotate("text",x=2.5,y=0.06,label=paste0("iab4-5p 7mer-a1 [",table(ddf$`iab4-5p_7mer-a1`),"] pValue = ",df["No_iab4_site","iab4-5p_7mer-a1"]),hjust = 1,color="#008744") +
      stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_7mer-m8`)),aes(logFC), geom = "step",color="#ffa700",size=1) + 
      annotate("text",x=2.5,y=0.09,label=paste0("iab4-5p 7mer-m8 [",table(ddf$`iab4-5p_7mer-m8`),"] pValue = ",df["No_iab4_site","iab4-5p_7mer-m8"]),hjust = 1,color="#ffa700") +
      stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_8mer`)),aes(logFC), geom = "step",color="#0058e7",size=1) + 
      annotate("text",x=2.5,y=0.12,label=paste0("iab4-5p 8mer [",table(ddf$`iab4-5p_8mer`),"] pValue = ",df["No_iab4_site","iab4-5p_8mer"]),hjust = 1,color="#0058e7") +
      xlim(-2.5,2.5) + labs(x="\nlog2[mRNA Fold Change]",y="Cumulative fraction \n",title=comparison) + 
      theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5))
    
    ECDF8 <- ggplot(ddf , aes(logFC)) + stat_ecdf(geom = "step",size=1) +  
      stat_ecdf(data = subset(ddf,!is.na(No_iab8_site)),aes(logFC), geom = "step",size=1) + 
      annotate("text",x=2.5,y=0.03,label=paste0("No iab8 target [",table(ddf$No_iab8_site),"]"),hjust = 1) +
      stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_7mer-a1`)),aes(logFC), geom = "step",color="#008744",size=1) + 
      annotate("text",x=2.5,y=0.06,label=paste0("iab8-5p 7mer-a1 [",table(ddf$`iab8-5p_7mer-a1`),"] pValue = ",df["No_iab8_site","iab8-5p_7mer-a1"]),hjust = 1,color="#008744") +
      stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_7mer-m8`)),aes(logFC), geom = "step",color="#ffa700",size=1) + 
      annotate("text",x=2.5,y=0.09,label=paste0("iab8-5p 7mer-m8 [",table(ddf$`iab8-5p_7mer-m8`),"] pValue = ",df["No_iab8_site","iab8-5p_7mer-m8"]),hjust = 1,color="#ffa700") +
      stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_8mer`)),aes(logFC), geom = "step",color="#0058e7",size=1) + 
      annotate("text",x=2.5,y=0.12,label=paste0("iab8-5p 8mer [",table(ddf$`iab8-5p_8mer`),"] pValue = ",df["No_iab8_site","iab8-5p_8mer"]),hjust = 1,color="#0058e7") +
      xlim(-2.5,2.5) + labs(x="\nlog2[mRNA Fold Change]",y="Cumulative fraction \n",title=comparison) + 
      theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5))
    
    if(name_DEG=="DEG_RPKM"){
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab4/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF4)
      dev.off()
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab8/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF8)
      dev.off()
    } else if(name_DEG=="DEG_1RPKM") {
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab4/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF4)
      dev.off()
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab8/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF8)
      dev.off()
    } else {
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab4/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF4)
      dev.off()
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab8/logFC/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(ECDF8)
      dev.off()
    }
  } 
}

rm(list=c(ls()[-c(which(ls()=="Predicted_Targets_score"),which(ls()=="Hth"),which(ls()=="seed"))]))

########## Comulative plot (ECDF) RPKM expression single experiment

RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]

for(name_RPKM in c("RPKM","RPKM1","RPKM5")){
  RPKM_annotate <- get(name_RPKM)
  
  for (name in names(RPKM_annotate)[grep("mean",names(RPKM_annotate))]){
  comparison <- gsub("_mean","",name)
  ddf <- cbind(RPKM_annotate[,c(1:2)],RPKM_annotate[,grep(comparison,names(RPKM_annotate),fixed=T)]) %>% 
    dplyr::rename("gene_ID"=1,"RPKM"=3,"SD"=4) %>% dplyr::mutate(logRPKM=log2(RPKM)) # %>% dplyr::filter(SD<0.05) %>% dplyr::filter(logFC<2,logFC>-2)
  
  for(condition in names(seed)){
    if(class(unique(seed[[condition]][,1])) == "character"){
      gene <- data.frame("gene_ID"=unique(seed[[condition]][,1])) %>% dplyr::mutate(!!condition :="Yes")
    } else {
      gene <- unique(seed[[condition]][,1]) %>% dplyr::mutate(!!condition :="Yes")
    }
    ddf <- merge(ddf,gene,by="gene_ID",all.x=T)
  }
  
  ddf <- ddf %>% dplyr::mutate(No_iab4_site = ifelse(is.na(`iab4-5p_7mer-a1`) & is.na(`iab4-5p_7mer-m8`) & is.na(`iab4-5p_8mer`),"Yes",NA))
  ddf <- ddf %>% dplyr::mutate(No_iab8_site = ifelse(is.na(`iab8-5p_7mer-a1`) & is.na(`iab8-5p_7mer-m8`) & is.na(`iab8-5p_8mer`),"Yes",NA))
  
  df <- data.frame(matrix(NA, nrow = 22, ncol = 22))
  colnames(df) <- c(names(ddf)[6:27])
  rownames(df) <- c(names(ddf)[6:27])
  
  for (x in 1:length(df)){
    for (y in 1:length(df)){
      xn <- colnames(df)[x]
      yn <- colnames(df)[y]
      pv <- format(ks.test((ddf %>% dplyr::filter(!is.na(get(yn))))$logRPKM, (ddf %>% dplyr::filter(!is.na(get(xn))))$logRPKM, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
      if (pv == "0e+00"){pv = "< 2.2e-16"}
      df[y,x] <- pv
    }
  }
  
  ECDF4 <- ggplot(ddf , aes(logRPKM)) + stat_ecdf(geom = "step",size=1) +  
    stat_ecdf(data = subset(ddf,!is.na(No_iab4_site)),aes(logRPKM), geom = "step",size=1) + 
    annotate("text",x=15,y=0.03,label=paste0("No iab4 target [",table(ddf$No_iab4_site),"]"),hjust = 1) +
    stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_7mer-a1`)),aes(logRPKM), geom = "step",color="#008744",size=1) + 
    annotate("text",x=15,y=0.06,label=paste0("iab4-5p 7mer-a1 [",table(ddf$`iab4-5p_7mer-a1`),"] pValue = ",df["No_iab4_site","iab4-5p_7mer-a1"]),hjust = 1,color="#008744") +
    stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_7mer-m8`)),aes(logRPKM), geom = "step",color="#ffa700",size=1) + 
    annotate("text",x=15,y=0.09,label=paste0("iab4-5p 7mer-m8 [",table(ddf$`iab4-5p_7mer-m8`),"] pValue = ",df["No_iab4_site","iab4-5p_7mer-m8"]),hjust = 1,color="#ffa700") +
    stat_ecdf(data = subset(ddf,!is.na(`iab4-5p_8mer`)),aes(logRPKM), geom = "step",color="#0058e7",size=1) + 
    annotate("text",x=15,y=0.12,label=paste0("iab4-5p 8mer [",table(ddf$`iab4-5p_8mer`),"] pValue = ",df["No_iab4_site","iab4-5p_8mer"]),hjust = 1,color="#0058e7") +
    xlim(-10,15) + labs(x="\nlog2[mRNA RPKM]",y="Cumulative fraction \n",title=comparison) + 
    theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5))
  
  ECDF8 <- ggplot(ddf , aes(logRPKM)) + stat_ecdf(geom = "step",size=1) +  
    stat_ecdf(data = subset(ddf,!is.na(No_iab8_site)),aes(logRPKM), geom = "step",size=1) + 
    annotate("text",x=15,y=0.03,label=paste0("No iab8 target [",table(ddf$No_iab8_site),"]"),hjust = 1) +
    stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_7mer-a1`)),aes(logRPKM), geom = "step",color="#008744",size=1) + 
    annotate("text",x=15,y=0.06,label=paste0("iab8-5p 7mer-a1 [",table(ddf$`iab8-5p_7mer-a1`),"] pValue = ",df["No_iab8_site","iab8-5p_7mer-a1"]),hjust = 1,color="#008744") +
    stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_7mer-m8`)),aes(logRPKM), geom = "step",color="#ffa700",size=1) + 
    annotate("text",x=15,y=0.09,label=paste0("iab8-5p 7mer-m8 [",table(ddf$`iab8-5p_7mer-m8`),"] pValue = ",df["No_iab8_site","iab8-5p_7mer-m8"]),hjust = 1,color="#ffa700") +
    stat_ecdf(data = subset(ddf,!is.na(`iab8-5p_8mer`)),aes(logRPKM), geom = "step",color="#0058e7",size=1) + 
    annotate("text",x=15,y=0.12,label=paste0("iab8-5p 8mer [",table(ddf$`iab8-5p_8mer`),"] pValue = ",df["No_iab8_site","iab8-5p_8mer"]),hjust = 1,color="#0058e7") +
    xlim(-10,15) + labs(x="\nlog2[RPKM]",y="Cumulative fraction \n",title=comparison) + 
    theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5))
  
  if(name_RPKM=="RPKM"){
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab4/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF4)
    dev.off()
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab8/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF8)
    dev.off()
  } else if(name_RPKM=="RPKM1") {
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab4/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF4)
    dev.off()
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab8/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF8)
    dev.off()
  } else {
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab4/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF4)
    dev.off()
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab8/RPKM/",comparison,".pdf"),width = 10,height = 10,useDingbats=FALSE)
    print(ECDF8)
    dev.off()
  
  }
  }
}

rm(list=c(ls()[-c(which(ls()=="Predicted_Targets_score"),which(ls()=="Hth"),which(ls()=="seed"))]))

########## Comulative plot (ECDF) RPKM expression comparison

RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM1 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/RPKM_1_per_genotype.csv", header = T)[,-c(3:10)]
RPKM5 <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/RPKM_5_per_genotype.csv", header = T)[,-c(3:10)]

for(name_RPKM in c("RPKM","RPKM1","RPKM5")){
  RPKM_annotate <- get(name_RPKM)
  
  ddf <- cbind(RPKM_annotate[,c(1:2)],RPKM_annotate[,grep("mean",names(RPKM_annotate),fixed=T)]) 
  colnames(ddf) <- c("gene_ID","external_gene_name",gsub("_mean","",colnames(ddf[,-c(1:2)])))
  
  for(condition in names(seed)){
    if(class(unique(seed[[condition]][,1])) == "character"){
      gene <- data.frame("gene_ID"=unique(seed[[condition]][,1])) %>% dplyr::mutate(!!condition :="Yes")
      } else {
        gene <- unique(seed[[condition]][,1]) %>% dplyr::mutate(!!condition :="Yes")
        }
    ddf <- merge(ddf,gene,by="gene_ID",all.x=T)
  }
  
  ddf <- ddf %>% dplyr::mutate(No_iab4_site = ifelse(is.na(`iab4-5p_7mer-a1`) & is.na(`iab4-5p_7mer-m8`) & is.na(`iab4-5p_8mer`),"Yes",NA))
  ddf <- ddf %>% dplyr::mutate(iab4_site = ifelse(!is.na(`iab4-5p_7mer-a1`) | !is.na(`iab4-5p_7mer-m8`) | !is.na(`iab4-5p_8mer`),"Yes",NA))
  ddf <- ddf %>% dplyr::mutate(No_iab8_site = ifelse(is.na(`iab8-5p_7mer-a1`) & is.na(`iab8-5p_7mer-m8`) & is.na(`iab8-5p_8mer`),"Yes",NA))
  ddf <- ddf %>% dplyr::mutate(iab8_site = ifelse(!is.na(`iab8-5p_7mer-a1`) | !is.na(`iab8-5p_7mer-m8`) | !is.na(`iab8-5p_8mer`),"Yes",NA))
  ddf[3:8] <- log2(ddf[3:8]+1)
    
  for (x in 3:8){
    y = x + 1
    while (y <= 8){
      data <- ddf %>% dplyr::select(1,2,"logX"=x,"logY"=y,9:32)
      xn <- colnames(ddf)[x]
      yn <- colnames(ddf)[y]
      
      ECDF4 <- ggplot() +  
        stat_ecdf(data = subset(data,!is.na(iab4_site)),aes(logX), geom = "step",color="#666666",size=1) + 
        stat_ecdf(data = subset(data,!is.na(iab4_site)),aes(logY), geom = "step",color="#dedede",size=1) + 
        annotate("text",x=15,y=0.055,label=paste0("iab4-5p [",table(data$iab4_site),"]"),hjust = 1,color="#000000") +
        annotate("text",x=15,y=0.03,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(iab4_site))$logX,y=subset(data,!is.na(iab4_site))$logY, 
                                                                          alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#000000",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_7mer-a1`)),aes(logX), geom = "step",color="#329f69",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_7mer-a1`)),aes(logY), geom = "step",color="#006c36",size=1) + 
        annotate("text",x=15,y=0.115,label=paste0("iab4-5p 7mer-a1 [",table(data$`iab4-5p_7mer-a1`),"]"),hjust = 1,color="#008744") +
        annotate("text",x=15,y=0.09,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab4-5p_7mer-a1`))$logX,y=subset(data,!is.na(`iab4-5p_7mer-a1`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#008744",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_7mer-m8`)),aes(logX), geom = "step",color="#ffb832",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_7mer-m8`)),aes(logY), geom = "step",color="#cc8500",size=1) + 
        annotate("text",x=15,y=0.175,label=paste0("iab4-5p 7mer-8m [",table(data$`iab4-5p_7mer-a1`),"]"),hjust = 1,color="#ffa700") +
        annotate("text",x=15,y=0.15,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab4-5p_7mer-m8`))$logX,y=subset(data,!is.na(`iab4-5p_7mer-m8`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#ffa700",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_8mer`)),aes(logX), geom = "step",color="#3279eb",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab4-5p_8mer`)),aes(logY), geom = "step",color="#0046b8",size=1) + 
        annotate("text",x=15,y=0.235,label=paste0("iab4-5p 8mer [",table(data$`iab4-5p_8mer`),"]"),hjust = 1,color="#0058e7") +
        annotate("text",x=15,y=0.21,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab4-5p_8mer`))$logX,y=subset(data,!is.na(`iab4-5p_8mer`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#0058e7",size=3) +
        xlim(-10,15) + labs(x="\nlog2[RPKM]",y="Cumulative fraction \n",title=paste0("Compare ",xn," and ",yn)) + 
        theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
              axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
              axis.text.x = element_text(hjust = 1),
              plot.title = element_text(size=20, face="bold",hjust=0.5))
      
      ECDF8 <- ggplot() +  
        stat_ecdf(data = subset(data,!is.na(iab8_site)),aes(logX), geom = "step",color="#666666",size=1) + 
        stat_ecdf(data = subset(data,!is.na(iab8_site)),aes(logY), geom = "step",color="#dedede",size=1) + 
        annotate("text",x=15,y=0.055,label=paste0("iab4-5p [",table(data$iab8_site),"]"),hjust = 1,color="#000000") +
        annotate("text",x=15,y=0.03,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(iab8_site))$logX,y=subset(data,!is.na(iab8_site))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#000000",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_7mer-a1`)),aes(logX), geom = "step",color="#329f69",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_7mer-a1`)),aes(logY), geom = "step",color="#006c36",size=1) + 
        annotate("text",x=15,y=0.115,label=paste0("iab8-5p 7mer-a1 [",table(data$`iab8-5p_7mer-a1`),"]"),hjust = 1,color="#008744") +
        annotate("text",x=15,y=0.09,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab8-5p_7mer-a1`))$logX,y=subset(data,!is.na(`iab8-5p_7mer-a1`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#008744",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_7mer-m8`)),aes(logX), geom = "step",color="#ffb832",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_7mer-m8`)),aes(logY), geom = "step",color="#cc8500",size=1) + 
        annotate("text",x=15,y=0.175,label=paste0("iab8-5p 7mer-8m [",table(data$`iab8-5p_7mer-a1`),"]"),hjust = 1,color="#ffa700") +
        annotate("text",x=15,y=0.15,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab8-5p_7mer-m8`))$logX,y=subset(data,!is.na(`iab8-5p_7mer-m8`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#ffa700",size=3) +
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_8mer`)),aes(logX), geom = "step",color="#3279eb",size=1) + 
        stat_ecdf(data = subset(data,!is.na(`iab8-5p_8mer`)),aes(logY), geom = "step",color="#0046b8",size=1) + 
        annotate("text",x=15,y=0.235,label=paste0("iab8-5p 8mer [",table(data$`iab8-5p_8mer`),"]"),hjust = 1,color="#0058e7") +
        annotate("text",x=15,y=0.21,label=paste0("pValue = ",format(ks.test(subset(data,!is.na(`iab8-5p_8mer`))$logX,y=subset(data,!is.na(`iab8-5p_8mer`))$logY, 
                                                                            alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)),hjust = 1,color="#0058e7",size=3) +
        xlim(-10,15) + labs(x="\nlog2[RPKM]",y="Cumulative fraction \n",title=paste0("Compare ",xn," and ",yn)) + 
        theme(aspect.ratio=1, panel.background = element_rect(fill = "white", colour = "black"),
              axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
              axis.text.x = element_text(hjust = 1),
              plot.title = element_text(size=20, face="bold",hjust=0.5))
      
      if(name_RPKM=="RPKM"){
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab4/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF4)
        dev.off()
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/iab8/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF8)
        dev.off()
      } else if(name_RPKM=="RPKM1") {
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab4/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF4)
        dev.off()
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/ECDF/iab8/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF8)
        dev.off()
      } else {
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab4/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF4)
        dev.off()
        pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/ECDF/iab8/RPKM_compared//Compare ",xn," and ",yn,".pdf"),width = 10,height = 10,useDingbats=FALSE)
        print(ECDF8)
        dev.off()
      }
      y = y + 1 
    }
  }
}

rm(list=c(ls()[-c(which(ls()=="Predicted_Targets_score"),which(ls()=="Hth"),which(ls()=="seed"))]))

