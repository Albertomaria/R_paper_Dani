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

iab8 <- unique(Predicted_Targets_score[grep("iab8",Predicted_Targets_score$miRNA),]$gene_ID)
hth_11 <- unique(Hth[grep("Confidence",Hth$Seed),]$gene_ID) 
hth_13 <- unique(Hth[grep("Confidence",Hth$Seed,invert = T),]$gene_ID) 

DEG_RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T) %>%
  dplyr::mutate(iab8 = ifelse(flybase_gene_id %in% iab8, "iab8","No_target" ), hth11 = ifelse(flybase_gene_id %in% hth_11, "hth11","No_target" ), hth13 = ifelse(flybase_gene_id %in% hth_13, "hth13","No_target" ))
DEG_1RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/1RPKM/DEG_filter_1RPKM.csv",header = T) %>% 
  dplyr::mutate(iab8 = ifelse(flybase_gene_id %in% iab8, "iab8","No_target" ), hth11 = ifelse(flybase_gene_id %in% hth_11, "hth11","No_target" ), hth13 = ifelse(flybase_gene_id %in% hth_13, "hth13","No_target" ))
DEG_5RPKM <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/5RPKM/DEG_filter_5RPKM.csv",header = T) %>% 
  dplyr::mutate(iab8 = ifelse(flybase_gene_id %in% iab8, "iab8","No_target" ), hth11 = ifelse(flybase_gene_id %in% hth_11, "hth11","No_target" ), hth13 = ifelse(flybase_gene_id %in% hth_13, "hth13","No_target" ))


colnames(DEG_RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_RPKM[,-c(1:2,39:41)])))),"iab8","hth[2011]","hth[2013")
colnames(DEG_1RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_1RPKM[,-c(1:2,39:41)])))),"iab8","hth[2011]","hth[2013")
colnames(DEG_5RPKM) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_5RPKM[,-c(1:2,39:41)])))),"iab8","hth[2011]","hth[2013")

DEG_RPKM_list <- unique(DEG_RPKM %>% tidyr::gather(key = "experiment",value="value",-c(1:2,39:41) ) %>% tidyr::separate(experiment,into=c("Statistic","Experiment"),sep="\\[") %>%
                          dplyr::mutate(Experiment = gsub("]","",Experiment)) %>%
                          tidyr::spread(key="Statistic",value="value") %>% tidyr::gather(key = "Condition",value="Target",-c(1:2,6:10)) %>% 
                          dplyr::select(-Condition)) %>%
  dplyr::filter(Experiment %in% c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8"))
DEG_RPKM_list$Target = factor(DEG_RPKM_list$Target,levels = c("No_target","iab8","hth11","hth13"))
DEG_RPKM_list$Experiment = factor(DEG_RPKM_list$Experiment,levels = c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8"))

DEG_1RPKM_list <- unique(DEG_1RPKM %>% tidyr::gather(key = "experiment",value="value",-c(1:2,39:41) ) %>% tidyr::separate(experiment,into=c("Statistic","Experiment"),sep="\\[") %>%
                           dplyr::mutate(Experiment = gsub("]","",Experiment)) %>%
                           tidyr::spread(key="Statistic",value="value") %>% tidyr::gather(key = "Condition",value="Target",-c(1:2,6:10)) %>% 
                           dplyr::select(-Condition)) %>%
  dplyr::filter(Experiment %in% c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8")) #%>%
DEG_1RPKM_list$Target = factor(DEG_1RPKM_list$Target,levels = c("No_target","iab8","hth11","hth13"))
DEG_1RPKM_list$Experiment = factor(DEG_1RPKM_list$Experiment,levels = c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8"))

DEG_5RPKM_list <- unique(DEG_5RPKM %>% tidyr::gather(key = "experiment",value="value",-c(1:2,39:41) ) %>% tidyr::separate(experiment,into=c("Statistic","Experiment"),sep="\\[") %>%
                           dplyr::mutate(Experiment = gsub("]","",Experiment)) %>%
                           tidyr::spread(key="Statistic",value="value") %>% tidyr::gather(key = "Condition",value="Target",-c(1:2,6:10)) %>% 
                           dplyr::select(-Condition)) %>%
  dplyr::filter(Experiment %in% c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8")) #%>%
DEG_5RPKM_list$Target = factor(DEG_5RPKM_list$Target,levels = c("No_target","iab8","hth11","hth13"))
DEG_5RPKM_list$Experiment = factor(DEG_5RPKM_list$Experiment,levels = c("miR_iab8 - CS_iab8","hthBS_iab8 - CS_iab8","miR_iab8 - hthBS_iab8"))


pdf("./boxplot_iab8_comparison_DEG_zoom2DEG.pdf",width = 20,height = 10,useDingbats = FALSE)
ggplot(subset(DEG_RPKM_list,PValue<0.05),aes(x=Experiment,y=logFC,fill=Target)) + geom_boxplot() + #ylim(-2,2) +
  labs(title="Distribution DEG [no Filter]") +
  theme(aspect.ratio=1/3, panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), 
        axis.title = element_text(size=20,face = "bold"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1,angle=45,face="bold"),
        plot.title = element_text(size=20, face="bold",hjust=0.5))

ggplot(subset(DEG_1RPKM_list,PValue<10),aes(x=Experiment,y=abs(logFC),fill=Target)) + geom_violin() + ylim(0,2) +
  labs(title="Distribution DEG [Filter RPKM >1]") +
  theme(aspect.ratio=1/3, panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), 
        axis.title = element_text(size=20,face = "bold"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1,angle=45,face="bold"),
        plot.title = element_text(size=20, face="bold",hjust=0.5))

ggplot(DEG_5RPKM_list,aes(x=Experiment,y=logFC,fill=Target)) + geom_boxplot() + ylim(-2,2) +
  labs(title="Distribution DEG [Filter RPKM >5]") +
  theme(aspect.ratio=1/3, panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), 
        axis.title = element_text(size=20,face = "bold"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1,angle=45,face="bold"),
        plot.title = element_text(size=20, face="bold",hjust=0.5))
dev.off()

