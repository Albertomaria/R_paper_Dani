setwd("E:/L3_RNAseq")

library(tidyverse)
library(edgeR)
library(biomaRt)
library(mixOmics)
library(RColorBrewer)
library(goSTAG)
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

to_check=c("Ubx","abd-A","Abd-B","hth","exd","Tdc2","Ilp7")

Predicted_Targets_score <- read_delim("E:/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")
Predicted_Targets_score <- Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(1:3,5:6,11) %>% dplyr::filter(miRNA %in% c("dme-miR-279-3p", "dme-miR-315-5p")) %>%
  dplyr::mutate(Site_Type=ifelse(Site_Type==1,"7mer-a1",ifelse(Site_Type==2,"7mer-m8","8mer")),miRNA = gsub("dme-miR-","",miRNA)) %>% 
  group_by(gene_ID, Gene_Symbol, Transcript_ID, miRNA) %>% dplyr::slice(which.min(weighted_context_score)) %>% 
  tidyr::unite("miRNA",c("miRNA","Site_Type")) %>% 
  tidyr::spread(key=miRNA,value=weighted_context_score)

miR279 <- unique(Predicted_Targets_score %>% dplyr::ungroup() %>% dplyr::select(gene_ID,Gene_Symbol,dplyr::starts_with("279"))) %>% dplyr::rename("7mer.a1"=3, "7mer.m8"=4, "8mer"=5) %>%
  dplyr::filter(!is.na(`7mer.a1`) | !is.na(`7mer.m8`) | !is.na(`8mer`)) %>% tidyr::gather(key="Seed",value="Score",-gene_ID, -Gene_Symbol) %>% dplyr::filter(!is.na(Score))
miR279$Seed <- factor(miR279$Seed, levels = c("8mer","7mer.m8","7mer.a1"))
miR279$gene_ID <- as.character(miR279$gene_ID)

miR315 <- unique(Predicted_Targets_score %>% dplyr::ungroup() %>% dplyr::select(gene_ID,Gene_Symbol,dplyr::starts_with("315"))) %>% dplyr::rename("7mer.a1"=3, "7mer.m8"=4, "8mer"=5) %>%
  dplyr::filter(!is.na(`7mer.a1`) | !is.na(`7mer.m8`) | !is.na(`8mer`)) %>% tidyr::gather(key="Seed",value="Score",-gene_ID, -Gene_Symbol) %>% dplyr::filter(!is.na(Score))
miR315$Seed <- factor(miR279$Seed, levels = c("8mer","7mer.m8","7mer.a1"))
miR315$gene_ID <- as.character(miR279$gene_ID)

rm(list=c("mart"))

########## 
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
load("./L3_sequencing.RData")
RNA_count <- as.data.frame(SEQrna_3L$counts)
names <- NULL 
for (i in seq(1:8)){
  names[i] <- strsplit(names(RNA_count),"[.]")[[i]][7]} 

colnames(RNA_count) <- names

RNA_count <- tx4fly %>% dplyr::group_by(flybase_gene_id) %>% 
  dplyr::summarise(external_gene_name=dplyr::first(external_gene_name),entrezgene=dplyr::first(entrezgene),transcript_length=max(transcript_length)) %>% 
  merge(.,RNA_count,by.x="flybase_gene_id",by.y="row.names",all.y=T)

miR279_count <- RNA_count %>% dplyr::select(flybase_gene_id,external_gene_name, transcript_length,dplyr::contains("miR279"))
colnames(miR279_count) <- c("gene_ID","Gene_Symbol","transcript_length",gsub("L3CNS_","",colnames(miR279_count[,-c(1:3)])))

miR315_count <- RNA_count %>% dplyr::select(flybase_gene_id,external_gene_name, transcript_length,dplyr::contains("miR315"))
colnames(miR315_count) <- c("gene_ID","Gene_Symbol","transcript_length",gsub("L3CNS_","",colnames(miR315_count[,-c(1:3)])))

rm(list=c("SEQrna_3L","SEQrna_UTR_3L","i","names","RNA_count"))

########## RPKM

rpkm_279 <- cbind(miR279_count[1:3],rpkm(miR279_count[,4:7],gene.length=miR279_count$transcript_length))
rpkm_279_l <- rpkm_279 %>% tidyr::gather(key="Experiment",value = "rpkm",-gene_ID,-Gene_Symbol,-transcript_length) %>% dplyr::mutate(group = gsub('.{3}$', '', Experiment))
rpkm_315 <- cbind(miR315_count[1:3],rpkm(miR315_count[,4:7],gene.length=miR315_count$transcript_length))
rpkm_315_l <- rpkm_315 %>% tidyr::gather(key="Experiment",value = "rpkm",-gene_ID,-Gene_Symbol,-transcript_length) %>% dplyr::mutate(group = gsub('.{3}$', '', Experiment))

rpkm_279_annotate <- merge(rpkm_279,miR279[,-2],by="gene_ID",all.x=T) %>% dplyr::select(1:3,8:9,dplyr::everything())
rpkm_315_annotate <- merge(rpkm_315,miR315[,-2],by="gene_ID",all.x=T) %>% dplyr::select(1:3,8:9,dplyr::everything())
rpkm_279_annotate$Seed <- as.character(rpkm_279_annotate$Seed)
rpkm_279_annotate[is.na(rpkm_279_annotate)] <- ""
rpkm_315_annotate[is.na(rpkm_315_annotate)] <- ""

miR279_group <- factor(rpkm_279_l$group)
miR315_group <- factor(rpkm_315_l$group)

rpkm_miR279_mean <- rpkm_279 %>% dplyr::mutate(miR279_KO_mean=rowMeans(.[,4:5]),miR279_Res_mean=rowMeans(.[,6:7]),
                                               miR279_KO_SD=rowMeans(as.matrix(.[,4:5])),miR279_Res_SD=rowMeans(as.matrix(.[,6:7]))) %>%
  dplyr::select(1,2,3,8,9,10,11)
rpkm_miR315_mean <- rpkm_315 %>% dplyr::mutate(miR315_KO_mean=rowMeans(.[,4:5]),miR315_Res_mean=rowMeans(.[,6:7]),
                                               miR315_KO_SD=rowMeans(as.matrix(.[,4:5])),miR315_Res_SD=rowMeans(as.matrix(.[,6:7]))) %>%
  dplyr::select(1,2,3,8,9,10,11)

########## Preliminary Graphs

palet_boxplot <- c(rep("#ff467e",2),rep("#f01f1f",2),rep("#467eff",2),rep("#0b25ee",2))

pdf("./Preliminary_graph/Density_plot_RPKM_total.pdf",width = 10,height = 8)
ggplot(rpkm_279_l, aes(x=log2(rpkm),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(head(palet_boxplot,n=4))) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype",title="miR-279_3p") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

ggplot(rpkm_279_l, aes(x=log2(rpkm+1),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(head(palet_boxplot,n=4))) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype",title="miR-279_3p") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

ggplot(rpkm_315_l, aes(x=log2(rpkm),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(tail(palet_boxplot,n=4))) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype",title="miR-315_5p") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

ggplot(rpkm_315_l, aes(x=log2(rpkm+1),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(tail(palet_boxplot,n=4))) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype",title="miR-315_3p") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
dev.off()

pdf("./Preliminary_graph/Density_plot_RPKM_for_Experiment.pdf",width = 21,height = 8)
rpkm_279_l %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = head(palet_boxplot,n=4)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

rpkm_279_l %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm+1),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = head(palet_boxplot,n=4)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

rpkm_315_l %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = tail(palet_boxplot,n=4)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))

rpkm_279_l %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm+1),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = tail(palet_boxplot,n=4)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
dev.off()

pdf("./Preliminary_graph/Boxplot_RPKM_single_experiment.pdf",width = 10,height = 8)
ggplot(rpkm_279_l, aes(x=Experiment,y=log2(rpkm), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = head(palet_boxplot,n=4)) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rpkm_279_l, aes(x=Experiment,y=log2(rpkm+1), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM+1]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = head(palet_boxplot,n=4)) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rpkm_315_l, aes(x=Experiment,y=log2(rpkm), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = tail(palet_boxplot,n=4)) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rpkm_315_l, aes(x=Experiment,y=log2(rpkm+1), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM+1]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = tail(palet_boxplot,n=4)) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

########## filtering and edgeR object

miR279_group <- factor(gsub('.{3}$', '', names(miR279_count[,4:ncol(miR279_count)])))
miR315_group <- factor(gsub('.{3}$', '', names(miR315_count[,4:ncol(miR315_count)])))

keep_279 <-rpkm_279[which(rowSums(rpkm_279[,4:7]>=1)>=2),1]
keep_315 <-rpkm_315[which(rowSums(rpkm_315[,4:7]>=1)>=2),1]

miR279_filtered <- miR279_count %>% dplyr::filter(gene_ID %in% keep_279)
miR315_filtered <- miR315_count %>% dplyr::filter(gene_ID %in% keep_315)

miR279_edger <- DGEList(counts = miR279_filtered[,4:ncol(miR279_filtered)],group=miR279_group,genes=miR279_filtered[,1:3])
miR279_edger <- calcNormFactors(miR279_edger)

miR315_edger <- DGEList(counts = miR315_filtered[,4:ncol(miR315_filtered)],group=miR315_group,genes=miR315_filtered[,1:3])
miR315_edger <- calcNormFactors(miR315_edger)

#########    EdgeR analysis - DEG

design_279 <- model.matrix(~0+miR279_group)
colnames(design_279) <- levels(miR279_group)
design_315 <- model.matrix(~0+miR315_group)
colnames(design_315) <- levels(miR315_group)

miR279_edger <- estimateDisp(miR279_edger, design_279)
miR279_fit <- exactTest(miR279_edger, pair=c("miR279_KO", "miR279_Res"))

miR315_edger <- estimateDisp(miR315_edger, design_315)
miR315_fit <- exactTest(miR315_edger, pair=c("mir315_KO", "mir315_Res"))

##### MMD 279

miR279_DEG <- cbind(miR279_fit$genes[,1],miR279_fit$table) %>% dplyr::rename("gene_ID"=1)
miR279_DEG <- merge(merge(miR279_DEG,rpkm_miR279_mean[,c(1,5)],by="gene_ID",all.x=T),miR279[,-2],by="gene_ID",all.x=T)
  
sc_279 <- ggplot(miR279_DEG,aes(y=logFC,x=log2(miR279_Res_mean+1),label=gene_ID)) + geom_point(alpha=0.3,color="#dedede") +
  geom_point(data=miR279_DEG%>% dplyr::filter(!is.na(Seed)),aes(y=logFC,x=log2(miR279_Res_mean+1), fill=Seed), alpha=1, size=3,shape=21) +
  scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
  labs(y = "logFC[miR279_KO - miR279_Res]",x="log2[RPKM+1]") +
  geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
  geom_vline(xintercept = 1, size = 0.5, color="#990000" ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
        axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size=20, face="bold",hjust=0.5),
        legend.position="bottom")
  
vp_279 <- miR279_DEG %>% dplyr::mutate(exp="miR279",`log2[RPKM+1]`=log2(miR279_Res_mean+1))
ggplot(vp_279,aes(x= exp, y=log2(miR279_Res_mean+1))) + geom_violin(fill="#dedede") +
  coord_flip() +
  geom_hline(yintercept = 1, size = 0.5, color="#990000" ) +
  geom_hline(yintercept = 2.585, size = 0.5, color="#009900" ) +
  annotate("text",x=1.50,y=1.7,label = paste0("< 5 RPKM\n",table(vp_279$`log2[RPKM+1]`>2.585)[1]), color="#009900",size=3) +
  annotate("text",x=1.50,y=3.3,label = paste0("> 5 RPKM\n",table(vp_279$`log2[RPKM+1]`>2.585)[2]), color="#009900",size=3) +
  annotate("text",x=1.45,y=0.3,label = paste0("< 1 RPKM\n",table(vp_279$`log2[RPKM+1]`>1)[1]), color="#990000",size=3) +
  annotate("text",x=1.45,y=1.7,label = paste0("> 1 RPKM\n",table(vp_279$`log2[RPKM+1]`>1)[2]), color="#990000",size=3) + 
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_blank(),axis.title = element_blank(),
          plot.title = element_blank())
  
Volp_279 <- ggplot(miR279_DEG,aes(y=-log10(PValue),x=logFC)) + geom_point(alpha=0.3,color="#dedede") +
  geom_point(data=miR279_DEG%>% dplyr::filter(!is.na(Seed)),aes(y=-log10(PValue),x=logFC, fill=Seed), alpha=1, size=3,shape=21) +
  scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
  labs(y = "logFC[miR279_KO - miR279_Res]",x="logFC") +
  geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
  geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
        axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size=20, face="bold",hjust=0.5),
        legend.position="bottom")

##### MMD 315

miR315_DEG <- cbind(miR315_fit$genes[,1],miR315_fit$table) %>% dplyr::rename("gene_ID"=1)
miR315_DEG <- merge(merge(miR315_DEG,rpkm_miR315_mean[,c(1,5)],by="gene_ID",all.x=T),miR315[,-2],by="gene_ID",all.x=T)

sc_315 <- ggplot(miR315_DEG,aes(y=logFC,x=log2(miR315_Res_mean+1),label=gene_ID)) + geom_point(alpha=0.3,color="#dedede") +
  geom_point(data=miR315_DEG%>% dplyr::filter(!is.na(Seed)),aes(y=logFC,x=log2(miR315_Res_mean+1), fill=Seed), alpha=1, size=3,shape=21) +
  scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
  labs(y = "logFC[miR315_KO - miR315_Res]",x="log2[RPKM+1]") +
  geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
  geom_vline(xintercept = 1, size = 0.5, color="#990000" ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
        axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size=20, face="bold",hjust=0.5),
        legend.position="bottom")

vp_315 <- miR315_DEG %>% dplyr::mutate(exp="miR315",`log2[RPKM+1]`=log2(miR315_Res_mean+1))
ggplot(vp_315,aes(x= exp, y=log2(miR315_Res_mean+1))) + geom_violin(fill="#dedede") +
  coord_flip() +
  geom_hline(yintercept = 1, size = 0.5, color="#990000" ) +
  geom_hline(yintercept = 2.585, size = 0.5, color="#009900" ) +
  annotate("text",x=1.50,y=1.7,label = paste0("< 5 RPKM\n",table(vp_315$`log2[RPKM+1]`>2.585)[1]), color="#009900",size=3) +
  annotate("text",x=1.50,y=3.3,label = paste0("> 5 RPKM\n",table(vp_315$`log2[RPKM+1]`>2.585)[2]), color="#009900",size=3) +
  annotate("text",x=1.45,y=0.3,label = paste0("< 1 RPKM\n",table(vp_315$`log2[RPKM+1]`>1)[1]), color="#990000",size=3) +
  annotate("text",x=1.45,y=1.7,label = paste0("> 1 RPKM\n",table(vp_315$`log2[RPKM+1]`>1)[2]), color="#990000",size=3) + 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_blank(),axis.title = element_blank(),
        plot.title = element_blank())

Volp_315 <- ggplot(miR315_DEG,aes(y=-log10(PValue),x=logFC)) + geom_point(alpha=0.3,color="#dedede") +
  geom_point(data=miR315_DEG%>% dplyr::filter(!is.na(Seed)),aes(y=-log10(PValue),x=logFC, fill=Seed), alpha=1, size=3,shape=21) +
  scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
  labs(y = "logFC[miR315_KO - miR315_Res]",x="logFC") +
  geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
  geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
        axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size=20, face="bold",hjust=0.5),
        legend.position="bottom")

  
########## Comulative plot (ECDF)

DEG_df <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/3L/DEG_total.csv",header = T)
DEG_df_miR <- merge(DEG_df,attribute_papers[,c(1,7:8)],by = "flybase_gene_id",all.x=T)
colnames(DEG_df_miR) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_df_miR[,-c(1:2,11:12)])))),"miR315","miR279")
# DEG_df_miR <- merge(DEG_df_miR,miR_315_279[,-2],by.x="flybase_gene_id",by.y="gene_ID",all.x=T)
# colnames(DEG_df_miR) <- c(colnames(DEG_df_miR[,1:42]),"miR-315-5p","miR-279-3p/286-3p/996-3p")
# 
# RPKM_annotate <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM_annotate_miR <- merge(rpkm_mean,attribute_papers[,c(1,7:8)],by = "flybase_gene_id",all.x=T)
colnames(RPKM_annotate_miR) <- c(colnames(RPKM_annotate_miR)[1:14],"Seed_iab4","Score_iab4","Seed_iab8","Score_iab8")
RPKM_annotate_miR <- merge(RPKM_annotate_miR,miR_315_279[,-2],by.x="flybase_gene_id",by.y="gene_ID",all.x=T) %>% dplyr::select(1:2,15:20,everything())
RPKM_annotate_miR <- RPKM_annotate_miR %>% dplyr::rename("miR-315-5p"=7,"miR-279-3p/286-3p/996-3p"=8)

No_279 <- subset(DEG_df,is.na(miR279))
No_315 <- subset(DEG_df,is.na(miR315))
miR315 <- subset(DEG_df,!is.na(miR315))
miR279 <- subset(DEG_df,!is.na(miR279))

ks_315 <- format(ks.test(miR315$logFC, No_315$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
if (ks_315 == "0e+00"){ks_315 = "< 2.2e-16"}
ks_279 <- format(ks.test(miR279$logFC, No_279$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
if (ks_279 == "0e+00"){ks_279 = "< 2.2e-16"}
  
pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/3L/ECDF/miRNA/miR315.pdf"),width = 10,height = 10)
print( ggplot(No_315, aes(logFC)) + stat_ecdf(geom = "step",size=1) + annotate("text",x=5,y=0.12,label = paste0("No Site (",nrow(No_site),")"), color="#000000",size=4,hjust = 1) +
         stat_ecdf(data = miR315,aes(logFC), geom = "step",color="#FF8744",size=1) + annotate("text",x=5,y=0.09,label = paste0("miR315 (",nrow(miR315),") pValue = ",ks_315), color="#FF8744",size=4,hjust = 1) +
         labs(x="\nlog2[mRNA Fold Change]",y="Cumulative fraction \n",title="miR-315 [Recue - KO]") + xlim(-5,5) + 
         theme(panel.background = element_rect(fill = "white", colour = "black"),
               axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
               axis.text.x = element_text(hjust = 1),
               plot.title = element_text(size=20, face="bold",hjust=0.5)))
dev.off()

pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/3L/ECDF/miRNA/miR279.pdf"),width = 10,height = 10)
print( ggplot(No_279, aes(logFC)) + stat_ecdf(geom = "step",size=1) + annotate("text",x=5,y=0.12,label = paste0("No Site (",nrow(No_279),")"), color="#000000",size=4,hjust = 1) +
         stat_ecdf(data = miR279,aes(logFC), geom = "step",color="#FF8744",size=1) + annotate("text",x=5,y=0.09,label = paste0("miR279 (",nrow(miR279),") pValue = ",ks_279), color="#FF8744",size=4,hjust = 1) +
         labs(x="\nlog2[mRNA Fold Change]",y="Cumulative fraction \n",title="miR-279 [Recue - KO]") + xlim(-5,5) + 
         theme(panel.background = element_rect(fill = "white", colour = "black"),
               axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
               axis.text.x = element_text(hjust = 1),
               plot.title = element_text(size=20, face="bold",hjust=0.5)))
dev.off()


No_279 <- subset(x,is.na(miR279))
No_315 <- subset(x,is.na(miR315))
miR315 <- subset(x,!is.na(miR315))
miR279 <- subset(x,!is.na(miR279))

ks_315 <- format(ks.test(miR315$logFC, No_315$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
if (ks_315 == "0e+00"){ks_315 = "< 2.2e-16"}
ks_279 <- format(ks.test(miR279$logFC, No_279$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
if (ks_279 == "0e+00"){ks_279 = "< 2.2e-16"}

  No_site <- subset(x,is.na(Seed_iab4)&is.na(Seed_iab8)&is.na(`miR-315-5p`)&is.na(`miR-279-3p/286-3p/996-3p`))
  iab4 <- subset(x,!is.na(Seed_iab4))
  iab8 <- subset(x,!is.na(Seed_iab8))
  m315 <- subset(x,!is.na(`miR-315-5p`)) 
  m279 <- subset(x,!is.na(`miR-279-3p/286-3p/996-3p`)) 
  ks1 <- format(ks.test(iab4$`log2[RPKM]`, No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks1 == "0e+00"){ks1 = "< 2.2e-16"}
  ks2 <- format(ks.test(iab8$`log2[RPKM]`, No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks2 == "0e+00"){ks2 = "< 2.2e-16"}
  ks3 <- format(ks.test(m315$`log2[RPKM]`,No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks3 == "0e+00"){ks3 = "< 2.2e-16"}
  ks4 <- format(ks.test(m279$`log2[RPKM]`, No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks4 == "0e+00"){ks4 = "< 2.2e-16"}
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/miRNA/RPKM/",name,".pdf"),width = 10,height = 10)
  print(
    ggplot(No_site, aes(`log2[RPKM]`)) + stat_ecdf(geom = "step",size=1) + annotate("text",x=10,y=0.12,label = paste0("No Site (",nrow(No_site),")"), color="#000000",size=4,hjust = 1) +
      stat_ecdf(data = iab4,aes(`log2[RPKM]`), geom = "step",color="#008744",size=1) + annotate("text",x=10,y=0.09,label = paste0("iab4_5p (",nrow(iab4),") pValue = ",ks1), color="#008744",size=4,hjust = 1) +
      stat_ecdf(data = iab8,aes(`log2[RPKM]`), geom = "step",color="#ffa700",size=1) + annotate("text",x=10,y=0.06,label = paste0("iab8_5p (",nrow(iab8),") pValue = ",ks2), color="#ffa700",size=4,hjust = 1) +
      stat_ecdf(data = m315,aes(`log2[RPKM]`), geom = "step",color="#0057e7",size=1) + annotate("text",x=10,y=0.03,label = paste0("miR315_5p (",nrow(m315),") pValue = ",ks3), color="#0058e7",size=4,hjust = 1) +
      stat_ecdf(data = m279,aes(`log2[RPKM]`), geom = "step",color="#d62d20",size=1) + annotate("text",x=10,y=0,label = paste0("miR279_5p (",nrow(m279),") pValue = ",ks4), color="#d62d20",size=4,hjust = 1) +
      labs(x="\nlog2[Expression]",y="Cumulative fraction \n",title=name) + xlim(-0.5,10) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5)))
  
  dev.off()
}

rm(list = c("iab4","iab8","DEG_df","DEG_df_miR","RPKM_annotate","RPKM_annotate_miR","cutoff","name","condition","x","No_site","m315","m279",ls()[grep("ks",ls())]))