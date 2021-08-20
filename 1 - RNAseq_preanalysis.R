setwd("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/")
#ann <- "E:/Drosophila/ucsc_dme_r6.21/Drosophila_melanogaster.BDGP6.94.chr.gff3"

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

Predicted_iab_Targets_score <- read_delim("E:/TargetScan/iab/Predicted_Targets_Context_Scores.default_predictions_iab.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_iab_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")
Predicted_iab_Targets <- Predicted_iab_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(1:3,5:6,11) %>% 
  dplyr::mutate(Site_Type=ifelse(Site_Type==1,"7mer-a1",ifelse(Site_Type==2,"7mer-m8","8mer")),miRNA = gsub("dme-miR-","",miRNA)) %>% dplyr::filter(miRNA != "iab-4-3p") %>% 
  group_by(gene_ID, Gene_Symbol, Transcript_ID, miRNA) %>% dplyr::slice(which.min(weighted_context_score)) %>% 
  tidyr::separate(miRNA,c("miR","N","P"),sep="-") %>% tidyr::unite("miRNA",c("miR","N"),sep="") %>% tidyr::unite("site",c("P","Site_Type"),sep=" - ") %>% tidyr::unite("miRNA",c("miRNA","site")) %>% 
  tidyr::spread(key=miRNA,value=weighted_context_score)

Predicted_315_Targets <- read_delim("E:/TargetScan/iab/Predicted_Targets_Info.default_predictions_315.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_315_Targets) <- c("miR","gene_ID","gene_name","Transcript_ID","Species_ID","UTR_start","UTR_end","MSA_start","MSA_end","Seed_match","PCT")
Predicted_315_Targets <- Predicted_315_Targets %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(2,3,1,10) %>% group_by(gene_ID,miR,gene_name) %>% 
  dplyr::summarise(Seed_match=paste(Seed_match,collapse=",")) %>% tidyr::spread(key=miR,value=Seed_match)

Predicted_279_Targets <- read_delim("E:/TargetScan/iab/Predicted_Targets_Info.default_predictions_279.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_279_Targets) <- c("miR","gene_ID","gene_name","Transcript_ID","Species_ID","UTR_start","UTR_end","MSA_start","MSA_end","Seed_match","PCT")
Predicted_279_Targets <- Predicted_279_Targets %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(2,3,1,10) %>% group_by(gene_ID,miR,gene_name) %>% 
  dplyr::summarise(Seed_match=paste(Seed_match,collapse=",")) %>% tidyr::spread(key=miR,value=Seed_match)

Slattery_2013 <- read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pgen.1003753.s001.XLSX",
                            sheet = "Hth ChIP Peaks", col_types = c("blank", "blank", "blank", "text", "blank", "blank", "blank", "blank", "blank", 
                                                                    "text", "blank", "blank", "blank", "blank", "blank", "text"), skip = 1)

colnames(Slattery_2013) <- c("Hth EA=W","Hth W>EA","Hth EA>W")
xx <- data.frame("gene_ID"=NULL,"Slattery_2013"=NULL)
for (col in (1:3)){
  x <- Slattery_2013[,col] %>% dplyr::mutate(Slattery_2013 = names(Slattery_2013)[col]) %>% dplyr::select(gene_ID=1,2)
  xx <- rbind(xx,x)
}

xx$Slattery_2013 <- factor(xx$Slattery_2013, levels = c("Hth EA=W","Hth EA>W","Hth W>EA"))
Slattery_2013 <- unique(unique(xx) %>% dplyr::group_by(gene_ID) %>% dplyr::summarise(Slattery_2013=dplyr::first(Slattery_2013)) %>% tidyr::separate_rows(gene_ID,sep=";"))

Slattery_2011_HC <- unique(read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s007.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
                             dplyr::select(gene_ID = 1, gene_name = 2) %>% tidyr::separate_rows(gene_ID,sep=";") %>% tidyr::separate_rows(gene_name,sep=";") %>% 
                             dplyr::mutate(gene_name=gsub('.{3}$', '', gene_name),Slattery_2011="Hth High-Confidence"))
Slattery_2011_MC <- unique(read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s009.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
                             dplyr::select(gene_ID = 1, gene_name = 2) %>% tidyr::separate_rows(gene_ID,sep=";") %>% tidyr::separate_rows(gene_name,sep=";") %>%
                             dplyr::mutate(gene_name=gsub('.{3}$', '', gene_name),Slattery_2011="Hth Medium-Confidence"))
Slattery_2011 <- merge(Slattery_2011_HC,Slattery_2011_MC,by="gene_ID", all= T) %>% dplyr::mutate(Slattery_2011.y=ifelse(is.na(Slattery_2011.x),Slattery_2011.y,Slattery_2011.x)) %>%
  dplyr::select(1,"gene_name"=4,"Slattery_2011"=5)

conversion <- read.delim("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/unknown_conversion.txt", header=T)

Slattery_2011 <- merge(conversion,Slattery_2011,by.x="Submitted_id",by.y="gene_ID",all.y=T) %>% select(gene_ID=2,gene_name=4,6)
Slattery_2013 <- unique(merge(conversion,Slattery_2013,by.x="Submitted_id",by.y="gene_ID",all.y=T) %>% select(gene_ID=2,gene_name=4,5)) %>% dplyr::group_by(gene_ID) %>%
  dplyr::summarise(gene_name=dplyr::first(gene_name),Slattery_2013=dplyr::first(Slattery_2013))

miR_315_279 <- merge(Predicted_315_Targets,Predicted_279_Targets[,-2],by="gene_ID",all=T)
attribute_papers_iab <- merge(unique(tx4fly[,1:2]),merge(merge(Slattery_2011[,c(1,3)],Slattery_2013[,c(1,3)],by="gene_ID",all=T),Predicted_iab_Targets[,c(1,4:9)],by="gene_ID",all=T),
                              by.x="flybase_gene_id",by.y="gene_ID",all.y=T) %>% dplyr::group_by(flybase_gene_id) 
# attribute_papers <- merge(attribute_papers,miR_315_279[,-2],by.x="flybase_gene_id",by.y="gene_ID",all=T)

# write.csv(attribute_papers_iab,"L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/attribute_paper_iab.csv",row.names = F)
# write.csv(miR_315_279,"L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/target_miR_315_279.csv",row.names = F)

rm(list=c("mart","x","xx","col","names","Slattery_2011_HC","Slattery_2011_MC","conversion","mart","Predicted_279_Targets","Predicted_315_Targets","Slattery_2013","Slattery_2011","Predicted_iab_Targets"))

########## dataframe with reads number for experiment

load("E:/Project_09228_B/PITT_0298/R_analysis/feature_count_DANI.RData")

RNA_count <- as.data.frame(SEQrna_tot$counts)

names <- NULL 
for (i in seq(1:18)){
  names[i] <- strsplit(names(RNA_count),"[.]")[[i]][6]} 

colnames(RNA_count) <- names

RNA_count <- tx4fly %>% dplyr::group_by(flybase_gene_id) %>% 
  dplyr::summarise(external_gene_name=dplyr::first(external_gene_name),entrezgene=dplyr::first(entrezgene),transcript_length=max(transcript_length)) %>% 
  merge(.,RNA_count,by.x="flybase_gene_id",by.y="row.names",all.y=T)

# write.csv(RNA_count,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RNA_count.csv",row.names = F)

rm(list=c("SEQrna_tot","SEQrna_UTR_tot","i","names"))

########## RPKM

RNA_count <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RNA_count.csv",header = T)

rpkm <- cbind(RNA_count[1:4],rpkm(RNA_count[,5:22],gene.length=RNA_count$transcript_length))
RNA_rpkm <- rpkm %>% tidyr::gather(key="Experiment",value = "rpkm",-flybase_gene_id,-external_gene_name,-entrezgene,-transcript_length) %>% dplyr::mutate(group = gsub('.{3}$', '', Experiment))

rpkm_annotate <- unique(merge(rpkm,attribute_papers_iab[,-2],by="flybase_gene_id",all.x=T) %>% dplyr::select(1:2,4,23:30,dplyr::everything()))
rpkm_annotate[is.na(rpkm_annotate)] <- ""

RNA_group <- factor(RNA_rpkm$group)

rpkm_mean <- rpkm[,1:2]
for (exp in levels(RNA_group)){
  x <- rpkm[,grep(exp,names(rpkm))] %>% dplyr::mutate(!!paste0(exp,"_mean"):=rowMeans(.[,1:3]),!!paste0(exp,"_SD"):=rowSds(as.matrix(.[,1:3])))
  rpkm_mean = cbind(rpkm_mean,x[,4:5])
}
rpkm_mean <- rpkm_mean %>% dplyr::select(1,2,dplyr::ends_with("mean"),dplyr::ends_with("SD")) %>% merge(.,attribute_papers_iab[,-2],by="flybase_gene_id",all.x=T) %>% dplyr::select(1:2,15:22,dplyr::everything())
rpkm_mean[is.na(rpkm_mean)] <- ""

# write.csv(rpkm_annotate,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_with_miR_hth.csv",row.names = F)
# write.csv(rpkm_mean,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv",row.names = F)


########## Preliminary Graphs

palet_boxplot <- c(rep("#fff866",3),rep("#fff400",3),rep("#ff467e",3),rep("#f01f1f",3),rep("#467eff",3),rep("#0b25ee",3))

pdf("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM//Density_plot_RPKM_total.pdf",width = 10,height = 8)
ggplot(RNA_rpkm, aes(x=log2(rpkm),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(palet_boxplot)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
ggplot(RNA_rpkm, aes(x=log2(rpkm+1),fill=group)) + geom_density(alpha=0.1) +
  scale_fill_manual(values = unique(palet_boxplot)) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
dev.off()

pdf("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Density_plot_RPKM_for_Experiment.pdf",width = 21,height = 8)
RNA_rpkm %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = palet_boxplot) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
RNA_rpkm %>% tidyr::separate(Experiment,into=c("x","y","replica"),sep='_',remove=F) %>%
  ggplot(., aes(x=log2(rpkm+1),fill=Experiment)) + geom_density(alpha=0.5) +
  facet_grid(replica~group) +
  scale_fill_manual(values = palet_boxplot) +
  labs(caption = "\n Density plot for RPKM", y="Density \n",x="\nlog2[RPKM+1]",fill = "Genotype") + 
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=24), axis.title = element_text(size=24,face = "bold"))
dev.off()

pdf("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Boxplot_RPKM_single_experiment.pdf",width = 10,height = 8)
ggplot(RNA_rpkm, aes(x=Experiment,y=log2(rpkm), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = palet_boxplot) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(RNA_rpkm, aes(x=Experiment,y=log2(rpkm+1), fill=Experiment)) + geom_boxplot() +
  labs(caption = "\n Single genotype RPKM distribution", y="log2[RPKM+1]\n",x="\nGenotype",fill = "Genotype") +
  scale_fill_manual(values = palet_boxplot) +
  theme(aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# rpkm_mean %>% dplyr::select(dplyr::contains("mean")) %>% tidyr::gather(key="Experiment",value="RPKM") %>%
#   ggplot(.,aes(x=Experiment,y=log2(RPKM),fill=Experiment)) + geom_boxplot()

### Scatter and R^2
data <- rpkm_mean %>% select(1,2,dplyr::ends_with("mean"))

mark <- data %>% dplyr::filter(external_gene_name %in% to_check)

df <- data.frame("X"=NULL,"Y"=NULL,"R"=NULL)
df1 <- data.frame(matrix(NA, nrow = length(data)-2, ncol = length(data)-2),row.names = names(data[,-c(1:2)]))
colnames(df1) <- rownames(df1)

gene_color <- c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")

x = 3
while (x <= ncol(data)){
  y = x+1
  while (y <= ncol(data)){
    x_name = names(data)[x]
    x_genotype = gsub("_mean","",names(data)[x])
    y_name = names(data)[y]
    y_genotype = gsub("_mean","",names(data)[y])
    col = "gene_id"
    R = data.frame("Label"=paste0("R^2 = ",round(summary(lm(data[,y] ~ data[,x]))$adj.r.squared,2)))
    graph <- data %>% dplyr::select(one_of("flybase_gene_id","external_gene_name",x_name,y_name)) %>% dplyr::select(gene_id=1,2,C1=3,C2=4)
    graph_mark <- mark %>% dplyr::select(one_of("flybase_gene_id","external_gene_name",x_name,y_name)) %>% dplyr::select(gene_id=1,2,C1=3,C2=4)
    sp <-   ggplot(graph,aes(x=log2(C1),y=log2(C2),label=external_gene_name)) + geom_point(aes(alpha=0.5),size=1) +
      geom_point(data=graph_mark,aes(x=log2(C1),y=log2(C2), color=external_gene_name), alpha=1, size=3) +
      scale_color_manual(values = gene_color) + xlim(-10,15) + ylim(-10,15) +
      scale_alpha(guide = 'none') + scale_size(guide = 'none') +
      labs(x=paste0("log2[",x_genotype,"]"),y=paste0("log2[",y_genotype,"]"), color= "Genes") +
      geom_abline(slope = 1, intercept = 0, color = "#990000", size=0.5) + 
      geom_hline(yintercept = 0, colour = "gray65",size = 0.5, linetype = "dotted") + 
      geom_vline(xintercept = 0, colour = "gray65",size = 0.5, linetype = "dotted") +
      geom_text( data = R,  mapping = aes(x = 10, y = -7.5, label=Label, alpha = 1), hjust   = -0.1, vjust   = -1) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5),
            legend.position="bottom") 
    spp <- ggplotly(sp,tooltip = "label")
    pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Scatter/",x_genotype,"-",y_genotype,".pdf"),height = 8,width = 10,useDingbats=FALSE)
    print(sp)
    dev.off()
    htmlwidgets::saveWidget(spp, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Scatter/",x_genotype,"-",y_genotype,".html"))
    df <- rbind(df,data.frame("X"=x_name,"Y"=y_name,"R"=round(summary(lm(data[,y] ~ data[,x]))$adj.r.squared,5)))
    df1[x_name,y_name] <- round(summary(lm(data[,y] ~ data[,x]))$adj.r.squared,5)
    y = y + 1
  }
  x = x+1
}


write.csv(df1,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Scatter/R_matrix.csv")
write.csv(df,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Scatter/R.csv")

rm(list = c("data","mark","df","df1","x","y","x_genotype","y_genotype","x_name","y_name","gene_color","col","R","graph","graph_mark","sp","spp"))

#########    EdgeR analysis - DEG

RNA_group <- factor(gsub('.{3}$', '', names(RNA_count[,5:ncol(RNA_count)])))

RNA_edger <- DGEList(counts = RNA_count[,5:ncol(RNA_count)],group=RNA_group,genes=RNA_count[,1:4])
RNA_edger <- calcNormFactors(RNA_edger)

design <- model.matrix(~0+RNA_group)
colnames(design) <- levels(RNA_group)

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

RNA_edger <- estimateDisp(RNA_edger, design)
RNA_fit <- glmQLFit(RNA_edger, design)

##### MMD total with 7 gene

DEG_df <- glmQLFTest(RNA_fit,contrast=my.contrasts[,1])$genes[,1:2]

DEG = list()
i=1

while(i<=dim(my.contrasts)[2]){
  comparison <- glmQLFTest(RNA_fit,contrast=my.contrasts[,i])
  condition <- my.contrasts[which(my.contrasts[,i]!=0),i]
  table <- cbind(comparison$genes[,1],comparison$table)
  colnames(table) <-c("flybase_gene_id",paste0(colnames(table[,-1]),"[",names(which(condition>0))," - ",names(which(condition<0)),"]"))
  DEG_df <- merge(DEG_df,table,by="flybase_gene_id",all=T) 
  match <- paste("flybase_gene_id",names(condition)[1],names(condition)[2],sep="|")
  mean_rpkm <- rpkm %>% dplyr::select(dplyr::matches(match)) %>% dplyr::mutate(mean = log2(rowMeans(.[,2:7])+1)) %>% dplyr::select(flybase_gene_id=1 ,`log2[RPKM+1]`=8) %>%
    dplyr::mutate(exp=paste(names(which(condition>0)),names(which(condition<0)),sep=" - "))
  DEG[[colnames(my.contrasts)[i]]] = merge(as.data.frame(cbind(comparison$genes[,1:2],comparison$table)),mean_rpkm,by="flybase_gene_id",all.x=T)
  mark = DEG[[i]] %>% dplyr::filter(external_gene_name %in% to_check)

  sc <- ggplot(DEG[[i]],aes(y=logFC,x=`log2[RPKM+1]`,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=DEG[[i]]%>% dplyr::filter(external_gene_name %in% to_check),aes(y=logFC,x=`log2[RPKM+1]`, fill=external_gene_name), alpha=1, size=3,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    labs(y = paste0("logFC[",names(which(condition>0))," - ",names(which(condition<0)),"]"),x="log2[RPKM+1]") +
    geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")

  sc1 <- sc + geom_vline(xintercept = 1, size = 0.5, color="#990000" )
  sc5 <- sc + geom_vline(xintercept = 2.585, size = 0.5, color="#009900" )

  vp <- ggplot(DEG[[i]],aes(x= exp, y=`log2[RPKM+1]`)) + geom_violin(fill="#dedede") +
    coord_flip() +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_blank(),axis.title = element_blank(),
          plot.title = element_blank())

  vp1 <- vp + geom_hline(yintercept = 1, size = 0.5, color="#990000" ) +
    annotate("text",x=1.48,y=0.3,label = paste0("< 1 RPKM\n",table(DEG[[i]]$`log2[RPKM+1]`>1)[1]), color="#990000",size=3) +
    annotate("text",x=1.48,y=1.7,label = paste0("> 1 RPKM\n",table(DEG[[i]]$`log2[RPKM+1]`>1)[2]), color="#990000",size=3)
  vp5 <- vp + geom_hline(yintercept = 2.585, size = 0.5, color="#009900" ) +
    annotate("text",x=1.48,y=1.7,label = paste0("< 5 RPKM\n",table(DEG[[i]]$`log2[RPKM+1]`>2.585)[1]), color="#009900",size=3) +
    annotate("text",x=1.48,y=3.3,label = paste0("> 5 RPKM\n",table(DEG[[i]]$`log2[RPKM+1]`>2.585)[2]), color="#009900",size=3)

  Volp <- ggplot(DEG[[i]],aes(y=-log10(PValue),x=logFC,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=DEG[[i]]%>% dplyr::filter(external_gene_name %in% to_check),aes(y=-log10(PValue),x=logFC, fill=external_gene_name), alpha=1, size=3,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    labs(y = "log10[pValue]",title = paste0(names(which(condition>0))," - ",names(which(condition<0))),x="logFC") +
    geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
    geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")

  ps1 <- ggplotly(sc1,tooltip = "label")
  ps5 <- ggplotly(sc5,tooltip = "label")
  pv <- ggplotly(Volp,tooltip = "label")

  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp1,sc1,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps1, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/",colnames(my.contrasts)[i],".html"))

  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp5,sc5,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps5, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/",colnames(my.contrasts)[i],".html"))


  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(Volp)
  dev.off()
  htmlwidgets::saveWidget(pv, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/",colnames(my.contrasts)[i],".html"))
  i=i+1
}

write.csv(DEG_df,"L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",row.names = F)
rm(list = c("comparison","condition","table","DEG_df","match","mean_rpkm","DEG","mark","sc","sc1","sc5","vp","vp1","vp5","Volp","ps1","ps5","pv"))

##### MMD total for miRNA

iab4 <- attribute_papers_iab %>% dplyr::select(flybase_gene_id,external_gene_name,dplyr::starts_with("iab4")) %>% dplyr::rename("7mer.a1"=3, "7mer.m8"=4, "8mer"=5) %>%
  dplyr::filter(!is.na(`7mer.a1`) | !is.na(`7mer.m8`) | !is.na(`8mer`)) %>% tidyr::gather(key="Seed",value="Score",-flybase_gene_id, -external_gene_name) %>% dplyr::filter(!is.na(Score))
iab4$Seed <- factor(iab4$Seed, levels = c("8mer","7mer.m8","7mer.a1"))
iab4$flybase_gene_id <- as.character(iab4$flybase_gene_id)
iab4$external_gene_name <- as.character(iab4$external_gene_name)

iab8 <- attribute_papers_iab %>% dplyr::select(flybase_gene_id,external_gene_name,dplyr::starts_with("iab8")) %>% dplyr::rename("7mer.a1"=3, "7mer.m8"=4, "8mer"=5) %>%
  dplyr::filter(!is.na(`7mer.a1`) | !is.na(`7mer.m8`) | !is.na(`8mer`)) %>% tidyr::gather(key="Seed",value="Score",-flybase_gene_id, -external_gene_name) %>% dplyr::filter(!is.na(Score))
iab8$Seed <- factor(iab8$Seed, levels = c("8mer","7mer.m8","7mer.a1"))
iab8$flybase_gene_id <- as.character(iab8$flybase_gene_id)
iab8$external_gene_name <- as.character(iab8$external_gene_name)

DEG_df <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_df_miR <- merge(merge(DEG_df,iab4[,-2],by = "flybase_gene_id",all.x=T),iab8[,-2],by = "flybase_gene_id",all.x=T)
colnames(DEG_df_miR) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_df_miR[,-c(1:2,39:42)])))),"Seed_iab4","Score_iab4","Seed_iab8","Score_iab8")

i=1

while(i<=dim(my.contrasts)[2]){
  condition <- my.contrasts[which(my.contrasts[,i]!=0),i]
  name <- paste0("[",names(which(condition>0))," - ",names(which(condition<0)),"]")
  table <- DEG_df_miR %>% select(flybase_gene_id,external_gene_name,dplyr::contains(name),Seed_iab4, Score_iab4, Seed_iab8, Score_iab8)
  colnames(table) <- c("flybase_gene_id","external_gene_name","logFC","logCPM","F","PValue","Seed_iab4","Score_iab4","Seed_iab8","Score_iab8")
  
  match <- paste("flybase_gene_id",names(condition)[1],names(condition)[2],sep="|")
  mean_rpkm <- rpkm %>% dplyr::select(dplyr::matches(match)) %>% dplyr::mutate(mean = log2(rowMeans(.[,2:7])+1)) %>% dplyr::select(flybase_gene_id=1 ,`log2[RPKM+1]`=8) %>%
    dplyr::mutate(exp=paste(names(which(condition>0)),names(which(condition<0)),sep=" - "))
  table = merge(table,mean_rpkm,by="flybase_gene_id",all.x=T)
  markD4 = table %>% dplyr::filter(is.na(Seed_iab4))
  markD8 = table %>% dplyr::filter(is.na(Seed_iab8))
  mark4 = table %>% dplyr::filter(!is.na(Seed_iab4))
  mark8 = table %>% dplyr::filter(!is.na(Seed_iab8))
  
  pt4 <- format(prop.test(x=c(table(mark4$logFC>0)[2],table(markD4$logFC>0)[2]),n=c(nrow(mark4),nrow(markD4)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt4_8mer <- format(prop.test(x=c(table((mark4%>%dplyr::filter(Seed_iab4=="8mer"))$logFC>0)[2],table(markD4$logFC>0)[2]),
                               n=c(nrow(mark4%>%dplyr::filter(Seed_iab4=="8mer")),nrow(markD4)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt4_7a1 <- format(prop.test(x=c(table((mark4%>%dplyr::filter(Seed_iab4=="7mer.a1"))$logFC>0)[2],table(markD4$logFC>0)[2]),
                               n=c(nrow(mark4%>%dplyr::filter(Seed_iab4=="7mer.a1")),nrow(markD4)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt4_7m8 <- format(prop.test(x=c(table((mark4%>%dplyr::filter(Seed_iab4=="7mer.m8"))$logFC>0)[2],table(markD4$logFC>0)[2]),
                               n=c(nrow(mark4%>%dplyr::filter(Seed_iab4=="7mer.m8")),nrow(markD4)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  
  pt8 <- format(prop.test(x=c(table(mark8$logFC>0)[2],table(markD8$logFC>0)[2]),n=c(nrow(mark8),nrow(markD8)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt8_8mer <- format(prop.test(x=c(table((mark8%>%dplyr::filter(Seed_iab8=="8mer"))$logFC>0)[2],table(markD8$logFC>0)[2]),
                               n=c(nrow(mark8%>%dplyr::filter(Seed_iab8=="8mer")),nrow(markD8)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt8_7a1 <- format(prop.test(x=c(table((mark8%>%dplyr::filter(Seed_iab8=="7mer.a1"))$logFC>0)[2],table(markD8$logFC>0)[2]),
                              n=c(nrow(mark8%>%dplyr::filter(Seed_iab8=="7mer.a1")),nrow(markD8)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt8_7m8 <- format(prop.test(x=c(table((mark8%>%dplyr::filter(Seed_iab8=="7mer.m8"))$logFC>0)[2],table(markD8$logFC>0)[2]),
                              n=c(nrow(mark8%>%dplyr::filter(Seed_iab8=="7mer.m8")),nrow(markD8)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  
  sc_4 <- ggplot(table,aes(y=logFC,x=`log2[RPKM+1]`,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=mark4 ,aes(y=logFC,x=`log2[RPKM+1]`, fill=Seed_iab4, alpha=Score_iab4*-1), size=2,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    scale_alpha(guide="none") +
    labs(y = paste0("logFC[",names(which(condition>0))," - ",names(which(condition<0)),"]"),x="log2[RPKM+1]") +
    geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")
  sc_8 <- ggplot(table,aes(y=logFC,x=`log2[RPKM+1]`,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=mark8 ,aes(y=logFC,x=`log2[RPKM+1]`, fill=Seed_iab8, alpha=Score_iab8*-1), size=2,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    scale_alpha(guide="none") +
    labs(y = paste0("logFC[",names(which(condition>0))," - ",names(which(condition<0)),"]"),x="log2[RPKM+1]") +
    geom_hline(yintercept = c(-0.59,0.59), color="#999999", size = 0.5, linetype="dotdash") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")
  sc14 <- sc_4 + geom_vline(xintercept = 1, size = 0.5, color="#990000" )
  sc54 <- sc_4 + geom_vline(xintercept = 2.585, size = 0.5, color="#009900" )
  sc18 <- sc_8 + geom_vline(xintercept = 1, size = 0.5, color="#990000" )
  sc58 <- sc_8 + geom_vline(xintercept = 2.585, size = 0.5, color="#009900" )
  
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
  
  Volp4 <- ggplot(table,aes(y=-log10(PValue),x=logFC,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=mark4,aes(y=-log10(PValue),x=logFC, fill=Seed_iab4, alpha=Score_iab4*-1), size=3,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    scale_alpha(guide="none") +
    labs(y = "log10[pValue]",title = paste0(names(which(condition>0))," - ",names(which(condition<0))),x="logFC") +
    annotate("text",x=-7,y=17,label=paste0("Total iab4 target pValue = ",pt4),hjust=0.1) +
    annotate("text",x=-7,y=16.5,label=paste0("mer7-a1 iab4 target pValue = ",pt4_7a1),color="#ffe119",hjust=0.1) +
    annotate("text",x=-7,y=16,label=paste0("mer7-m8 iab4 target pValue = ",pt4_7m8),color="#3cd44b",hjust=0.1) +
    annotate("text",x=-7,y=15.5,label=paste0("mer8 iab4 target pValue = ",pt4_8mer),color="#e6194b",hjust=0.1) +
    geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
    geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")
  Volp8 <- ggplot(table,aes(y=-log10(PValue),x=logFC,label=external_gene_name)) + geom_point(alpha=0.3,color="#dedede") +
    geom_point(data=mark8,aes(y=-log10(PValue),x=logFC, fill=Seed_iab8, alpha=Score_iab8*-1), size=3,shape=21) +
    scale_fill_manual(values = c("#e6194b","#3cd44b","#ffe119","#4363d8","#f58231","#e6beff","#000075")) +
    scale_alpha(guide="none") +
    labs(y = "log10[pValue]",title = paste0(names(which(condition>0))," - ",names(which(condition<0))),x="logFC") +
    annotate("text",x=-7,y=17,label=paste0("Total iab8 target pValue = ",pt8),hjust=0.1) +
    annotate("text",x=-7,y=16.5,label=paste0("mer7-a1 iab8 target pValue = ",pt8_7a1),color="#ffe119",hjust=0.1) +
    annotate("text",x=-7,y=16,label=paste0("mer7-m8 iab8 target pValue = ",pt8_7m8),color="#3cd44b",hjust=0.1) +
    annotate("text",x=-7,y=15.5,label=paste0("mer8 iab8 target pValue = ",pt8_8mer),color="#e6194b",hjust=0.1) +
    geom_hline(yintercept = 1.30103, color="#990000", size = 0.5) +
    geom_vline(xintercept = c(-0.5849625,0.5849625), size = 0.5, color="#999999", linetype="dotdash" ) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"),
          axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size=20, face="bold",hjust=0.5),
          legend.position="bottom")
  
  ps14 <- ggplotly(sc14,tooltip = "label")
  ps18 <- ggplotly(sc18,tooltip = "label")
  ps54 <- ggplotly(sc54,tooltip = "label")
  ps58 <- ggplotly(sc58,tooltip = "label")
  pv4 <- ggplotly(Volp4,tooltip = "label")
  pv8 <- ggplotly(Volp8,tooltip = "label")
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/miRNA/iab4/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp1,sc14,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps14, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/miRNA/iab4/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/miRNA/iab8/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp1,sc18,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps18, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/miRNA/iab8/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/miRNA/iab4/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp5,sc54,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps54, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/miRNA/iab4/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/miRNA/iab8/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp5,sc58,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps58, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/miRNA/iab8/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/miRNA/iab4/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(Volp4)
  dev.off()
  htmlwidgets::saveWidget(pv4, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/miRNA/iab4/",colnames(my.contrasts)[i],".html"))

  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/miRNA/iab8/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(Volp8)
  dev.off()
  htmlwidgets::saveWidget(pv8, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/miRNA/iab8/",colnames(my.contrasts)[i],".html"))
  
  i=i+1
}

rm(list = c("i","iab4","iab8","DEG_df_miR","condition","name","table","DEG_df","match","mean_rpkm","mark4","mark8","markD4","markD8",ls()[grep("sc",ls())],ls()[grep("vp",ls())],
            ls()[grep("Volp",ls())],ls()[grep("ps",ls())],ls()[grep("pv",ls())],ls()[grep("pt",ls())]))

########## Hth percentile
ecdf_fun <- function(x,perc) ecdf(x)(perc)

Predicted_Targets_score <- read_delim("E:/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")

iab4 <- unlist(unique((Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::filter(miRNA == "dme-miR-iab-4-5p"))$gene_ID))
iab8 <- unlist(unique((Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::filter(miRNA == "dme-miR-iab-8-5p"))$gene_ID))

DEG_df <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
colnames(DEG_df) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_df[,-c(1:2)])))))

DEG_percentile <- DEG_df %>% tidyr::gather(key="Comparison","Value",-flybase_gene_id,-external_gene_name) %>% tidyr::separate(Comparison,c("Stat","Comparison"),sep="\\[") %>% 
  dplyr::mutate(Comparison = gsub("]","",Comparison)) %>% tidyr::spread(key=Stat,value=Value) %>% dplyr::filter(PValue < 0.05,logFC > 0)

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
    
    percent <- round(ecdf_fun(subset(DEG_per_percentile,Comparison == comp)$logFC,Hth_point) * 100,2)
    
    graph <- if(length(percent) != 0 ) {
      ggplot(xd, aes(x, y)) + 
        geom_area(data = subset(xd, x >= Hth_x), fill = "#baffc9") + geom_area(data = subset(xd, x < Hth_x), fill = "#c0c5ce") +
        geom_line(size = 0.5 ) + labs(x="logFC",y="density\n",title="") +
        annotate("segment", x = Hth_x + lx, xend = Hth_x, y =Hth_y + ly , yend = Hth_y, colour = "#4363d8", size=1.5, alpha=1, arrow=arrow(angle=30,length = unit(0.1, "inches"), ends = "last", type = "closed")) +
        annotate("text", x = Hth_x + ax, y =Hth_y + ay,label = paste0("Hth (percentile =",percent,"%)"),size=5,hjust=0) +
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
    
    if(percentile=="DEG_percentile"){
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Percentile/Hth_total/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(graph)
      dev.off()
    } else if(percentile=="DEG_iab4_percentile") {
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Percentile/Hth_iab4/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(graph)
      dev.off()
    } else {
      pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Percentile/Hth_iab8/",comp,".pdf"),width = 10,height = 10,useDingbats=FALSE)
      print(graph)
      dev.off()
    }
  }
}


rm(list = c("ecdf_fun","DEG_df","DEG_per_percentile","xd","Hth_point","Hth_x","Hth_y","percent","comp"))

##### MMD total for Hth


Slattery_2013 <- unique(unique(attribute_papers_iab[,c(1:2,4)] %>% tidyr::separate_rows(Slattery_2013,sep=",")) %>% dplyr::filter(!is.na(Slattery_2013)) %>% 
                          dplyr::mutate(Interaction=ifelse(Slattery_2013=="Hth EA>W",1,ifelse(Slattery_2013=="Hth EA=W",2,3))) %>%
                          dplyr::group_by(flybase_gene_id) %>% dplyr::summarise(external_gene_name=first(external_gene_name), Interaction=first(Interaction)) %>%
                          dplyr::mutate(Interaction=ifelse(Interaction==1,"Hth EA>W",ifelse(Interaction==2,"Hth EA=W","Hth W>EA"))) %>%
                          dplyr::filter(!is.na(external_gene_name)))

Slattery_2011 <- unique(unique(attribute_papers_iab[,c(1:2,3)] %>% tidyr::separate_rows(Slattery_2011,sep=",")) %>% dplyr::filter(!is.na(Slattery_2011)) %>% 
                          dplyr::mutate(Interaction=ifelse(Slattery_2011=="Hth High-Confidence",1,2)) %>%
                          dplyr::group_by(flybase_gene_id) %>% dplyr::summarise(external_gene_name=first(external_gene_name), Interaction=first(Interaction)) %>%
                          dplyr::mutate(Interaction=ifelse(Interaction==1,"Hth High-Confidence","Hth Medium-Confidence")) %>%
                          dplyr::filter(!is.na(external_gene_name)))

DEG_df <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_df_miR <- merge(merge(DEG_df,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)
colnames(DEG_df_miR) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_df_miR[,-c(1:2,39:42)])))),"Slattery_2011","Slattery_2013")

i=1

while(i<=dim(my.contrasts)[2]){
  condition <- my.contrasts[which(my.contrasts[,i]!=0),i]
  name <- paste0("[",names(which(condition>0))," - ",names(which(condition<0)),"]")
  table <- DEG_df_miR %>% select(flybase_gene_id,external_gene_name,dplyr::contains(name),Slattery_2011,Slattery_2013)
  colnames(table) <- c("flybase_gene_id","external_gene_name","logFC","logCPM","F","PValue","Slattery_2011","Slattery_2013")
  
  match <- paste("flybase_gene_id",names(condition)[1],names(condition)[2],sep="|")
  mean_rpkm <- rpkm %>% dplyr::select(dplyr::matches(match)) %>% dplyr::mutate(mean = log2(rowMeans(.[,2:7])+1)) %>% dplyr::select(flybase_gene_id=1 ,`log2[RPKM+1]`=8) %>%
    dplyr::mutate(exp=paste(names(which(condition>0)),names(which(condition<0)),sep=" - "))
  table = merge(table,mean_rpkm,by="flybase_gene_id",all.x=T)
  markD11 = table %>% dplyr::filter(is.na(Slattery_2011))
  markD13 = table %>% dplyr::filter(is.na(Slattery_2013))
  mark11 = table %>% dplyr::filter(!is.na(Slattery_2011))
  mark13 = table %>% dplyr::filter(!is.na(Slattery_2013))
  
  pt11 <- format(prop.test(x=c(table(mark11$logFC>0)[2],table(markD11$logFC>0)[2]),n=c(nrow(mark11),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt11_HC <- format(prop.test(x=c(table((mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence"))$logFC>0)[2],table(markD11$logFC>0)[2]),
                               n=c(nrow(mark11%>%dplyr::filter(Slattery_2011=="Hth High-Confidence")),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt11_MC <- format(prop.test(x=c(table((mark11%>%dplyr::filter(Slattery_2011=="Hth Medium-Confidence"))$logFC>0)[2],table(markD11$logFC>0)[2]),
                              n=c(nrow(mark11%>%dplyr::filter(Slattery_2011=="Hth Medium-Confidence")),nrow(markD11)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  
  pt13 <- format(prop.test(x=c(table(mark13$logFC>0)[2],table(markD13$logFC>0)[2]),n=c(nrow(mark13),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt13_EA <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W"))$logFC>0)[2],table(markD13$logFC>0)[2]),
                               n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA>W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt13_EW <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W"))$logFC>0)[2],table(markD13$logFC>0)[2]),
                              n=c(nrow(mark13%>%dplyr::filter(Slattery_2013=="Hth EA=W")),nrow(markD13)),alternative = "greater",correct = T)$p.value,scientific=T,digits = 4)
  pt13_W <- format(prop.test(x=c(table((mark13%>%dplyr::filter(Slattery_2013=="Hth W>EA"))$logFC>0)[2],table(markD13$logFC>0)[2]),
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
    annotate("text",x=-7,y=17,label=paste0("Total Hth target [Slattery 2011] pValue = ",pt11),hjust=0.1) +
    annotate("text",x=-7,y=16.5,label=paste0("High confidence Hth target [Slattery 2011] pValue = ",pt11_HC),color="#e6194b",hjust=0.1) +
    annotate("text",x=-7,y=16,label=paste0("Low confidence Hth target [Slattery 2011] pValue = ",pt11_MC),color="#3cd44b",hjust=0.1) +
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
    annotate("text",x=-7,y=17,label=paste0("Total Hth target [Slattery 2013] pValue = ",pt13),hjust=0.1) +
    annotate("text",x=-7,y=16.5,label=paste0("Hth target [Slattery 2013 - W > EA] pValue = ",pt13_W),color="#ffe119",hjust=0.1) +
    annotate("text",x=-7,y=16,label=paste0("Hth target [Slattery 2013 - EA > W] pValue = ",pt13_EA),color="#3cd44b",hjust=0.1) +
    annotate("text",x=-7,y=15.5,label=paste0("Hth target [Slattery 2013 - EA = W] pValue = ",pt13_EW),color="#e6194b",hjust=0.1) +
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
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp1,sc111,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps111, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp1,sc113,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps113, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/1RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp5,sc511,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps511, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(cowplot::plot_grid(vp5,sc513,align = "v",axis = "l",nrow = 2,rel_heights = c(0.2,1)))
  dev.off()
  htmlwidgets::saveWidget(ps513, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/MMD/5RPKM/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/Hth/Slattery_2011/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(Volp11)
  dev.off()
  htmlwidgets::saveWidget(pv11, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/Hth/Slattery_2011/",colnames(my.contrasts)[i],".html"))
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/Hth/Slattery_2013/",colnames(my.contrasts)[i],".pdf"),width = 10,height = 10,useDingbats=FALSE)
  print(Volp13)
  dev.off()
  htmlwidgets::saveWidget(pv13, paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/Volcano/Hth/Slattery_2013/",colnames(my.contrasts)[i],".html"))
  
  i=i+1
}

rm(list = c("i","Slattery_2011","Slattery_2013","DEG_df_miR","condition","name","table","DEG_df","match","mean_rpkm","mark11","mark13","markD11","markD13",ls()[grep("sc",ls())],ls()[grep("vp",ls())],
            ls()[grep("Volp",ls())],ls()[grep("ps",ls())],ls()[grep("pv",ls())],ls()[grep("pt",ls())]))

########## Comulative plot (ECDF)

Slattery_2013 <- unique(unique(attribute_papers_iab[,c(1:2,4)] %>% tidyr::separate_rows(Slattery_2013,sep=",")) %>% dplyr::filter(!is.na(Slattery_2013)) %>% 
                          dplyr::mutate(Interaction=ifelse(Slattery_2013=="Hth EA>W",1,ifelse(Slattery_2013=="Hth EA=W",2,3))) %>%
                          dplyr::group_by(flybase_gene_id) %>% dplyr::summarise(external_gene_name=first(external_gene_name), Interaction=first(Interaction)) %>%
                          dplyr::mutate(Interaction=ifelse(Interaction==1,"Hth EA>W",ifelse(Interaction==2,"Hth EA=W","Hth W>EA"))) %>%
                          dplyr::filter(!is.na(external_gene_name)))

Slattery_2011 <- unique(unique(attribute_papers_iab[,c(1:2,3)] %>% tidyr::separate_rows(Slattery_2011,sep=",")) %>% dplyr::filter(!is.na(Slattery_2011)) %>% 
                          dplyr::mutate(Interaction=ifelse(Slattery_2011=="Hth High-Confidence",1,2)) %>%
                          dplyr::group_by(flybase_gene_id) %>% dplyr::summarise(external_gene_name=first(external_gene_name), Interaction=first(Interaction)) %>%
                          dplyr::mutate(Interaction=ifelse(Interaction==1,"Hth High-Confidence","Hth Medium-Confidence")) %>%
                          dplyr::filter(!is.na(external_gene_name)))

DEG_df <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/DEG_total.csv",header = T)
DEG_df_miR <- merge(merge(DEG_df,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)
colnames(DEG_df_miR) <- c("flybase_gene_id","external_gene_name",gsub("[.]","[",gsub("[.]$","]",gsub("[.]{3}"," - ",colnames(DEG_df_miR[,-c(1:2,39:42)])))),"Slattery_2011","Slattery_2013")
#DEG_df_miR <- merge(DEG_df_miR,miR_315_279[,-2],by.x="flybase_gene_id",by.y="gene_ID",all.x=T)
#colnames(DEG_df_miR) <- c(colnames(DEG_df_miR[,1:42]),"miR-315-5p","miR-279-3p/286-3p/996-3p")

RPKM_annotate <- read.csv("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/RPKM_per_genotype.csv", header = T)[,-c(3:10)] 
RPKM_annotate_miR <- merge(merge(RPKM_annotate,Slattery_2011[,-2],by = "flybase_gene_id",all.x=T),Slattery_2013[,-2],by = "flybase_gene_id",all.x=T)
colnames(RPKM_annotate_miR) <- c(colnames(RPKM_annotate_miR)[1:14],"Slattery_2011","Slattery_2013")
#RPKM_annotate_miR <- merge(RPKM_annotate_miR,miR_315_279[,-2],by.x="flybase_gene_id",by.y="gene_ID",all.x=T) %>% dplyr::select(1:2,15:20,everything())
#RPKM_annotate_miR <- RPKM_annotate_miR %>% dplyr::rename("miR-315-5p"=7,"miR-279-3p/286-3p/996-3p"=8)

cutoff <- 10

for (name in names(DEG_df_miR)[grep("logFC",names(DEG_df_miR))]){
  comparison <- gsub("logFC","",name)
  x <- cbind(DEG_df_miR[,c(1:2,39:40)],DEG_df_miR[,grep(comparison,names(DEG_df_miR),fixed=T)]) %>% 
    dplyr::rename("gene_ID"=1,"logFC"=5,"logCPM"=6,"F"=7,"PValue"=8) #%>% dplyr::filter(PValue<cutoff) %>% dplyr::filter(logFC<2,logFC>-2)
  
  No_site <- subset(x,is.na(Slattery_2011)&is.na(Slattery_2013))
  S11 <- subset(x,!is.na(Slattery_2011))
  S13 <- subset(x,!is.na(Slattery_2013))
  ks1 <- format(ks.test(S11$logFC, No_site$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks1 == "0e+00"){ks1 = "< 2.2e-16"}
  ks2 <- format(ks.test(S13$logFC, No_site$logFC, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks2 == "0e+00"){ks2 = "< 2.2e-16"}
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/Hth/logFC/",comparison,".pdf"),width = 10,height = 10)
  print(
    ggplot(No_site, aes(logFC)) + stat_ecdf(geom = "step",size=1) + annotate("text",x=7,y=0.12,label = paste0("No Site (",nrow(No_site),")"), color="#000000",size=4,hjust = 1) +
      stat_ecdf(data = S11,aes(logFC), geom = "step",color="#008744",size=1) + annotate("text",x=7,y=0.09,label = paste0("Slattere_2011 (",nrow(S11),") pValue = ",ks1), color="#008744",size=4,hjust = 1) +
      stat_ecdf(data = S13,aes(logFC), geom = "step",color="#ffa700",size=1) + annotate("text",x=7,y=0.06,label = paste0("Slattere_2011 (",nrow(S13),") pValue = ",ks2), color="#ffa700",size=4,hjust = 1) +
      labs(x="\nlog2[mRNA Fold Change]",y="Cumulative fraction \n",title=comparison) + xlim(-5,5) + 
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5)))
  dev.off()
}

CV_cutoff <- 100000

for (name in names(RPKM_annotate_miR)[grep("_mean",names(RPKM_annotate_miR))]){
  comparison <- gsub("_mean","",name)
  x <- cbind(RPKM_annotate_miR[,c(1:2,15:16)],RPKM_annotate_miR[,grep(comparison,names(RPKM_annotate_miR),fixed=T)]) %>% dplyr::rename("mean"=5,"SD"=6) %>% 
    dplyr::mutate(CV=SD/mean,`log2[RPKM]`=log2(mean)) #%>% dplyr::filter(CV<CV_cutoff) 
  
  No_site <- subset(x,is.na(Slattery_2011)&is.na(Slattery_2013))
  S11 <- subset(x,!is.na(Slattery_2011))
  S13 <- subset(x,!is.na(Slattery_2013))
  ks1 <- format(ks.test(S11$`log2[RPKM]`, No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks1 == "0e+00"){ks1 = "< 2.2e-16"}
  ks2 <- format(ks.test(S13$`log2[RPKM]`, No_site$`log2[RPKM]`, alternative = "two.sided",exact = T)$p.value,scientific=T,digits = 4)
  if (ks2 == "0e+00"){ks2 = "< 2.2e-16"}
  
  pdf(paste0("L:/dani/RNAseq/Project_09228_B/PDF/01_23_2019/RPKM/ECDF/miRNA/RPKM/",name,".pdf"),width = 10,height = 10)
  print(
    ggplot(No_site, aes(`log2[RPKM]`)) + stat_ecdf(geom = "step",size=1) + annotate("text",x=10,y=0.12,label = paste0("No Site (",nrow(No_site),")"), color="#000000",size=4,hjust = 1) +
      stat_ecdf(data = S11,aes(`log2[RPKM]`), geom = "step",color="#008744",size=1) + annotate("text",x=10,y=0.09,label = paste0("Slattery_2011 (",nrow(S11),") pValue = ",ks1), color="#008744",size=4,hjust = 1) +
      stat_ecdf(data = S13,aes(`log2[RPKM]`), geom = "step",color="#ffa700",size=1) + annotate("text",x=10,y=0.06,label = paste0("Slattery_2013 (",nrow(S13),") pValue = ",ks2), color="#ffa700",size=4,hjust = 1) +
      labs(x="\nlog2[Expression]",y="Cumulative fraction \n",title=name) + xlim(-0.5,10) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text = element_text(size=12), axis.title = element_text(size=20,face = "bold"), 
            axis.text.x = element_text(hjust = 1),
            plot.title = element_text(size=20, face="bold",hjust=0.5)))
  
  dev.off()
}

rm(list = c("iab4","iab8","DEG_df","DEG_df_miR","RPKM_annotate","RPKM_annotate_miR","cutoff","name","condition","x","No_site","m315","m279",ls()[grep("ks",ls())]))
