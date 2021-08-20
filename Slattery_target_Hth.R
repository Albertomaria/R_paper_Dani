##########
### Script to colect and convert gene_ID for Hth target predicted in Slattery et al. 2011 and Slattery et al. 2013
##########

setwd("E:/Project_09228_B/PITT_0298/R_analysis/")

library(tidyverse)
Predicted_iab_Targets_score <- read_delim("E:/TargetScan/iab/Predicted_Targets_Context_Scores.default_predictions_iab.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_iab_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")
Predicted_iab_Targets <- Predicted_iab_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(1:3,5:6,11) %>% 
  dplyr::mutate(Site_Type=ifelse(Site_Type==1,"7mer-a1",ifelse(Site_Type==2,"7mer-m8","8mer")),miRNA = gsub("dme-miR-","",miRNA)) %>% dplyr::filter(miRNA != "iab-4-3p") %>% 
  group_by(gene_ID, Gene_Symbol, Transcript_ID, miRNA) %>% dplyr::slice(which.min(weighted_context_score)) %>% 
  tidyr::separate(miRNA,c("miR","N","P"),sep="-") %>% tidyr::unite("miRNA",c("miR","N"),sep="") %>% tidyr::unite("site",c("P","Site_Type"),sep=" - ") %>% tidyr::unite("miRNA",c("miRNA","site")) %>% 
  tidyr::spread(key=miRNA,value=weighted_context_score)

mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
#View(listAttributes(mart))
tx4fly <- getBM(attributes = c("external_gene_name","flybase_gene_id","transcript_length","entrezgene"),mart=mart)

to_check=c("Ubx","abd-A","Abd-B","hth","exd","Tdc2","Ilp7")

Predicted_Targets_score <- read_delim("E:/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(Predicted_Targets_score) <- c("gene_ID","Gene_Symbol","Transcript_ID","Species_ID","miRNA","Site_Type","UTR_start","UTR_end","context_score","context_score_percentile","weighted_context_score","weighted_context_score_percentile")
Predicted_Targets_score <- Predicted_Targets_score %>% dplyr::filter(Species_ID==7227) %>% dplyr::select(1:3,5:6,11) %>% 
  dplyr::filter(miRNA %in% c("dme-miR-279-3p", "dme-miR-315-5p","dme-miR-iab-4-3p","dme-miR-iab-4-5p","dme-miR-iab-8-5p")) %>%
  dplyr::mutate(miRNA = gsub("iab-4","iab4",miRNA),miRNA = gsub("iab-8","iab8",miRNA),miRNA = gsub("dme-miR-","",miRNA),Site_Type=ifelse(Site_Type==1,"7mer-a1",ifelse(Site_Type==2,"7mer-m8","8mer"))) %>% 
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


Slattery_2013 <- read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pgen.1003753.s001.XLSX",
                            sheet = "Hth ChIP Peaks", col_types = c("blank", "blank", "blank", "text", "blank", "blank", "blank", "blank", "blank", 
                                                                    "text", "blank", "blank", "blank", "blank", "blank", "text"), skip = 1)

colnames(Slattery_2013) <- c("Hth EA=W","Hth W>EA","Hth EA>W")
xx <- data.frame("gene_ID"=NULL,"Slattery_2013"=NULL)
for (col in (1:3)){
  x <- Slattery_2013[,col] %>% dplyr::mutate(Slattery_2013 = names(Slattery_2013)[col]) %>% dplyr::select(gene_ID=1,2)
  xx <- rbind(xx,x)
}
Slattery_2013 <- unique(unique(xx) %>% dplyr::group_by(gene_ID) %>% dplyr::summarise(Slattery_2013=paste(Slattery_2013,collapse=",")) %>% tidyr::separate_rows(gene_ID,sep=";"))

Slattery_2011_HC <- unique(read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s007.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
                             dplyr::select(gene_ID = 1, gene_name = 2) %>% tidyr::separate_rows(gene_ID,sep=";") %>% tidyr::separate_rows(gene_name,sep=";") %>% 
                             dplyr::mutate(gene_name=gsub('.{3}$', '', gene_name),Slattery_2011="Hth High-Confidence"))
Slattery_2011_MC <- unique(read_excel("L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/journal.pone.0014686.s009.XLS", col_types = c("blank", "blank", "blank", "text", "text"), skip = 1) %>%
                             dplyr::select(gene_ID = 1, gene_name = 2) %>% tidyr::separate_rows(gene_ID,sep=";") %>% tidyr::separate_rows(gene_name,sep=";") %>%
                             dplyr::mutate(gene_name=gsub('.{3}$', '', gene_name),Slattery_2011="Hth Medium-Confidence"))
Slattery_2011 <- rbind(Slattery_2011_HC,Slattery_2011_MC) %>% dplyr::group_by(gene_ID) %>% dplyr::summarise(gene_name=first(gene_name),Slattery_2011=paste(Slattery_2011,collapse = ","))

write.csv(data.frame("gene_ID"=unique(c(Slattery_2013$gene_ID,Slattery_2011$gene_ID))),"L:/dani/RNAseq/Project_09228_B/HTH and Ubx Chip papers/unknown_gene.csv")
