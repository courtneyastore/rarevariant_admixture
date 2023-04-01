library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library("ggpubr")
library(plyr)
library(stringr)
library(gridExtra)
library(ggrepel)
library(viridis)
library(hrbrthemes)
options(ggrepel.max.overlaps = Inf)

sazonovs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/Formatted_sazonov_significant_variants.txt"))
sazonovs_df <- sazonovs_df[(sazonovs_df$FigureStatus == "Known causal candidate") | (sazonovs_df$FigureStatus == "New variant in known locus") | (sazonovs_df$FigureStatus == "New Locus"), ]
sazonovs_df$varID <- paste(sazonovs_df$Chr,"_",sazonovs_df$Pos,sep="")
sazonovs_df$varID <- gsub(" ","",sazonovs_df$varID)
sazonovs_df$varID <- gsub(",","",sazonovs_df$varID)

sazonovs_df$EUR_case_counts <- sazonovs_df$`Case_AF1(TWIST.SAIGE)` * (2*6109)
sazonovs_df$EUR_control_counts <- sazonovs_df$`Control_AF1(TWIST.SAIGE)` * (2*6064)

sazonovs_df$EUR_case_counts_Nextera <- sazonovs_df$`Case_AF1(Nextera.SAIGE)` * (2*11125)
sazonovs_df$EUR_control_counts_Nextera <- sazonovs_df$`Control_AF1(Nextera.SAIGE)` * (2*25145)

sazonovs_df <- sazonovs_df[c('EUR_case_counts','EUR_control_counts','EUR_case_counts_Nextera','EUR_control_counts_Nextera','varID','gene','AminoAcidSub')]

afr_carrier_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/AFREmory_Sazonovs_GT.tsv"))
afr_carrier_df$ID <- paste(afr_carrier_df$ID1,"_",afr_carrier_df$ID2,sep="")
#afr_carrier_df <- afr_carrier_df[c('ID','SNP','var_type','label')]

afr_have_var_df <- afr_carrier_df[(afr_carrier_df$var_type == "het" | afr_carrier_df$var_type == "hom_alt"), ]
afr_have_var_df <- afr_have_var_df[c('ID','varID','A1','A2','label')]
afr_have_var_df <- unique(afr_have_var_df)

case_afr_have_var_df <- afr_have_var_df[afr_have_var_df$label == "Case", ]
afr_case_counts <- count(case_afr_have_var_df, 'varID')
colnames(afr_case_counts) <- c('varID','AFR_case_counts')

control_afr_have_var_df <- afr_have_var_df[afr_have_var_df$label == "Control", ]
afr_control_counts <- count(control_afr_have_var_df, 'varID')
colnames(afr_control_counts) <- c('varID','AFR_control_counts')

afr_df <- merge(afr_case_counts,afr_control_counts,by="varID")

afr_eur_df <- merge(afr_df,sazonovs_df,by="varID")

twist_fishers_df <- NULL

for (i in 1:nrow(afr_eur_df)){
  variant <- afr_eur_df[i,1]
  total_afr_case <- 1774
  total_afr_control <- 1644
  AFR_case_counts <- afr_eur_df[i,2]
  
  total_eur_case <- 6109
  total_eur_control <- 6064
  EUR_case_counts <- afr_eur_df[i,4]
  
  gene <- afr_eur_df[i,8]
  aa_sub <- afr_eur_df[i,9]
  
  fisher_result <- fisher.test(rbind(c(AFR_case_counts,total_afr_case),c(EUR_case_counts,total_eur_case)))
  p = fisher_result$p.value
  or = as.numeric(gsub("odds ratio","",fisher_result$estimate))
  lo_ci_or = fisher_result$conf.int[1]
  hi_ci_or = fisher_result$conf.int[2]
  
  twist_fishers_df = rbind(twist_fishers_df, data.frame(variant,lo_ci_or,hi_ci_or,or,p,aa_sub,gene))  
  
}

write.table(twist_fishers_df,file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/AA_vs_EURTwist_fishers_test.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

nextera_fishers_df <- NULL

for (i in 1:nrow(afr_eur_df)){
  variant <- afr_eur_df[i,1]
  total_afr_case <- 1774
  total_afr_control <- 1644
  AFR_case_counts <- afr_eur_df[i,2]
  
  total_eur_case <- 11125 
  total_eur_control <- 25145 
  EUR_case_counts <- afr_eur_df[i,6]
  
  gene <- afr_eur_df[i,8]
  aa_sub <- afr_eur_df[i,9]
  
  fisher_result <- fisher.test(rbind(c(AFR_case_counts,total_afr_case),c(EUR_case_counts,total_eur_case)))
  p = fisher_result$p.value
  or = as.numeric(gsub("odds ratio","",fisher_result$estimate))
  lo_ci_or = fisher_result$conf.int[1]
  hi_ci_or = fisher_result$conf.int[2]
  
  nextera_fishers_df = rbind(nextera_fishers_df, data.frame(variant,lo_ci_or,hi_ci_or,or,p,aa_sub,gene))  
  
}
write.table(nextera_fishers_df,file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/AA_vs_EURNextera_fishers_test.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


