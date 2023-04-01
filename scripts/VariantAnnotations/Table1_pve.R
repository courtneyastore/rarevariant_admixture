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
library(ggpubr)
options(ggrepel.max.overlaps = Inf)

# Extract OR for EUR IBD
sazonov_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/sazonovs/sazonov_significant_variants.txt"))
sazonov_df <- sazonov_df[(sazonov_df$FigureStatus == "Known causal candidate") | (sazonov_df$FigureStatus == "New variant in known locus") | (sazonov_df$FigureStatus == "New Locus"), ]

sazonov_df$variant_id <- paste(sazonov_df$Chr,sazonov_df$Pos,sazonov_df$A0,sazonov_df$A1,sep="-")
extract_df <- sazonov_df[c('variant_id','FigureStatus')]
sazonov_pruned_df <- sazonov_df[c('variant_id','or_meta')]
colnames(sazonov_pruned_df) <- c('variant_id','or_eur')

# Get variant annotations
sazonov_annotations_df <- sazonov_df[c('variant_id','gene','AminoAcidSub','consequence','FigureStatus','or_meta')]

risk_sazonov_df <- sazonov_annotations_df[sazonov_annotations_df$or_meta > 1, ]
risk_sazonov_df$risk_protective <- "Risk"
protective_sazonov_df <- sazonov_annotations_df[sazonov_annotations_df$or_meta < 1, ]
protective_sazonov_df$risk_protective <- "Protective"

final_sazonov_annotations_df <- rbind(risk_sazonov_df,protective_sazonov_df)

# Extract OR for AFR IBD
fishers_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/Fishers_test_results/V2MASTER_FishersExactTest_AcrossIBDControl_Cohorts.txt"))
afr_fishers_df <- fishers_df[fishers_df$label == "AFR Emory", ]
afr_fishers_df <- afr_fishers_df[c('variant_id','or')]
colnames(afr_fishers_df) <- c('variant_id','or_afr')

# Merge OR for AFR IBD and EUR IBD on variant_id
final_or_df <- merge(sazonov_pruned_df,afr_fishers_df,by="variant_id",all.x = TRUE)

# Extract AF for AFR and EUR control (GNOMAD NFE and GNOMAD AFR)
control_gnomad_af_df <- sazonov_df[c('variant_id','gnomad_AFR_AF','gnomad_NFE_AF')]

# Merge OR with control AFs from GNOMAD
final_or_af_df <- merge(final_or_df,control_gnomad_af_df,by="variant_id")

# Calculate PVE AFR
final_or_af_df$afr_pve <- 200*final_or_af_df$gnomad_AFR_AF*(1-final_or_af_df$gnomad_AFR_AF) * (log(final_or_af_df$or_afr))^2

# Calculate PVE NFE
final_or_af_df$eur_pve <- 200*final_or_af_df$gnomad_NFE_AF*(1-final_or_af_df$gnomad_NFE_AF) * (log(final_or_af_df$or_eur))^2

# Merge annotations to final table
final_or_af_annotations_df <- merge(final_sazonov_annotations_df,final_or_af_df,by="variant_id")

final_or_af_annotations_df$or_eur <- round(final_or_af_annotations_df$or_eur,digits = 4)
final_or_af_annotations_df$or_afr <- round(final_or_af_annotations_df$or_afr,digits = 4)
final_or_af_annotations_df$afr_pve <- round(final_or_af_annotations_df$afr_pve,digits = 4)
final_or_af_annotations_df$eur_pve <- round(final_or_af_annotations_df$eur_pve,digits = 4)

write.table(final_or_af_annotations_df,"/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/tables/V2table1.tsv",sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

