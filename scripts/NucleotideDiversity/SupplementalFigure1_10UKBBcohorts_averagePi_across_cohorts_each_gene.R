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

gene_map_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/annotations/sazonovs_gene_coordinates.lst"))

gene_cds_length_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/annotations/total_CDS_region_lengths_sazonovs_cds_regions.lst"))

gen_pos_map_df <- NULL

# Loop over all variants and their counts for case and controls 
for(i in 1:nrow(gene_map_df)) 
{
  CHROM <- gene_map_df[i,1]
  POS_START <- gene_map_df[i,2]
  POS_END <- gene_map_df[i,3]
  ENSG <- gene_map_df[i,4]
  SYMBOL <- gene_map_df[i,5]
  
  gene_seq_range_df <- as.data.frame(seq.int(POS_START,POS_END))
  colnames(gene_seq_range_df) <- c("POS")
  gene_seq_range_df$CHROM <- CHROM
  gene_seq_range_df$ENSG <- ENSG
  gene_seq_range_df$SYMBOL <- SYMBOL
  gene_seq_range_df$CHROM_POS <- paste(gene_seq_range_df$CHROM,"-",gene_seq_range_df$POS)
  gene_seq_range_df$CHROM_POS <- gsub(" ","",gene_seq_range_df$CHROM_POS)
    
  gen_pos_map_df = rbind(gen_pos_map_df, gene_seq_range_df)
}

# UKBB IBD 1
ibd1_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort1_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd1_ukbb_df <- ibd1_ukbb_df[complete.cases(ibd1_ukbb_df), ]
ibd1_ukbb_df$CHROM_POS <- paste(ibd1_ukbb_df$CHROM,"-",ibd1_ukbb_df$POS)
ibd1_ukbb_df$CHROM_POS <- gsub(" ","",ibd1_ukbb_df$CHROM_POS)
merge_ibd1_ukbb_df <- merge(ibd1_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd1_ukbb_df <- merge_ibd1_ukbb_df[merge_ibd1_ukbb_df$PI!=0, ]
avgPi_merge_ibd1_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd1_ukbb_df, mean )
n_sites_per_gene_ibd1_ukbb_df <- count(merge_ibd1_ukbb_df, 'SYMBOL')
avgPi_merge_ibd1_ukbb_df <- merge(avgPi_merge_ibd1_ukbb_df,n_sites_per_gene_ibd1_ukbb_df,by="SYMBOL")
avgPi_merge_ibd1_ukbb_df$label <- "UKBB IBD Cohort 1"
avgPi_merge_ibd1_ukbb_df$color <- "palevioletred4"

# UKBB IBD 2
ibd2_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort2_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd2_ukbb_df <- ibd2_ukbb_df[complete.cases(ibd2_ukbb_df), ]
ibd2_ukbb_df$CHROM_POS <- paste(ibd2_ukbb_df$CHROM,"-",ibd2_ukbb_df$POS)
ibd2_ukbb_df$CHROM_POS <- gsub(" ","",ibd2_ukbb_df$CHROM_POS)
merge_ibd2_ukbb_df <- merge(ibd2_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd2_ukbb_df <- merge_ibd2_ukbb_df[merge_ibd2_ukbb_df$PI!=0, ]
avgPi_merge_ibd2_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd2_ukbb_df, mean )
n_sites_per_gene_ibd2_ukbb_df <- count(merge_ibd2_ukbb_df, 'SYMBOL')
avgPi_merge_ibd2_ukbb_df <- merge(avgPi_merge_ibd2_ukbb_df,n_sites_per_gene_ibd2_ukbb_df,by="SYMBOL")
avgPi_merge_ibd2_ukbb_df$label <- "UKBB IBD Cohort 2"
avgPi_merge_ibd2_ukbb_df$color <- "palevioletred4"

# UKBB IBD 3
ibd3_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort3_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd3_ukbb_df <- ibd3_ukbb_df[complete.cases(ibd3_ukbb_df), ]
ibd3_ukbb_df$CHROM_POS <- paste(ibd3_ukbb_df$CHROM,"-",ibd3_ukbb_df$POS)
ibd3_ukbb_df$CHROM_POS <- gsub(" ","",ibd3_ukbb_df$CHROM_POS)
merge_ibd3_ukbb_df <- merge(ibd3_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd3_ukbb_df <- merge_ibd3_ukbb_df[merge_ibd3_ukbb_df$PI!=0, ]
avgPi_merge_ibd3_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd3_ukbb_df, mean )
n_sites_per_gene_ibd3_ukbb_df <- count(merge_ibd3_ukbb_df, 'SYMBOL')
avgPi_merge_ibd3_ukbb_df <- merge(avgPi_merge_ibd3_ukbb_df,n_sites_per_gene_ibd3_ukbb_df,by="SYMBOL")
avgPi_merge_ibd3_ukbb_df$label <- "UKBB IBD Cohort 3"
avgPi_merge_ibd3_ukbb_df$color <- "palevioletred4"

# UKBB IBD 4
ibd4_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort4_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd4_ukbb_df <- ibd4_ukbb_df[complete.cases(ibd4_ukbb_df), ]
ibd4_ukbb_df$CHROM_POS <- paste(ibd4_ukbb_df$CHROM,"-",ibd4_ukbb_df$POS)
ibd4_ukbb_df$CHROM_POS <- gsub(" ","",ibd4_ukbb_df$CHROM_POS)
merge_ibd4_ukbb_df <- merge(ibd4_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd4_ukbb_df <- merge_ibd4_ukbb_df[merge_ibd4_ukbb_df$PI!=0, ]
avgPi_merge_ibd4_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd4_ukbb_df, mean )
n_sites_per_gene_ibd4_ukbb_df <- count(merge_ibd4_ukbb_df, 'SYMBOL')
avgPi_merge_ibd4_ukbb_df <- merge(avgPi_merge_ibd4_ukbb_df,n_sites_per_gene_ibd4_ukbb_df,by="SYMBOL")
avgPi_merge_ibd4_ukbb_df$label <- "UKBB IBD Cohort 4"
avgPi_merge_ibd4_ukbb_df$color <- "palevioletred4"

# UKBB IBD 5
ibd5_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort5_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd5_ukbb_df <- ibd5_ukbb_df[complete.cases(ibd5_ukbb_df), ]
ibd5_ukbb_df$CHROM_POS <- paste(ibd5_ukbb_df$CHROM,"-",ibd5_ukbb_df$POS)
ibd5_ukbb_df$CHROM_POS <- gsub(" ","",ibd5_ukbb_df$CHROM_POS)
merge_ibd5_ukbb_df <- merge(ibd5_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd5_ukbb_df <- merge_ibd5_ukbb_df[merge_ibd5_ukbb_df$PI!=0, ]
avgPi_merge_ibd5_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd5_ukbb_df, mean )
n_sites_per_gene_ibd5_ukbb_df <- count(merge_ibd5_ukbb_df, 'SYMBOL')
avgPi_merge_ibd5_ukbb_df <- merge(avgPi_merge_ibd5_ukbb_df,n_sites_per_gene_ibd5_ukbb_df,by="SYMBOL")
avgPi_merge_ibd5_ukbb_df$label <- "UKBB IBD Cohort 5"
avgPi_merge_ibd5_ukbb_df$color <- "palevioletred4"
# UKBB IBD 6
ibd6_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort6_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd6_ukbb_df <- ibd6_ukbb_df[complete.cases(ibd6_ukbb_df), ]
ibd6_ukbb_df$CHROM_POS <- paste(ibd6_ukbb_df$CHROM,"-",ibd6_ukbb_df$POS)
ibd6_ukbb_df$CHROM_POS <- gsub(" ","",ibd6_ukbb_df$CHROM_POS)
merge_ibd6_ukbb_df <- merge(ibd6_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd6_ukbb_df <- merge_ibd6_ukbb_df[merge_ibd6_ukbb_df$PI!=0, ]
avgPi_merge_ibd6_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd6_ukbb_df, mean )
n_sites_per_gene_ibd6_ukbb_df <- count(merge_ibd6_ukbb_df, 'SYMBOL')
avgPi_merge_ibd6_ukbb_df <- merge(avgPi_merge_ibd6_ukbb_df,n_sites_per_gene_ibd6_ukbb_df,by="SYMBOL")
avgPi_merge_ibd6_ukbb_df$label <- "UKBB IBD Cohort 6"
avgPi_merge_ibd6_ukbb_df$color <- "palevioletred4"

# UKBB IBD 7
ibd7_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort7_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd7_ukbb_df <- ibd7_ukbb_df[complete.cases(ibd7_ukbb_df), ]
ibd7_ukbb_df$CHROM_POS <- paste(ibd7_ukbb_df$CHROM,"-",ibd7_ukbb_df$POS)
ibd7_ukbb_df$CHROM_POS <- gsub(" ","",ibd7_ukbb_df$CHROM_POS)
merge_ibd7_ukbb_df <- merge(ibd7_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd7_ukbb_df <- merge_ibd7_ukbb_df[merge_ibd7_ukbb_df$PI!=0, ]
avgPi_merge_ibd7_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd7_ukbb_df, mean )
n_sites_per_gene_ibd7_ukbb_df <- count(merge_ibd7_ukbb_df, 'SYMBOL')
avgPi_merge_ibd7_ukbb_df <- merge(avgPi_merge_ibd7_ukbb_df,n_sites_per_gene_ibd7_ukbb_df,by="SYMBOL")
avgPi_merge_ibd7_ukbb_df$label <- "UKBB IBD Cohort 7"
avgPi_merge_ibd7_ukbb_df$color <- "palevioletred4"

# UKBB IBD 8
ibd8_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort8_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd8_ukbb_df <- ibd8_ukbb_df[complete.cases(ibd8_ukbb_df), ]
ibd8_ukbb_df$CHROM_POS <- paste(ibd8_ukbb_df$CHROM,"-",ibd8_ukbb_df$POS)
ibd8_ukbb_df$CHROM_POS <- gsub(" ","",ibd8_ukbb_df$CHROM_POS)
merge_ibd8_ukbb_df <- merge(ibd8_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd8_ukbb_df <- merge_ibd8_ukbb_df[merge_ibd8_ukbb_df$PI!=0, ]
avgPi_merge_ibd8_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd8_ukbb_df, mean )
n_sites_per_gene_ibd8_ukbb_df <- count(merge_ibd8_ukbb_df, 'SYMBOL')
avgPi_merge_ibd8_ukbb_df <- merge(avgPi_merge_ibd8_ukbb_df,n_sites_per_gene_ibd8_ukbb_df,by="SYMBOL")
avgPi_merge_ibd8_ukbb_df$label <- "UKBB IBD Cohort 8"
avgPi_merge_ibd8_ukbb_df$color <- "palevioletred4"

# UKBB IBD 9
ibd9_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort9_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd9_ukbb_df <- ibd9_ukbb_df[complete.cases(ibd9_ukbb_df), ]
ibd9_ukbb_df$CHROM_POS <- paste(ibd9_ukbb_df$CHROM,"-",ibd9_ukbb_df$POS)
ibd9_ukbb_df$CHROM_POS <- gsub(" ","",ibd9_ukbb_df$CHROM_POS)
merge_ibd9_ukbb_df <- merge(ibd9_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd_ukbb_df <- merge_ibd9_ukbb_df[merge_ibd9_ukbb_df$PI!=0, ]
avgPi_merge_ibd9_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd9_ukbb_df, mean )
n_sites_per_gene_ibd9_ukbb_df <- count(merge_ibd9_ukbb_df, 'SYMBOL')
avgPi_merge_ibd9_ukbb_df <- merge(avgPi_merge_ibd9_ukbb_df,n_sites_per_gene_ibd9_ukbb_df,by="SYMBOL")
avgPi_merge_ibd9_ukbb_df$label <- "UKBB IBD Cohort 9"
avgPi_merge_ibd9_ukbb_df$color <- "palevioletred4"

# UKBB IBD 10
ibd10_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort10_1774_IBDCases_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
ibd10_ukbb_df <- ibd10_ukbb_df[complete.cases(ibd10_ukbb_df), ]
ibd10_ukbb_df$CHROM_POS <- paste(ibd10_ukbb_df$CHROM,"-",ibd10_ukbb_df$POS)
ibd10_ukbb_df$CHROM_POS <- gsub(" ","",ibd10_ukbb_df$CHROM_POS)
merge_ibd10_ukbb_df <- merge(ibd10_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd10_ukbb_df <- merge_ibd10_ukbb_df[merge_ibd10_ukbb_df$PI!=0, ]
avgPi_merge_ibd10_ukbb_df <- aggregate( PI ~ SYMBOL, merge_ibd10_ukbb_df, mean )
n_sites_per_gene_ibd10_ukbb_df <- count(merge_ibd10_ukbb_df, 'SYMBOL')
avgPi_merge_ibd10_ukbb_df <- merge(avgPi_merge_ibd10_ukbb_df,n_sites_per_gene_ibd10_ukbb_df,by="SYMBOL")
avgPi_merge_ibd10_ukbb_df$label <- "UKBB IBD Cohort 10"
avgPi_merge_ibd10_ukbb_df$color <- "palevioletred4"

# UKBB Control 1
cont1_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort1_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont1_ukbb_df <- cont1_ukbb_df[complete.cases(cont1_ukbb_df), ]
cont1_ukbb_df$CHROM_POS <- paste(cont1_ukbb_df$CHROM,"-",cont1_ukbb_df$POS)
cont1_ukbb_df$CHROM_POS <- gsub(" ","",cont1_ukbb_df$CHROM_POS)
merge_cont1_ukbb_df <- merge(cont1_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont1_ukbb_df <- merge_cont1_ukbb_df[merge_cont1_ukbb_df$PI!=0, ]
avgPi_merge_cont1_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont1_ukbb_df, mean )
n_sites_per_gene_cont1_ukbb_df <- count(merge_cont1_ukbb_df, 'SYMBOL')
avgPi_merge_cont1_ukbb_df <- merge(avgPi_merge_cont1_ukbb_df,n_sites_per_gene_cont1_ukbb_df,by="SYMBOL")
avgPi_merge_cont1_ukbb_df$label <- "UKBB Control Cohort 1"
avgPi_merge_cont1_ukbb_df$color <- "palevioletred3"

# UKBB Control 2
cont2_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort2_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont2_ukbb_df <- cont2_ukbb_df[complete.cases(cont2_ukbb_df), ]
cont2_ukbb_df$CHROM_POS <- paste(cont2_ukbb_df$CHROM,"-",cont2_ukbb_df$POS)
cont2_ukbb_df$CHROM_POS <- gsub(" ","",cont2_ukbb_df$CHROM_POS)
merge_cont2_ukbb_df <- merge(cont2_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont2_ukbb_df <- merge_cont2_ukbb_df[merge_cont2_ukbb_df$PI!=0, ]
avgPi_merge_cont2_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont2_ukbb_df, mean )
n_sites_per_gene_cont2_ukbb_df <- count(merge_cont2_ukbb_df, 'SYMBOL')
avgPi_merge_cont2_ukbb_df <- merge(avgPi_merge_cont2_ukbb_df,n_sites_per_gene_cont2_ukbb_df,by="SYMBOL")
avgPi_merge_cont2_ukbb_df$label <- "UKBB Control Cohort 2"
avgPi_merge_cont2_ukbb_df$color <- "palevioletred3"

# UKBB Control 3
cont3_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort3_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont3_ukbb_df <- cont3_ukbb_df[complete.cases(cont3_ukbb_df), ]
cont3_ukbb_df$CHROM_POS <- paste(cont3_ukbb_df$CHROM,"-",cont3_ukbb_df$POS)
cont3_ukbb_df$CHROM_POS <- gsub(" ","",cont3_ukbb_df$CHROM_POS)
merge_cont3_ukbb_df <- merge(cont3_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont3_ukbb_df <- merge_cont3_ukbb_df[merge_cont3_ukbb_df$PI!=0, ]
avgPi_merge_cont3_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont3_ukbb_df, mean )
n_sites_per_gene_cont3_ukbb_df <- count(merge_cont3_ukbb_df, 'SYMBOL')
avgPi_merge_cont3_ukbb_df <- merge(avgPi_merge_cont3_ukbb_df,n_sites_per_gene_cont3_ukbb_df,by="SYMBOL")
avgPi_merge_cont3_ukbb_df$label <- "UKBB Control Cohort 3"
avgPi_merge_cont3_ukbb_df$color <- "palevioletred3"

# UKBB Control 4
cont4_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort4_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont4_ukbb_df <- cont4_ukbb_df[complete.cases(cont4_ukbb_df), ]
cont4_ukbb_df$CHROM_POS <- paste(cont4_ukbb_df$CHROM,"-",cont4_ukbb_df$POS)
cont4_ukbb_df$CHROM_POS <- gsub(" ","",cont4_ukbb_df$CHROM_POS)
merge_cont4_ukbb_df <- merge(cont4_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont4_ukbb_df <- merge_cont4_ukbb_df[merge_cont4_ukbb_df$PI!=0, ]
avgPi_merge_cont4_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont4_ukbb_df, mean )
n_sites_per_gene_cont4_ukbb_df <- count(merge_cont4_ukbb_df, 'SYMBOL')
avgPi_merge_cont4_ukbb_df <- merge(avgPi_merge_cont4_ukbb_df,n_sites_per_gene_cont4_ukbb_df,by="SYMBOL")
avgPi_merge_cont4_ukbb_df$label <- "UKBB Control Cohort 4"
avgPi_merge_cont4_ukbb_df$color <- "palevioletred3"

# UKBB Control 5
cont5_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort5_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont5_ukbb_df <- cont5_ukbb_df[complete.cases(cont5_ukbb_df), ]
cont5_ukbb_df$CHROM_POS <- paste(cont5_ukbb_df$CHROM,"-",cont5_ukbb_df$POS)
cont5_ukbb_df$CHROM_POS <- gsub(" ","",cont5_ukbb_df$CHROM_POS)
merge_cont5_ukbb_df <- merge(cont5_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont5_ukbb_df <- merge_cont5_ukbb_df[merge_cont5_ukbb_df$PI!=0, ]
avgPi_merge_cont5_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont5_ukbb_df, mean )
n_sites_per_gene_cont5_ukbb_df <- count(merge_cont5_ukbb_df, 'SYMBOL')
avgPi_merge_cont5_ukbb_df <- merge(avgPi_merge_cont5_ukbb_df,n_sites_per_gene_cont5_ukbb_df,by="SYMBOL")
avgPi_merge_cont5_ukbb_df$label <- "UKBB Control Cohort 5"
avgPi_merge_cont5_ukbb_df$color <- "palevioletred3"

# UKBB Control 6
cont6_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort6_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont6_ukbb_df <- cont6_ukbb_df[complete.cases(cont6_ukbb_df), ]
cont6_ukbb_df$CHROM_POS <- paste(cont6_ukbb_df$CHROM,"-",cont6_ukbb_df$POS)
cont6_ukbb_df$CHROM_POS <- gsub(" ","",cont6_ukbb_df$CHROM_POS)
merge_cont6_ukbb_df <- merge(cont6_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont6_ukbb_df <- merge_cont6_ukbb_df[merge_cont6_ukbb_df$PI!=0, ]
avgPi_merge_cont6_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont6_ukbb_df, mean )
n_sites_per_gene_cont6_ukbb_df <- count(merge_cont6_ukbb_df, 'SYMBOL')
avgPi_merge_cont6_ukbb_df <- merge(avgPi_merge_cont6_ukbb_df,n_sites_per_gene_cont6_ukbb_df,by="SYMBOL")
avgPi_merge_cont6_ukbb_df$label <- "UKBB Control Cohort 6"
avgPi_merge_cont6_ukbb_df$color <- "palevioletred3"

# UKBB Control 7
cont7_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort7_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont7_ukbb_df <- cont7_ukbb_df[complete.cases(cont7_ukbb_df), ]
cont7_ukbb_df$CHROM_POS <- paste(cont7_ukbb_df$CHROM,"-",cont7_ukbb_df$POS)
cont7_ukbb_df$CHROM_POS <- gsub(" ","",cont7_ukbb_df$CHROM_POS)
merge_cont7_ukbb_df <- merge(cont7_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont7_ukbb_df <- merge_cont7_ukbb_df[merge_cont7_ukbb_df$PI!=0, ]
avgPi_merge_cont7_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont7_ukbb_df, mean )
n_sites_per_gene_cont7_ukbb_df <- count(merge_cont7_ukbb_df, 'SYMBOL')
avgPi_merge_cont7_ukbb_df <- merge(avgPi_merge_cont7_ukbb_df,n_sites_per_gene_cont7_ukbb_df,by="SYMBOL")
avgPi_merge_cont7_ukbb_df$label <- "UKBB Control Cohort 7"
avgPi_merge_cont7_ukbb_df$color <- "palevioletred3"

# UKBB Control 8
cont8_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort8_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont8_ukbb_df <- cont8_ukbb_df[complete.cases(cont8_ukbb_df), ]
cont8_ukbb_df$CHROM_POS <- paste(cont8_ukbb_df$CHROM,"-",cont8_ukbb_df$POS)
cont8_ukbb_df$CHROM_POS <- gsub(" ","",cont8_ukbb_df$CHROM_POS)
merge_cont8_ukbb_df <- merge(cont8_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont8_ukbb_df <- merge_cont8_ukbb_df[merge_cont8_ukbb_df$PI!=0, ]
avgPi_merge_cont8_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont8_ukbb_df, mean )
n_sites_per_gene_cont8_ukbb_df <- count(merge_cont8_ukbb_df, 'SYMBOL')
avgPi_merge_cont8_ukbb_df <- merge(avgPi_merge_cont8_ukbb_df,n_sites_per_gene_cont8_ukbb_df,by="SYMBOL")
avgPi_merge_cont8_ukbb_df$label <- "UKBB Control Cohort 8"
avgPi_merge_cont8_ukbb_df$color <- "palevioletred3"

# UKBB Control 9
cont9_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort9_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont9_ukbb_df <- cont9_ukbb_df[complete.cases(cont9_ukbb_df), ]
cont9_ukbb_df$CHROM_POS <- paste(cont9_ukbb_df$CHROM,"-",cont9_ukbb_df$POS)
cont9_ukbb_df$CHROM_POS <- gsub(" ","",cont9_ukbb_df$CHROM_POS)
merge_cont9_ukbb_df <- merge(cont9_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont9_ukbb_df <- merge_cont9_ukbb_df[merge_cont9_ukbb_df$PI!=0, ]
avgPi_merge_cont9_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont9_ukbb_df, mean )
n_sites_per_gene_cont9_ukbb_df <- count(merge_cont9_ukbb_df, 'SYMBOL')
avgPi_merge_cont9_ukbb_df <- merge(avgPi_merge_cont9_ukbb_df,n_sites_per_gene_cont9_ukbb_df,by="SYMBOL")
avgPi_merge_cont9_ukbb_df$label <- "UKBB Control Cohort 9"
avgPi_merge_cont9_ukbb_df$color <- "palevioletred3"

# UKBB Control 10
cont10_ukbb_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/PI_VCF_Sazonovs_CDSregions_Cohort10_1644_IBDControls_k50_k51_UKBB_Caucasian_SexChecked.pheno.vcf.sites.pi"))
cont10_ukbb_df <- cont10_ukbb_df[complete.cases(cont10_ukbb_df), ]
cont10_ukbb_df$CHROM_POS <- paste(cont10_ukbb_df$CHROM,"-",cont10_ukbb_df$POS)
cont10_ukbb_df$CHROM_POS <- gsub(" ","",cont10_ukbb_df$CHROM_POS)
merge_cont10_ukbb_df <- merge(cont10_ukbb_df,gen_pos_map_df,by="CHROM_POS")
merge_cont10_ukbb_df <- merge_cont10_ukbb_df[merge_cont10_ukbb_df$PI!=0, ]
avgPi_merge_cont10_ukbb_df <- aggregate( PI ~ SYMBOL, merge_cont10_ukbb_df, mean )
n_sites_per_gene_cont10_ukbb_df <- count(merge_cont10_ukbb_df, 'SYMBOL')
avgPi_merge_cont10_ukbb_df <- merge(avgPi_merge_cont10_ukbb_df,n_sites_per_gene_cont10_ukbb_df,by="SYMBOL")
avgPi_merge_cont10_ukbb_df$label <- "UKBB Control Cohort 10"
avgPi_merge_cont10_ukbb_df$color <- "palevioletred3"

# AFR IBD
ibd_afr_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/V2_Pi_IBDCase_sazonovs_genes_CDSregions_big_daly_v3_AFR_WGS.sites.pi"))
ibd_afr_df <- ibd_afr_df[complete.cases(ibd_afr_df), ]
ibd_afr_df$CHROM_POS <- paste(ibd_afr_df$CHROM,"-",ibd_afr_df$POS)
ibd_afr_df$CHROM_POS <- gsub(" ","",ibd_afr_df$CHROM_POS)
merge_ibd_afr_df <- merge(ibd_afr_df,gen_pos_map_df,by="CHROM_POS")
merge_ibd_afr_df <- merge_ibd_afr_df[merge_ibd_afr_df$PI!=0, ]
avgPi_merge_ibd_afr_df <- aggregate( PI ~ SYMBOL, merge_ibd_afr_df, mean )
n_sites_per_gene_ibd_afr_df <- count(merge_ibd_afr_df, 'SYMBOL')
avgPi_merge_ibd_afr_df <- merge(avgPi_merge_ibd_afr_df,n_sites_per_gene_ibd_afr_df,by="SYMBOL")
avgPi_merge_ibd_afr_df$label <- "AFR Emory IBD"
avgPi_merge_ibd_afr_df$color <- "blue4"

# AFR Control
control_afr_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/nucleotide_diversity/V2_Pi_IBDControl_sazonovs_genes_CDSregions_big_daly_v3_AFR_WGS.sites.pi"))
control_afr_df <- control_afr_df[complete.cases(control_afr_df), ]
control_afr_df$CHROM_POS <- paste(control_afr_df$CHROM,"-",control_afr_df$POS)
control_afr_df$CHROM_POS <- gsub(" ","",control_afr_df$CHROM_POS)
merge_control_afr_df <- merge(control_afr_df,gen_pos_map_df,by="CHROM_POS")
merge_control_afr_df <- merge_control_afr_df[merge_control_afr_df$PI!=0, ]
avgPi_merge_control_afr_df<- aggregate( PI ~ SYMBOL, merge_control_afr_df, mean )
n_sites_per_gene_icontrol_afr_df <- count(merge_control_afr_df, 'SYMBOL')
avgPi_merge_control_afr_df <- merge(avgPi_merge_control_afr_df,n_sites_per_gene_icontrol_afr_df,by="SYMBOL")
avgPi_merge_control_afr_df$label <- "AFR Emory Control"
avgPi_merge_control_afr_df$color <- "royalblue1"

final_df <- rbind(avgPi_merge_ibd_afr_df,avgPi_merge_control_afr_df,avgPi_merge_ibd1_ukbb_df,avgPi_merge_cont1_ukbb_df,avgPi_merge_ibd2_ukbb_df,avgPi_merge_cont2_ukbb_df,avgPi_merge_ibd3_ukbb_df,avgPi_merge_cont3_ukbb_df,avgPi_merge_ibd4_ukbb_df,avgPi_merge_cont4_ukbb_df,avgPi_merge_ibd5_ukbb_df,avgPi_merge_cont5_ukbb_df,avgPi_merge_ibd6_ukbb_df,avgPi_merge_cont6_ukbb_df,avgPi_merge_ibd7_ukbb_df,avgPi_merge_cont7_ukbb_df,avgPi_merge_ibd8_ukbb_df,avgPi_merge_cont8_ukbb_df,avgPi_merge_ibd9_ukbb_df,avgPi_merge_cont9_ukbb_df,avgPi_merge_ibd10_ukbb_df,avgPi_merge_cont10_ukbb_df)

colors  <- c("AFR Emory IBD" = "blue4",
             "AFR Emory Control" = "royalblue1",
             "UKBB IBD Cohort 1" = "palevioletred4",
             "UKBB Control Cohort 1" = "palevioletred3",
             "UKBB IBD Cohort 2" = "palevioletred4",
             "UKBB Control Cohort 2" = "palevioletred3",
             "UKBB IBD Cohort 3" = "palevioletred4",
             "UKBB Control Cohort 3" = "palevioletred3",
             "UKBB IBD Cohort 4" = "palevioletred4",
             "UKBB Control Cohort 4" = "palevioletred3",
             "UKBB IBD Cohort 5" = "palevioletred4",
             "UKBB Control Cohort 5" = "palevioletred3",
             "UKBB IBD Cohort 6" = "palevioletred4",
             "UKBB Control Cohort 6" = "palevioletred3",
             "UKBB IBD Cohort 7" = "palevioletred4",
             "UKBB Control Cohort 7" = "palevioletred3",
             "UKBB IBD Cohort 8" = "palevioletred4",
             "UKBB Control Cohort 8" = "palevioletred3",
             "UKBB IBD Cohort 9" = "palevioletred4",
             "UKBB Control Cohort 9" = "palevioletred3",
             "UKBB IBD Cohort 10" = "palevioletred4",
             "UKBB Control Cohort 10" = "palevioletred3")


#final_df1$label2 <- factor(final_df1$label,levels = c("AFR Emory IBD", "AFR Emory Control", "GNOMAD AFR", "UKBB IBD (Genebass)","UKBB Control (Genebass)","EUR Emory Pediatric IBD","GNOMAD NFE"))
final_df$label2 <- factor(final_df$label,levels = c("AFR Emory IBD", "AFR Emory Control","UKBB IBD Cohort 1","UKBB Control Cohort 1","UKBB IBD Cohort 2","UKBB Control Cohort 2","UKBB IBD Cohort 3","UKBB Control Cohort 3","UKBB IBD Cohort 4","UKBB Control Cohort 4","UKBB IBD Cohort 5","UKBB Control Cohort 5","UKBB IBD Cohort 6","UKBB Control Cohort 6","UKBB IBD Cohort 7","UKBB Control Cohort 7","UKBB IBD Cohort 8","UKBB Control Cohort 8","UKBB IBD Cohort 9","UKBB Control Cohort 9","UKBB IBD Cohort 10","UKBB Control Cohort 10"))
#final_df2 <- na.omit(final_df2, cols = c("AF", "SYMBOL"))

# Merge CDS region length 
final_df <- merge(gene_cds_length_df,final_df,by="SYMBOL")
final_df$adjusted_pi <- ((final_df$PI * final_df$freq) + (final_df$cds_region_length - final_df$freq)*0) / final_df$cds_region_length

# Get sazonovs annotations 
sazonov_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/sazonovs/sazonov_significant_variants.txt"))
sazonov_df <- sazonov_df[(sazonov_df$FigureStatus == "Known causal candidate") | (sazonov_df$FigureStatus == "New variant in known locus") | (sazonov_df$FigureStatus == "New Locus"), ]
sazonov_df$variant_id <- paste(sazonov_df$Chr,sazonov_df$Pos,sazonov_df$A0,sazonov_df$A1,sep="-")
sazonov_df <- sazonov_df[c('variant_id','or_meta','gene')]
colnames(sazonov_df) <- c('variant_id','or_meta','SYMBOL')

risk_sazonov_df <- sazonov_df[sazonov_df$or_meta > 1, ]
risk_sazonov_df <- risk_sazonov_df[c('SYMBOL')]
risk_sazonov_df <- unique(risk_sazonov_df)

protective_sazonov_df <- sazonov_df[sazonov_df$or_meta < 1, ]
protective_sazonov_df <- protective_sazonov_df[c('SYMBOL')]
protective_sazonov_df <- unique(protective_sazonov_df)

final_risk_sazonov_df <- merge(risk_sazonov_df,final_df,by="SYMBOL")
final_protective_sazonov_df <- merge(protective_sazonov_df,final_df,by="SYMBOL")
final_risk_sazonov_df$new_SYMBOL <- factor(final_risk_sazonov_df$SYMBOL, levels=c("PTAFR", "RELA","SDF2L1","CCR7", "HGFAC","IL10RA","PDLIM5","SLC39A8","DOK2","NOD2"))
final_protective_sazonov_df$new_SYMBOL <- factor(final_protective_sazonov_df$SYMBOL, levels=c("TAGAP","IL23R","CARD9","TYK2"))

risk_adjpi_bar_plt <- ggplot(final_risk_sazonov_df, aes(x = as.factor(label2), y = adjusted_pi,fill=as.factor(label2))) +
  geom_bar(position = 'dodge', stat = 'identity') + ylim(0,0.0011) + 
  labs(title="Risk genes",y="Adjusted average pi") +
  scale_fill_manual(name="Cohort",values = colors) + theme_classic() + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  facet_wrap(~new_SYMBOL,scale="free",nrow = 2)
risk_adjpi_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/V6_SupplementalFigure1A_riskgenes_adjusted_avg_pi_ukbb_downsamples.pdf",limitsize = FALSE, height = 12, width = 30, dpi=300)

protective_adjpi_bar_plt <- ggplot(final_protective_sazonov_df, aes(x = as.factor(label2), y = adjusted_pi,fill=as.factor(label2))) +
  geom_bar(position = 'dodge', stat = 'identity') + ylim(0,0.0011) + 
  labs(title="Protective genes",y="Adjusted average pi") +
  scale_fill_manual(name="Cohort",values = colors) + theme_classic() + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  facet_wrap(~new_SYMBOL,scale="free",nrow=1)
protective_adjpi_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/V6_SupplementalFigure1B_protectivegenes_adjusted_avg_pi_ukbb_downsamples.pdf",limitsize = FALSE, height = 8, width = 25, dpi=300)


# Formatting data for Table 2
final_df2 <- final_df[c('SYMBOL','adjusted_pi','label')]

afr_ibd_df <- final_df2[final_df2$label=="AFR Emory IBD", ]
afr_ibd_df <- afr_ibd_df[c('SYMBOL','adjusted_pi')]
colnames(afr_ibd_df) <- c('SYMBOL','AFR_Emory_IBD_adjusted_pi')

afr_cont_df <- final_df2[final_df2$label=="AFR Emory Control", ]
afr_cont_df <- afr_cont_df[c('SYMBOL','adjusted_pi')]
colnames(afr_cont_df) <- c('SYMBOL','AFR_Emory_Control_adjusted_pi')

ukb_ibd_df1 <- final_df2[final_df2$label=="UKBB IBD Cohort 1", ]
ukb_ibd_df1 <- ukb_ibd_df1[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df1) <- c('SYMBOL','UKBB_IBD_Cohort_1_adjusted_pi')

ukb_cont_df1 <- final_df2[final_df2$label=="UKBB Control Cohort 1", ]
ukb_cont_df1 <- ukb_cont_df1[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df1) <- c('SYMBOL','UKBB_Control_Cohort_1_adjusted_pi')

ukb_ibd_df2 <- final_df2[final_df2$label=="UKBB IBD Cohort 2", ]
ukb_ibd_df2 <- ukb_ibd_df2[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df2) <- c('SYMBOL','UKBB_IBD_Cohort_2_adjusted_pi')

ukb_cont_df2 <- final_df2[final_df2$label=="UKBB Control Cohort 2", ]
ukb_cont_df2 <- ukb_cont_df2[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df2) <- c('SYMBOL','UKBB_Control_Cohort_2_adjusted_pi')

ukb_ibd_df3 <- final_df2[final_df2$label=="UKBB IBD Cohort 3", ]
ukb_ibd_df3 <- ukb_ibd_df3[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df3) <- c('SYMBOL','UKBB_IBD_Cohort_3_adjusted_pi')

ukb_cont_df3<- final_df2[final_df2$label=="UKBB Control Cohort 3", ]
ukb_cont_df3 <- ukb_cont_df3[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df3) <- c('SYMBOL','UKBB_Control_Cohort_3_adjusted_pi')

ukb_ibd_df4 <- final_df2[final_df2$label=="UKBB IBD Cohort 4", ]
ukb_ibd_df4 <- ukb_ibd_df4[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df4) <- c('SYMBOL','UKBB_IBD_Cohort_4_adjusted_pi')

ukb_cont_df4<- final_df2[final_df2$label=="UKBB Control Cohort 4", ]
ukb_cont_df4 <- ukb_cont_df4[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df4) <- c('SYMBOL','UKBB_Control_Cohort_4_adjusted_pi')

ukb_ibd_df5 <- final_df2[final_df2$label=="UKBB IBD Cohort 5", ]
ukb_ibd_df5 <- ukb_ibd_df5[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df5) <- c('SYMBOL','UKBB_IBD_Cohort_5_adjusted_pi')

ukb_cont_df5<- final_df2[final_df2$label=="UKBB Control Cohort 5", ]
ukb_cont_df5 <- ukb_cont_df5[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df5) <- c('SYMBOL','UKBB_Control_Cohort_5_adjusted_pi')

ukb_ibd_df6 <- final_df2[final_df2$label=="UKBB IBD Cohort 6", ]
ukb_ibd_df6 <- ukb_ibd_df6[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df6) <- c('SYMBOL','UKBB_IBD_Cohort_6_adjusted_pi')

ukb_cont_df6<- final_df2[final_df2$label=="UKBB Control Cohort 6", ]
ukb_cont_df6 <- ukb_cont_df6[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df6) <- c('SYMBOL','UKBB_Control_Cohort_6_adjusted_pi')

ukb_ibd_df7 <- final_df2[final_df2$label=="UKBB IBD Cohort 7", ]
ukb_ibd_df7 <- ukb_ibd_df7[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df7) <- c('SYMBOL','UKBB_IBD_Cohort_7_adjusted_pi')

ukb_cont_df7<- final_df2[final_df2$label=="UKBB Control Cohort 7", ]
ukb_cont_df7 <- ukb_cont_df7[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df7) <- c('SYMBOL','UKBB_Control_Cohort_7_adjusted_pi')

ukb_ibd_df8 <- final_df2[final_df2$label=="UKBB IBD Cohort 8", ]
ukb_ibd_df8 <- ukb_ibd_df8[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df8) <- c('SYMBOL','UKBB_IBD_Cohort_8_adjusted_pi')

ukb_cont_df8<- final_df2[final_df2$label=="UKBB Control Cohort 8", ]
ukb_cont_df8 <- ukb_cont_df8[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df8) <- c('SYMBOL','UKBB_Control_Cohort_8_adjusted_pi')

ukb_ibd_df9 <- final_df2[final_df2$label=="UKBB IBD Cohort 9", ]
ukb_ibd_df9 <- ukb_ibd_df9[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df9) <- c('SYMBOL','UKBB_IBD_Cohort_9_adjusted_pi')

ukb_cont_df9<- final_df2[final_df2$label=="UKBB Control Cohort 9", ]
ukb_cont_df9 <- ukb_cont_df9[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df9) <- c('SYMBOL','UKBB_Control_Cohort_9_adjusted_pi')

ukb_ibd_df10 <- final_df2[final_df2$label=="UKBB IBD Cohort 10", ]
ukb_ibd_df10 <- ukb_ibd_df10[c('SYMBOL','adjusted_pi')]
colnames(ukb_ibd_df10) <- c('SYMBOL','UKBB_IBD_Cohort_10_adjusted_pi')

ukb_cont_df10<- final_df2[final_df2$label=="UKBB Control Cohort 10", ]
ukb_cont_df10 <- ukb_cont_df10[c('SYMBOL','adjusted_pi')]
colnames(ukb_cont_df10) <- c('SYMBOL','UKBB_Control_Cohort_10_adjusted_pi')

final_df_list <- list(afr_ibd_df,afr_cont_df,ukb_ibd_df1,ukb_cont_df1,ukb_ibd_df2,ukb_cont_df2,ukb_ibd_df3,ukb_cont_df3,ukb_ibd_df4,ukb_cont_df4,ukb_ibd_df5,ukb_cont_df5,ukb_ibd_df6,ukb_cont_df6,ukb_ibd_df7,ukb_cont_df7,ukb_ibd_df8,ukb_cont_df8,ukb_ibd_df9,ukb_cont_df9,ukb_ibd_df10,ukb_cont_df10)

final_df3 <- final_df_list %>% reduce(full_join, by='SYMBOL')
write.table(final_df3, file='/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/tables/Table2.tsv', quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)


#--------------------------------------------
bar_plt <- ggplot(final_df, aes(x = as.factor(label2), y = PI,fill=as.factor(label2))) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(title="Average nucleotide diversity (pi) for each Sazonovs et.al. genes CDS regions across IBD/control cohorts",y="Average pi") +
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5),legend.position="bottom",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  facet_wrap(~SYMBOL,scale="free")
bar_plt
ggsave("/Users/courtneyastore/Desktop/EUR_AFR_IBD_RareVar_Comparison_project/TajimasD_NucleotideDiversity_SazonovsGene/avg_pi_ukbb_downsamples.png",limitsize = FALSE, height = 35, width = 42, dpi=300)

