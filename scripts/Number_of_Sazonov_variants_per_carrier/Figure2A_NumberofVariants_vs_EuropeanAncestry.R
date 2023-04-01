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

variant_25counts_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/V225SazonovsVariants_Counts.txt"))
variant_allcounts_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/V2AllSazonovsVariants_Counts.txt"))

ancestry_fractions_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/finalAncestryEstimates.2.Q"))
ancestry_fractions_df$ID <- paste(ancestry_fractions_df$ID1,ancestry_fractions_df$ID2,sep="_")

merge_25_df <- merge(variant_25counts_df,ancestry_fractions_df,by="ID")
merge_25_df$sumrange_str <- as.factor(merge_25_df$sumrange)
merge_25_df$cohort <- ifelse(merge_25_df$label == "Case", "AFR Emory IBD","AFR Emory Control")
merge_25_df$disease_status <- ifelse(merge_25_df$label == "Case", 1,0)

merge_all_df <- merge(variant_allcounts_df,ancestry_fractions_df,by="ID")
#merge_all_df$sumrange_str <- as.factor(merge_all_df$sumrange)
merge_all_df$cohort <- ifelse(merge_all_df$label == "Case", "AFR Emory IBD","AFR Emory Control")
merge_all_df$disease_status <- ifelse(merge_all_df$label == "Case", 1,0)
#merge_all_df[nrow(merge_all_df) + 1,] = c("blah",8,"Case","bla","h",0,0,0,0,0.9,0.9,"AFR Emory IBD",1)
merge_all_df$sumrange_str <- as.factor(merge_all_df$sumrange)
#merge_all_df$sumrange <- as.numeric(merge_all_df$sumrange)
#lm_all <- lm(European~sumrange,data = merge_all_df)
#lm_25 <- lm(European~sumrange,data = merge_25_df)

#lm_all_disease_status <- lm(European~sumrange + disease_status,data = merge_all_df)
#lm_25_disease_status <- lm(European~sumrange + disease_status,data = merge_25_df)


colors  <- c("AFR Emory IBD" = "blue4",
             "AFR Emory Control" = "royalblue1")

eur_n_sazonovs_25_box_plt <- ggplot(merge_25_df, aes(x=sumrange_str, y=European,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Fraction of European Ancestry",x="Number of Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
eur_n_sazonovs_25_box_plt
ggsave(file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/Figures/V5NumberofSazonovsVariants_vs_FractionEuropean.pdf",eur_n_sazonovs_25_box_plt,limitsize = FALSE, height = 8, width = 8, dpi=300)


eur_n_sazonovs_all_box_plt <- ggplot(merge_all_df) + geom_boxplot(aes(x = sumrange_str, y = European, fill = cohort),position = position_dodge(preserve = "single")) +
  theme_classic() + labs(y="Fraction of European Ancestry",x="Number of Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors,drop=FALSE) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
eur_n_sazonovs_all_box_plt

final_plt <- grid.arrange(eur_n_sazonovs_25_box_plt,eur_n_sazonovs_all_box_plt,nrow = 1)
ggsave(file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/Figures/V4NumberofSazonovsVariants_vs_FractionEuropean.pdf",final_plt,limitsize = FALSE, height = 6, width = 12.5, dpi=300)


#----------------------------------------
afr_n_sazonovs_25_box_plt <- ggplot(merge_25_df, aes(x=sumrange_str, y=African,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Fraction of African Ancestry",x="Number of Sazonovs variants",title="Sazonovs variants explaining the most variance") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
afr_n_sazonovs_25_box_plt

afr_n_sazonovs_all_box_plt <- ggplot(merge_all_df, aes(x=sumrange_str, y=African,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Fraction of African Ancestry",x="Number of Sazonovs variants",title="All Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
afr_n_sazonovs_all_box_plt

eur_n_sazonovs_all_box_plt <- ggplot(merge_all_df, aes(x=sumrange_str, y=European,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Fraction of European Ancestry",x="Number of Sazonovs variants",title="All Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors,drop=FALSE) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
eur_n_sazonovs_all_box_plt


n_sazonvs_25_box_2plt <- ggplot(merge_25_df, aes(y=sumrange,x=cohort,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Number of Sazonovs variants",title="Sazonovs variants explaining the most variance") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
n_sazonvs_25_box_2plt

n_sazonvs_all_hist_plt <- ggplot(merge_all_df, aes(x=sumrange,fill=cohort)) +
  geom_density(alpha=0.7) + theme_classic() + labs(y="Density",x="Number of Sazonovs variants",title="All Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
n_sazonvs_all_hist_plt

