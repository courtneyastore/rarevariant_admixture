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

casecontrol_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/AllSamplesGnomixData/Comphetinput_MASTER_case_control_sample_annotations.fam"))
casecontrol_df$ID <- gsub("\\..*","",casecontrol_df$ID)
casecontrol_df <- unique(casecontrol_df)

ancestry_fractions_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/finalAncestryEstimates.2.Q"))
ancestry_fractions_df$ID <- paste(ancestry_fractions_df$ID1,ancestry_fractions_df$ID2,sep="_")

sini_afr_ancestry_fractions_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/Afr_Prop_AAPRS_AAweights.txt"))
sini_afr_ancestry_fractions_df$ID <- paste(sini_afr_ancestry_fractions_df$FID,sini_afr_ancestry_fractions_df$IID,sep="_")
afr_proportions_df <- sini_afr_ancestry_fractions_df[c('ID','africanProportion')]
afr_proportions_df$africanProportion <- afr_proportions_df$africanProportion / 100
sini_afr_ancestry_fractions_df <- sini_afr_ancestry_fractions_df[c('ID','PRS')]
colnames(sini_afr_ancestry_fractions_df) <- c('ID','PRS')
sini_afr_ancestry_fractions_df$label <- "African"
sini_afr_ancestry_fractions_df <- merge(sini_afr_ancestry_fractions_df,afr_proportions_df,by="ID")

afr_lm <- lm(africanProportion~PRS,data=sini_afr_ancestry_fractions_df)
summary(afr_lm)

sini_eur_ancestry_fractions_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/Afr_Prop_AAPRS_EURweights.txt"))
sini_eur_ancestry_fractions_df$ID <- paste(sini_eur_ancestry_fractions_df$FID,sini_eur_ancestry_fractions_df$IID,sep="_")
sini_eur_ancestry_fractions_df <- sini_eur_ancestry_fractions_df[c('ID','PRS')]
colnames(sini_eur_ancestry_fractions_df) <- c('ID','PRS')
sini_eur_ancestry_fractions_df$label <- "European"
sini_eur_ancestry_fractions_df <- merge(sini_eur_ancestry_fractions_df,afr_proportions_df,by="ID")

eur_lm <- lm(africanProportion~PRS,data=sini_eur_ancestry_fractions_df)
summary(eur_lm)$r.squared

final_df <- rbind(sini_afr_ancestry_fractions_df,sini_eur_ancestry_fractions_df)

#merge_df <- merge(final_df,afr_proportions_df,by="ID")
# expression("Fraction of African Ancestry ~ African weights; p-value < 2e-16, r = -0.42,"~ R^2~"= 0.1887"~"\nFraction of African Ancestry ~ European weights; p-value < 2e-16, r = -0.16,"~ R^2~"= 0.0251")) + 

colors  <- c("African" = "blue4",
             "European" = "palevioletred4")
# "Fraction of African Ancestry ~ African weights; p-value < 2e-16, r = -0.42, R^2 = 0.1887 Fraction of African Ancestry ~ European weights; p-value < 2e-16, r = -0.16, R^2 =0.0251"
plt <- ggplot(final_df, aes(y=PRS, x=africanProportion,color=label,fill=label)) +
  geom_point(size=3,alpha=0.2) + theme_classic() + 
  labs(x="Fraction of African Ancestry",y="PRS") + 
  scale_color_manual(name="Ancestry weights",values = colors) + geom_smooth(method='lm', formula= y~x) + 
  scale_fill_manual(name="Ancestry weights",values=colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
plt
ggsave(file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/Figures/V3_Figure3B.pdf",plt,limitsize = FALSE, height = 7, width = 8, dpi=300)



















merge_df <- merge(ancestry_fractions_df,sini_afr_ancestry_fractions_df,by="ID")




prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/ResultFiles/PRSdeLange_HWE1e_9_VariantMissingness0.05_big_daly_v3.sscore"))
prs_df$ID <- paste(prs_df$FID,prs_df$IID,sep="_")
merge_df <- merge(prs_df,ancestry_fractions_df,by="ID")
merge_df$PRS <- scale(merge_df$SCORE1_SUM)

merge_df2 <- merge(merge_df,casecontrol_df,by="ID")

merge_df2$cohort <- ifelse(merge_df2$Pop == "Case", "AFR Emory IBD","AFR Emory Control")
merge_df2$disease_status <- ifelse(merge_df2$Pop == "Case", 1,0)


model <- lm(European~PRS, data = merge_df2)
model2 <- lm(African~PRS, data = merge_df2)

colors  <- c("AFR Emory IBD" = "blue4",
             "AFR Emory Control" = "royalblue1")

plt <- ggplot(merge_df2, aes(x=PRS, y=European,color=cohort)) +
  geom_point(size=3,alpha=0.5) + theme_classic() + labs(y="Fraction of European Ancestry",x="de Lange PRS") + 
  scale_color_manual(name="Cohort",values = colors) + geom_smooth(method='lm', formula= y~x) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
plt


plt2 <- ggplot(merge_df2, aes(y=PRS, x=European)) +
  geom_point(size=3) + theme_classic() + labs(x="Fraction of European Ancestry",y="de Lange PRS",title="Fraction of European ancestry ~ de Lange PRS\ncorrelation=0.013, pval=2.11e-11, OR=1.01") + 
  geom_smooth(method='lm', formula= y~x) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
plt2
ggsave(file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/Figures/V1.FractionEUR_vs_deLange_PRS.pdf",plt2,limitsize = FALSE, height = 8, width = 8, dpi=300)

# PRS vs number of Sazonovs variants
merge_df3 <- merge(merge_df2,variant_25counts_df,by="ID",all.x=TRUE)
merge_df3$carrier_status <- ifelse(is.na(merge_df3$sumrange), "Non-carrier","Carrier")

my_comparisons <- list( c("Carrier","Non-carrier") )
colors_carriers  <- c("Non-carrier" = "purple",
             "Carrier" = "green")
plt3 <- ggplot(merge_df3, aes(x=carrier_status, y=PRS)) +
  geom_boxplot(aes(fill=cohort)) + 
  scale_fill_manual(name="Cohort",values = colors) + 
  facet_grid(~cohort) +
  stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test") + 
  theme_classic() + labs(y="de Lange PRS",x="Rare variant carrier status",title="Sazonovs variants explaining the most variance") + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
plt3











afr_n_sazonovs_all_box_plt <- ggplot(merge_all_df, aes(x=sumrange_str, y=African,fill=cohort)) +
  geom_boxplot() + theme_classic() + labs(y="Fraction of African Ancestry",x="Number of Sazonovs variants",title="All Sazonovs variants") + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom",axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
afr_n_sazonovs_all_box_plt


final_plt <- grid.arrange(eur_n_sazonovs_all_box_plt,eur_n_sazonovs_25_box_plt,nrow = 1)
ggsave(file="/Users/courtneyastore/Dropbox (GaTech)/IBDLocalAncestry/Figures/V3NumberofSazonovsVariants_vs_FractionEuropean.pdf",final_plt,limitsize = FALSE, height = 6, width = 12.5, dpi=300)


#----------------------------------------
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

