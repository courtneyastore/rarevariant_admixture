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

sazonov_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/sazonovs/sazonov_significant_variants.txt"))
sazonov_df <- sazonov_df[(sazonov_df$FigureStatus == "Known causal candidate") | (sazonov_df$FigureStatus == "New variant in known locus") | (sazonov_df$FigureStatus == "New Locus"), ]

other_af_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/sazonovs/Sazonovs_CausalVariants_GNOMAD_AFRSubPop.txt"))

ukbb_ibd_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/AFs/UKBB_IBD.frq"))
ukbb_control_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/AFs/UKBB_Control.frq"))

afr_emory_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/processing_files/AFs/mapID_format_AFR_Sazonovs_a2_allele_supplied_AFs.frq.cc"))

ukbb_ibd_df$variant_id <- gsub(":","-",ukbb_ibd_df$SNP)
ukbb_control_df$variant_id <- gsub(":","-",ukbb_control_df$SNP)

my_cols <- c("Chr", "Pos", "A0","A1")    
sazonov_df$variant_id <- do.call(paste, c(sazonov_df[my_cols], sep = "-"))
sazonov_df$variant_id <- gsub(" ","",sazonov_df$variant_id)
sazonov_df$variant_id <- gsub(",","",sazonov_df$variant_id)

gnomad_df <- other_af_df[c("varID","gnomad_AFR_AF","gnomad_NFE_AF")]
colnames(gnomad_df) <- c("variant_id","gnomad_AFR_AF","gnomad_NFE_AF")

risk_sazonov_df <- sazonov_df[sazonov_df$or_meta > 1, ]
risk_sazonov_df$risk_protective <- "Risk"

protective_sazonov_df <- sazonov_df[sazonov_df$or_meta < 1, ]
protective_sazonov_df$risk_protective <- "Protective"

common_sazonov_df <- sazonov_df[sazonov_df$gnomad_NFE_AF >= 0.05, ]
common_sazonov_df$common_rare<- "Common"

low_freq_sazonov_df <- sazonov_df[(sazonov_df$gnomad_NFE_AF < 0.05) & (sazonov_df$gnomad_NFE_AF >= 0.005), ]
low_freq_sazonov_df$common_rare <- "Low frequency"

rare_freq_sazonov_df <- sazonov_df[(sazonov_df$gnomad_NFE_AF < 0.005) & (sazonov_df$gnomad_NFE_AF >= 0.001), ]
rare_freq_sazonov_df$common_rare <- "Rare"

ultra_rare_freq_sazonov_df <- sazonov_df[(sazonov_df$gnomad_NFE_AF < 0.001), ]
ultra_rare_freq_sazonov_df$common_rare <- "Ultra rare"

common_rare_sazonov_df <- rbind(common_sazonov_df,low_freq_sazonov_df,rare_freq_sazonov_df,ultra_rare_freq_sazonov_df)
common_rare_sazonov_df <- common_rare_sazonov_df[c('variant_id','common_rare')]

protective_risk_sazonov_df <- rbind(risk_sazonov_df,protective_sazonov_df)
protective_risk_sazonov_df <- protective_risk_sazonov_df[c("variant_id","consequence","gene","FigureStatus","risk_protective")]

final_sazonov_df <- merge(protective_risk_sazonov_df,common_rare_sazonov_df,by="variant_id")


gnomad_nfe_af_df <- gnomad_df[c("variant_id",'gnomad_NFE_AF')]
colnames(gnomad_nfe_af_df) <- c("variant_id","AF")
gnomad_nfe_af_df$label <- "GNOMAD NFE"
gnomad_nfe_af_df$color <- "springgreen4"

gnomad_afr_af_df <- gnomad_df[c("variant_id",'gnomad_AFR_AF')]
colnames(gnomad_afr_af_df) <- c("variant_id","AF")
gnomad_afr_af_df$label <- "GNOMAD AFR"
gnomad_afr_af_df$color <- "skyblue"

ukbb_ibd_af_df <- ukbb_ibd_df[c('variant_id','MAF')]
colnames(ukbb_ibd_af_df) <- c('variant_id','AF')
ukbb_ibd_af_df$label <- "UKBB IBD"
ukbb_ibd_af_df$color <- "palevioletred4"

ukbb_control_af_df <- ukbb_control_df[c('variant_id','MAF')]
colnames(ukbb_control_af_df) <- c('variant_id','AF')
ukbb_control_af_df$label <- "UKBB Control"
ukbb_control_af_df$color <- "palevioletred3"

nextera_ibd_af_df <- sazonov_df[c("variant_id",'Case_AF1(Nextera.SAIGE)')]
colnames(nextera_ibd_af_df) <- c("variant_id","AF")
nextera_ibd_af_df$label <- "Nextera IBD (Sazonovs et.al.)"
nextera_ibd_af_df$color <- "salmon4"

nextera_control_af_df <- sazonov_df[c("variant_id",'Control_AF1(Nextera.SAIGE)')]
colnames(nextera_control_af_df) <- c("variant_id","AF")
nextera_control_af_df$label <- "Nextera Control (Sazonovs et.al.)"
nextera_control_af_df$color <- "sandybrown"

twist_ibd_af_df <- sazonov_df[c("variant_id",'Case_AF1(TWIST.SAIGE)')]
colnames(twist_ibd_af_df) <- c("variant_id","AF")
twist_ibd_af_df$label <- "TWIST IBD (Sazonovs et.al.)"
twist_ibd_af_df$color <- "orangered3"

twist_control_af_df <- sazonov_df[c("variant_id",'Control_AF1(TWIST.SAIGE)')]
colnames(twist_control_af_df) <- c("variant_id","AF")
twist_control_af_df$label <- "TWIST Control (Sazonovs et.al.)"
twist_control_af_df$color <- "#lightcoral"

# African IBD emory data
afr_emory_IBD_af_df <- afr_emory_df[c("variant_id",'AF_AFR_Emory_het_case')]
colnames(afr_emory_IBD_af_df) <- c("variant_id","AF")
afr_emory_IBD_af_df$label <- "AFR Emory IBD"
afr_emory_IBD_af_df$color <- "blue4"
# African control emory data
afr_emory_control_af_df <- afr_emory_df[c("variant_id",'AF_AFR_Emory_het_control')]
colnames(afr_emory_control_af_df) <- c("variant_id","AF")
afr_emory_control_af_df$label <- "AFR Emory Control"
afr_emory_control_af_df$color <- "royalblue1"

final_df <- rbind(gnomad_nfe_af_df,gnomad_afr_af_df,nextera_ibd_af_df,nextera_control_af_df,twist_ibd_af_df,twist_control_af_df,ukbb_ibd_af_df,ukbb_control_af_df,gnomad_afr_af_df,afr_emory_IBD_af_df,afr_emory_control_af_df)

final_df1 <- merge(final_df,final_sazonov_df,on="variant_id")

final_df1 <- 

colors  <- c("AFR Emory IBD" = "blue4",
             "AFR Emory Control" = "royalblue1",
             "GNOMAD AFR" = "skyblue",
             "UKBB IBD" = "palevioletred4",
             "UKBB Control" = "palevioletred3",
             "TWIST IBD (Sazonovs et.al.)" = "orangered3",
             "TWIST Control (Sazonovs et.al.)" = "lightcoral",
             "Nextera IBD (Sazonovs et.al.)" = "salmon4",
             "Nextera Control (Sazonovs et.al.)" = "sandybrown",
             "GNOMAD NFE" = "springgreen4")

final_df1$label2 <- factor(final_df1$label,levels = c("AFR Emory IBD", "AFR Emory Control", "GNOMAD AFR", "UKBB IBD","UKBB Control","TWIST IBD (Sazonovs et.al.)","TWIST Control (Sazonovs et.al.)","Nextera IBD (Sazonovs et.al.)","Nextera Control (Sazonovs et.al.)","GNOMAD NFE"))
final_df1 <- final_df1[complete.cases(final_df1), ]
final_df1$AF <- as.numeric(final_df1$AF)
final_df1 <- final_df1[final_df1$AF != 0, ]

bar_plt <- ggplot(final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(title="AF of Sazonovs et.al. 25 significant variants (known causal, new variants in known locus, and new locus) across IBD & control cohorts", x="Variant IDs") +
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5),legend.position="bottom",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")

leg <- get_legend(bar_plt)
leg_plt <- as_ggplot(leg)
leg_plt

common_nod2_final_df1 <- final_df1[(final_df1$gene == "NOD2") & (final_df1$common_rare == "Low frequency") & (final_df1$variant_id != "16-50712018-C-T"), ] 
common_nod2_final_df1$gene <- "NOD2 common variants"
nod2_common_bar_plt <- ggplot(common_nod2_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
nod2_common_bar_plt

rare_nod2_final_df1 <- final_df1[(final_df1$gene == "NOD2") & (final_df1$common_rare == "Rare" | final_df1$common_rare == "Ultra rare" | final_df1$variant_id == "16-50712018-C-T") , ] 
rare_nod2_final_df1$gene <- "NOD2 rare variants"
nod2_rare_bar_plt <- ggplot(rare_nod2_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
nod2_rare_bar_plt

hgfac_final_df1 <- final_df1[final_df1$gene == "HGFAC", ]
hgfac_bar_plt <- ggplot(hgfac_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.1) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
hgfac_bar_plt

slc39A8_final_df1 <- final_df1[final_df1$gene == "SLC39A8", ]
slc39A8_bar_plt <- ggplot(slc39A8_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.1) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
slc39A8_bar_plt

ccr7_final_df1 <- final_df1[final_df1$gene == "CCR7", ]
ccr7_bar_plt <- ggplot(ccr7_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.06) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
ccr7_bar_plt

dok2_final_df1 <- final_df1[final_df1$gene == "DOK2", ]
dok2_bar_plt <- ggplot(dok2_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.06) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
dok2_bar_plt

sdf2l1_final_df1 <- final_df1[final_df1$gene == "SDF2L1", ]
sdf2l1_bar_plt <- ggplot(sdf2l1_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.06) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
sdf2l1_bar_plt

rela_final_df1 <- final_df1[final_df1$gene == "RELA", ]
rela_bar_plt <- ggplot(rela_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.01) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
rela_bar_plt

ptafr_final_df1 <- final_df1[final_df1$gene == "PTAFR", ]
ptafr_bar_plt <- ggplot(ptafr_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.01) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
ptafr_bar_plt

pdlim5_final_df1 <- final_df1[final_df1$gene == "PDLIM5", ]
pdlim5_bar_plt <- ggplot(pdlim5_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.01) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
pdlim5_bar_plt

il10ra_final_df1 <- final_df1[final_df1$gene == "IL10RA", ]
il10ra_bar_plt <- ggplot(il10ra_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + theme_classic() + ylim(0,0.01) + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
il10ra_bar_plt

# Protective variants 
il23r_final_df1 <- final_df1[final_df1$gene == "IL23R", ]
il23r_bar_plt <- ggplot(il23r_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.07) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
il23r_bar_plt

tyk2_final_df1 <- final_df1[final_df1$gene == "TYK2", ]
tyk2_bar_plt <- ggplot(tyk2_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.07) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
tyk2_bar_plt

tagap_final_df1 <- final_df1[final_df1$gene == "TAGAP", ]
tagap_bar_plt <- ggplot(tagap_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.03) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
tagap_bar_plt

card9_final_df1 <- final_df1[final_df1$gene == "CARD9", ]
card9_bar_plt <- ggplot(card9_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="",y="") + ylim(0,0.008) + theme_classic() + 
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
card9_bar_plt

# Risk
risk_plt1 <- ggarrange(hgfac_bar_plt, slc39A8_bar_plt,nod2_common_bar_plt, nod2_rare_bar_plt,ccr7_bar_plt,dok2_bar_plt,sdf2l1_bar_plt,ncol=4,nrow=2)
risk_plt1
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/V2Figure1_risk1.pdf",limitsize = FALSE, height = 10, width = 18, dpi=300)

risk_plt2 <- ggarrange(rela_bar_plt, ptafr_bar_plt,pdlim5_bar_plt, il10ra_bar_plt,ncol=4,nrow=1)
risk_plt2
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/V2Figure1_risk2.pdf",limitsize = FALSE, height = 5, width = 18, dpi=300)

# Protective
protective_plot <- ggarrange(il23r_bar_plt,tyk2_bar_plt, tagap_bar_plt, card9_bar_plt,nrow=1,ncol = 4)
protective_plot
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/V2Figure1_protective.pdf",limitsize = FALSE, height = 5, width = 18, dpi=300)





final_plt <- ggarrange(risk_plt1,risk_plt2,protective_plot,leg_plt,nrow=4,ncol=1)
final_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1_V4.pdf",limitsize = FALSE, height = 30, width = 18, dpi=300)




risk_final_df1 <- final_df1[(final_df1$risk_protective == "Risk") & (final_df1$gene != "NOD2"), ]

# Figure 1A
common_risk_final_df1 <- risk_final_df1[risk_final_df1$common_rare == "Common", ]
risk_common_bar_plt <- ggplot(common_risk_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="Variant IDs",y="AF",title="Common risk variants") +
  scale_fill_manual(name="Cohort",values = colors) + ylim(0,0.1) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
risk_common_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1A_commonrisk_AFSazonovs.pdf",limitsize = FALSE, height = 6, width = 10, dpi=300)

# Figure 1B
rare_risk_final_df1 <- risk_final_df1[risk_final_df1$common_rare != "Common", ]
risk_rare_bar_plt <- ggplot(rare_risk_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="Variant IDs",y="AF",title="Rare risk variants") +
  scale_fill_manual(name="Cohort",values = colors) + ylim(0,0.06) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
risk_rare_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1B_rarerisk_AFSazonovs.pdf",limitsize = FALSE, height = 14, width = 12, dpi=300)

# Figure 1C
nod2_risk_final_df1 <- final_df1[(final_df1$risk_protective == "Risk") & (final_df1$gene == "NOD2"), ]
nod2_risk_bar_plt <- ggplot(nod2_risk_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="Variant IDs",y="AF",title="NOD2 risk variants") +
  scale_fill_manual(name="Cohort",values = colors) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~FigureStatus,scale="free")
nod2_risk_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1C_NOD2risk_AFSazonovs.pdf",limitsize = FALSE, height = 6, width = 10, dpi=300)

# Figure 1D
protective_final_df1 <- final_df1[final_df1$risk_protective == "Protective", ]
protective_bar_plt <- ggplot(protective_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="Variant IDs",y="AF",title="Protective variants") +
  scale_fill_manual(name="Cohort",values = colors) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45,hjust=1)) +
  facet_wrap(~gene,scale="free",nrow=1)
protective_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1D_protective_AFSazonovs.pdf",limitsize = FALSE, height = 6, width = 16, dpi=300)


leg <- get_legend(protective_bar_plt)
leg_plt <- as_ggplot(leg)
leg_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1_legend.pdf",height = 8, width = 18,dpi = 300)




#---------------
risk_bar_plt <- ggplot(risk_final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x="Variant IDs",y="AF") +
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5),legend.position="bottom",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
risk_bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1A_risk_AFSazonovs.pdf",limitsize = FALSE, height = 15, width = 15, dpi=300)

bar_plt <- ggplot(final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(title="AF of Sazonovs et.al. 25 significant variants (known causal, new variants in known locus, and new locus) across IBD & control cohorts", x="Variant IDs") +
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5),legend.position="bottom",axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_wrap(~gene,scale="free")
bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1_all_AFSazonovs.pdf",limitsize = FALSE, height = 35, width = 42, dpi=300)



bar_plt <- ggplot(final_df1, aes(x = as.factor(variant_id), y = AF,fill=as.factor(label2))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(title="AF of Sazonovs et.al. 25 significant variants (known causal, new variants in known locus, and new locus) across IBD & control cohorts", x="Variant IDs") +
  scale_fill_manual(name="Cohort",values = colors) + 
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5),legend.position="bottom",axis.text.x = element_text(angle = 45,hjust=1)) + 
   facet_wrap(~gene,scale="free")
bar_plt
ggsave("/Users/courtneyastore/Desktop/final_IBD_AFR_EUR/figures/Figure1_AFSazonovs.pdf",limitsize = FALSE, height = 35, width = 42, dpi=300)
