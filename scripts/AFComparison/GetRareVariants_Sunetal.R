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

# Add disease group map from ST2
sun_et_al_diseasemap_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/Sun_etal_ST2_disease_summary.txt"))
sun_et_al_diseasemap_df <- sun_et_al_diseasemap_df[c('Disease cluster ID','Disease group')]
sun_et_al_diseasemap_df <- unique(sun_et_al_diseasemap_df)

# Add CWAS associations from ST3
sun_et_al_sigvar_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/Sun_etal_ST3_CWAS_associations.txt"))

# Merge CWAS with disease group information
merge_diseasegroup_sun_et_al_sigvar_df <- merge(sun_et_al_sigvar_df,sun_et_al_diseasemap_df,by="Disease cluster ID")
drop <- c("Disease cluster ID")
merge_diseasegroup_sun_et_al_sigvar_df = merge_diseasegroup_sun_et_al_sigvar_df[,!(names(merge_diseasegroup_sun_et_al_sigvar_df) %in% drop)]
merge_diseasegroup_sun_et_al_sigvar_df <- unique(merge_diseasegroup_sun_et_al_sigvar_df)

#---------------------Get numbers-------------------------------
# Number of unique significant associations
nrow(sun_et_al_sigvar_df)

# Number of unique disease groups
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$`Disease group`))

# Number of unique diseases finn gen and ukbb 
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$`Disease description (UKB)`))
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$`Disease description (FG)`))

# Number of unique variants
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$`Variant ID`))
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$rsID))

# Number of unique variant consequences
length(unique(merge_diseasegroup_sun_et_al_sigvar_df$`Most severe consequence`))

#-----------------------------------------------------------------

# Get MAF from MAF % for UKB and FG
sigvar_df <- sun_et_al_sigvar_df[c('Variant ID','Disease description (UKB)','Disease description (FG)','MAF (UKB)','MAF (FG)','Disease cluster ID')]
sigvar_df <- unique(sigvar_df)
colnames(sigvar_df) <- c("varID",'UKB_disease','FinnGen_disease','MAF_UKBB','MAF_FG','Disease cluster ID')
sigvar_df$MAF_UKBB <- as.numeric(gsub("%","",sigvar_df$MAF_UKBB)) / 100
sigvar_df$MAF_FG <- as.numeric(gsub("%","",sigvar_df$MAF_FG)) / 100

# Number of UKB rare associations MAF <0.01 
nrow(unique(sigvar_df[sigvar_df$MAF_UKBB < 0.01, ]))
nrow(unique(sigvar_df[sigvar_df$MAF_FG < 0.01, ]))
nrow(unique(sigvar_df[(sigvar_df$MAF_FG < 0.01) | (sigvar_df$MAF_UKBB < 0.01), ]))

rare_sigvar_df <- unique(sigvar_df[(sigvar_df$MAF_FG < 0.01) | (sigvar_df$MAF_UKBB < 0.01), ])

#-------------------------------------------------------------------

# Get VEP AFs

#vep_df <- as.data.frame(fread("/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/format_VEP_output_Sun_etal_ST3_CWAS_associations.txt.vcf"))
vep_df1 <- as.data.frame(fread("/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/format_v2VEP_output_Sun_etal_ST3_CWAS_associations.txt.vcf"))
vep_df2 <- as.data.frame(fread("/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/format_vep_output_missing_Sunetal_variants.vcf"))

vep_df1 <- vep_df1[c('#Uploaded_variation','REF_ALLELE','Allele','gnomADg_NFE_AF','gnomADg_AFR_AF')]
vep_df1$varID <- paste(vep_df1$`#Uploaded_variation`,vep_df1$REF_ALLELE,vep_df1$Allele,sep=":")
vep_df1 <- vep_df1[c('varID','REF_ALLELE','Allele','gnomADg_NFE_AF','gnomADg_AFR_AF')]
vep_df1 <- unique(vep_df1)

vep_df2 <- vep_df2[c('#Uploaded_variation','REF_ALLELE','Allele','gnomADg_NFE_AF','gnomADg_AFR_AF')]
colnames(vep_df2) <- c("varID",'REF_ALLELE',"Allele",'gnomADg_NFE_AF','gnomADg_AFR_AF')
#vep_df2$varID <- paste(vep_df2$`#Uploaded_variation`,vep_df2$REF_ALLELE,vep_df2$Allele,sep=":")
vep_df2 <- unique(vep_df2)

vep_df <- rbind(vep_df1,vep_df2)
varID_sigvar_df <- sigvar_df[c("varID")]
varID_sigvar_df <- unique(varID_sigvar_df)
final_vep_df <- merge(varID_sigvar_df,vep_df,by="varID")

# Number of variants with GNOMAD NFE AF > 0.5 
nrow(final_vep_df[final_vep_df$gnomADg_NFE_AF > 0.5, ])

# Correct alleles
major_allele_vep_df <- final_vep_df[final_vep_df$gnomADg_NFE_AF > 0.5, ]
major_allele_vep_df$gnomADg_NFE_AF <- 1-major_allele_vep_df$gnomADg_NFE_AF
major_allele_vep_df$gnomADg_AFR_AF <- 1-major_allele_vep_df$gnomADg_AFR_AF
major_allele_vep_df$corrected_varID <- paste(major_allele_vep_df$`#Uploaded_variation`,major_allele_vep_df$Allele,major_allele_vep_df$REF_ALLELE,sep=":")
#for_mapping_major_allele_vep_df <- major_allele_vep_df

major_allele_vep_df <- major_allele_vep_df[c('varID','corrected_varID','gnomADg_NFE_AF','gnomADg_AFR_AF')]

minor_allele_vep_df <- final_vep_df[final_vep_df$gnomADg_NFE_AF <= 0.5, ]
minor_allele_vep_df <- minor_allele_vep_df[c('varID','gnomADg_NFE_AF','gnomADg_AFR_AF')]
minor_allele_vep_df$corrected_varID <- minor_allele_vep_df$varID
#for_mapping_minor_allele_vep_df <- minor_allele_vep_df

corrected_alleles_df <- rbind(minor_allele_vep_df,major_allele_vep_df)

# Number of rare variants based on GNOMAD NDW
nrow(corrected_alleles_df[corrected_alleles_df$gnomADg_NFE_AF < 0.01, ])

# Merge corrected allele information with significant associations
final_df <- merge(sigvar_df,corrected_alleles_df,by="varID")

# Number of significant disease-rare variant associations
nrow(unique(final_df[final_df$gnomADg_NFE_AF < 0.01, ]))
rare_variant_disease_sigvar_df <- unique(final_df[final_df$gnomADg_NFE_AF < 0.01, ])

# Number of unique disease groups
diseaseinfo_rare_variant_disease_sigvar_df <- merge(rare_variant_disease_sigvar_df,sun_et_al_diseasemap_df,by="Disease cluster ID",all.x=TRUE)
drop <- c("Disease cluster ID")
diseaseinfo_rare_variant_disease_sigvar_df = diseaseinfo_rare_variant_disease_sigvar_df[,!(names(diseaseinfo_rare_variant_disease_sigvar_df) %in% drop)]
diseaseinfo_rare_variant_disease_sigvar_df <- unique(diseaseinfo_rare_variant_disease_sigvar_df)
length(unique(diseaseinfo_rare_variant_disease_sigvar_df$`Disease group`))

# Export final corrected allele, rare (MAF < 0.01) variant disease associations. 
write.table(diseaseinfo_rare_variant_disease_sigvar_df,"/Users/courtneyastore/Desktop/Sun_etal_analysis/processing_files/Sun_etal_corrected_alleles_VEP_sigvardisease_variant_associations_GNOMADgNFEAF0.01.tsv",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

