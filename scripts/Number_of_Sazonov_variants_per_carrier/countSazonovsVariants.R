library(tidyverse)
#library(dplyr)
library(data.table)
#library("ggpubr")
#library(plyr)


sazonovs_df <- as.data.frame(fread("/storage/home/hcoda1/1/castore3/p-ggibson3-0/sazonov/Formatted_sazonov_significant_variants.txt"))
#sazonovs_df <- sazonovs_df[(sazonovs_df$FigureStatus == "Known causal candidate") | (sazonovs_df$FigureStatus == "New variant in known locus") | (sazonovs_df$FigureStatus == "New Locus"), ]

sazonovs_df$variant_id <- paste(sazonovs_df$Chr,"_",sazonovs_df$Pos,sep="")
sazonovs_df$variant_id <- gsub(" ","",sazonovs_df$variant_id)
sazonovs_df$variant_id <- gsub(",","",sazonovs_df$variant_id)

sazonovs_df <- sazonovs_df[sazonovs_df$variant_id != "16_50729867", ]
sazonovs_df = sazonovs_df[c('variant_id','gene')]

variant_id_map_df <- as.data.frame(fread("/storage/home/hcoda1/1/castore3/p-ggibson3-0/sazonov/MAP_AFR_Emory_sazonov_sig_var.txt"))

sazonovs_variant_id_map_df <- merge(variant_id_map_df,sazonovs_df,by="variant_id")

ukbb_df <- as.data.frame(fread("/storage/home/hcoda1/1/castore3/p-ggibson3-0/get_genotype_for_variants/AFREmory_Sazonovs_GT.tsv"))
ukbb_df <- ukbb_df[c('ID1','ID2','SNP','varID','var_type','label')]
ukbb_df$ID <- paste(ukbb_df$ID1,"_",ukbb_df$ID2,sep="")

ukbb_df <- ukbb_df[c('ID','SNP','var_type','label')]

id_label_df <- ukbb_df[c('ID','label')]
id_label_df <- unique(id_label_df)
#have_var_df <- ukbb_df

have_var_df <- ukbb_df[(ukbb_df$var_type == "het" | ukbb_df$var_type == "hom_alt"), ]

final_df <- merge(have_var_df,sazonovs_variant_id_map_df,by="SNP")

head(final_df)
nrow(final_df)

#variant_id_lst <- as.list(final_df$variant_id)

binarize_df <-  final_df %>% pivot_wider(ID,names_from = variant_id, values_from = variant_id, values_fn = list(variant_id = length),values_fill = list(variant_id = 0))


final_binarize_df <- binarize_df %>% rowwise() %>% mutate(sumrange = rowSums(across(where(is.numeric))))

final_binarize_df <- final_binarize_df %>% arrange(desc(sumrange))

counts_df <- final_binarize_df[c('ID','sumrange')]
counts_df <- merge(counts_df,id_label_df,by="ID")
head(counts_df)

write.table(counts_df,file="V2AllSazonovsVariants_Counts.txt",row.names=FALSE,sep="\t",quote = FALSE)
head(counts_df)
nrow(counts_df)
quit()
write.table(final_binarize_df, file = "V2_BinaryMatrix_AFREmoryIBD_25Sazonovs_significant_variants.tsv", row.names=FALSE, sep="\t",quote = FALSE)
