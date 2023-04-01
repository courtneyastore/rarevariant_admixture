#!/usr/bin/env python3
# You must run vcftools on the VCF file with regions you would like to assess first:
# vcftools --vcf <your VCF> --out <your out VCF> --site-pi
import argparse
import re
import os
import shutil
import pandas as pd
import subprocess as sp

# Read GRCh38 reference gene coordinates 
def read_gene_pos_file(f):
    d = pd.read_csv(f,sep="\t",low_memory=False)
    df = pd.DataFrame(d)

    print(df)
  
    return df

def get_nuc_div_case(gene_pos_df,d):
    case_file_lst = [i for i in os.listdir(d) if "Case_CDSRegionSazonovs_big_daly_v3" in i and "sites.pi" in i]

    master_df = pd.DataFrame()
    for i in gene_pos_df.index:
        gene = gene_pos_df["SYMBOL"][i]
        chr_ = gene_pos_df["chr"][i]
        start = gene_pos_df["start"][i]
        end = gene_pos_df["end"][i]

        all_sites = list(range(start,end))

        gene_coordinate_df = pd.DataFrame(all_sites)
        gene_coordinate_df.columns = ['POS']
        gene_coordinate_df['CHROM'] = int(chr_)
        gene_coordinate_df['ID'] = gene_coordinate_df['CHROM'].astype(str) + ":" + gene_coordinate_df['POS'].astype(str)

        nucleotide_diversity_file_i = [i for i in case_file_lst if gene in i][0]
        nucleotide_diversity_file_i = d + nucleotide_diversity_file_i
        
        nucleitodie_diversity_gene_i_df = pd.read_csv(nucleotide_diversity_file_i,sep="\t",low_memory=False)
        nucleitodie_diversity_gene_i_df['ID'] = nucleitodie_diversity_gene_i_df['CHROM'].astype(str) + ":" + nucleitodie_diversity_gene_i_df['POS'].astype(str)

        nucleitodie_diversity_gene_i_df = nucleitodie_diversity_gene_i_df[['ID','PI']]
        

        merge_i_df = pd.merge(nucleitodie_diversity_gene_i_df,gene_coordinate_df,on="ID")

        merge_i_df = merge_i_df[['CHROM','POS','PI']]
        master_df = pd.concat([master_df,merge_i_df])

    print(master_df)
    master_df.to_csv("V2_Pi_IBDCase_sazonovs_genes_CDSregions_big_daly_v3_AFR_WGS.sites.pi",sep="\t",index=False)
    return True 


def get_nuc_div_control(gene_pos_df,d):
    case_file_lst = [i for i in os.listdir(d) if "Control_CDSRegionSazonovs_big_daly_v3" in i and "sites.pi" in i]

    master_df = pd.DataFrame()
    for i in gene_pos_df.index:
        gene = gene_pos_df["SYMBOL"][i]
        chr_ = gene_pos_df["chr"][i]
        start = gene_pos_df["start"][i]
        end = gene_pos_df["end"][i]

        all_sites = list(range(start,end))

        gene_coordinate_df = pd.DataFrame(all_sites)
        gene_coordinate_df.columns = ['POS']
        gene_coordinate_df['CHROM'] = int(chr_)
        gene_coordinate_df['ID'] = gene_coordinate_df['CHROM'].astype(str) + ":" + gene_coordinate_df['POS'].astype(str)

        nucleotide_diversity_file_i = [i for i in case_file_lst if gene in i][0]
        nucleotide_diversity_file_i = d + nucleotide_diversity_file_i
        
        nucleitodie_diversity_gene_i_df = pd.read_csv(nucleotide_diversity_file_i,sep="\t",low_memory=False)
        nucleitodie_diversity_gene_i_df['ID'] = nucleitodie_diversity_gene_i_df['CHROM'].astype(str) + ":" + nucleitodie_diversity_gene_i_df['POS'].astype(str)

        nucleitodie_diversity_gene_i_df = nucleitodie_diversity_gene_i_df[['ID','PI']]
        

        merge_i_df = pd.merge(nucleitodie_diversity_gene_i_df,gene_coordinate_df,on="ID")

        merge_i_df = merge_i_df[['CHROM','POS','PI']]
        master_df = pd.concat([master_df,merge_i_df])

    print(master_df)
    master_df.to_csv("V2_Pi_IBDControl_sazonovs_genes_CDSregions_big_daly_v3_AFR_WGS.sites.pi",sep="\t",index=False)
    return True 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help = "Please input tab sep file with gene coordinates")
    parser.add_argument("-d", help = "Please input directory to your nucleotide diversity for keep samples each gene and case control")
    args = parser.parse_args()

    geneStartEnd = read_gene_pos_file(args.f)
    caseNucDiv = get_nuc_div_case(geneStartEnd,args.d)
    controlNucDiv = get_nuc_div_control(geneStartEnd,args.d)

if __name__ == "__main__":
    main()
