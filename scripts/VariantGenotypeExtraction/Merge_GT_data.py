#!/usr/bin/env python3
# python3 Step1_extract_genes_from_multisample_VCF.py -f ../example_gene_file.lst -d /storage/home/hcoda1/1/castore3/scratch/Emory_269_all/
import argparse
import re
import os
import pandas as pd
import subprocess as sp


def assess_genotypes(vcf_chr_dir,map_file):
    chr_dir_lst = os.listdir(vcf_chr_dir)
    chr_dir_lst = [i for i in chr_dir_lst if ".ped" in i]

    map_d = pd.read_csv(map_file,sep="\t",names=['varID','SNP'])
    map_df = pd.DataFrame(map_d)

    #map_df['variant_id'] = map_df['Location'].astype(str).str.replace(":","-") + "-" + map_df['REF_ALLELE'].astype(str) + "-" + map_df['Allele'].astype(str)
    #map_df = map_df[['variant_id','#Uploaded_variation']]
    #map_df = map_df.drop_duplicates()
    
    sazonov = pd.read_csv("/storage/home/hcoda1/1/castore3/p-ggibson3-0/sazonov/sazonov_variant_id_info.txt",sep=",")
    sazonov_df = pd.DataFrame(sazonov)
    sazonov_df['varID'] = sazonov_df['Chr'].astype(str) + "_" + sazonov_df['Pos'].astype(str)
    sazonov_df.columns = ['chr','pos','A1','A2','varID']

    #sazonov_df['variant_id'] =  sazonov_df['Chr'].astype(str) + "-" + sazonov_df['Pos'].astype(str) + "-" + sazonov_df['A0'].astype(str) + "-" + sazonov_df['A1'].astype(str)

    sazonov_map_df = pd.merge(sazonov_df,map_df,on="varID",how="left")
    
    print(sazonov_map_df)

    # Loop over folders in dir for each chromosome
    # 'ID1', 'ID2', 'fam1', 'fam2', 'fam3', 'fam4', 'GT_A1', 'GT_A2', 'GT','SNP', 'chr', 'pos', 'A1', 'A2', 'varID', 'var_type', 'label'
    master_df = pd.DataFrame(columns=['ID1', 'ID2', 'fam1', 'fam2', 'fam3', 'fam4', 'GT_A1', 'GT_A2', 'GT','SNP', 'chr', 'pos', 'A1', 'A2', 'varID', 'var_type', 'label']) 
    for i in chr_dir_lst:
        variant_id = str(i).replace(".txt.ped","")
        i_file = vcf_chr_dir + str(i)
        d = pd.read_csv(i_file,sep=" ",names=['ID1','ID2','fam1','fam2','fam3','fam4','GT_A1','GT_A2'])
        df = pd.DataFrame(d)
        varID = str(i).replace(".txt.ped","").replace("AFR_GT_","")
        df['GT'] = df['GT_A1'] + df['GT_A2']
        df['SNP'] = varID

        SNP_df = pd.merge(df,sazonov_map_df,on="SNP")
        print(SNP_df)
        het_df = SNP_df[SNP_df['GT_A1'] != SNP_df['GT_A2']]
        het_df['var_type'] = "het"
        ref_hom_df = SNP_df[(SNP_df['GT_A1'] == SNP_df['GT_A2']) & (SNP_df['GT_A1'] == SNP_df['A1'])]
        ref_hom_df['var_type'] = "hom_ref"
        alt_hom_df = SNP_df[(SNP_df['GT_A1'] == SNP_df['GT_A2']) & (SNP_df['GT_A1'] == SNP_df['A2'])]
        alt_hom_df['var_type'] = "hom_alt"

        final_df = pd.concat([het_df,ref_hom_df,alt_hom_df])
        
        case_df = final_df[final_df['fam4'] == 2]
        case_df['label'] = "Case"

        control_df = final_df[final_df['fam4'] == 1]
        control_df['label'] = "Control"

        case_control_final_df = pd.concat([case_df,control_df])

        master_df = pd.concat([master_df,case_control_final_df])
        
    print(master_df)
    master_df.to_csv("AFREmory_Sazonovs_GT.tsv",sep="\t",index=False)
        
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", help = "Please input directory to your individual chr directories")
    parser.add_argument("-v", help = "Please input VEP output")
    args = parser.parse_args()

    assessGenotypes = assess_genotypes(args.d,args.v)
    
if __name__ == "__main__":
    main()
