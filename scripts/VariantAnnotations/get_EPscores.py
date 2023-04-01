#!/usr/bin/env python3

import argparse
import re
import os
import pandas as pd
import subprocess as sp

def read_input_var_file(f):
    d = pd.read_csv(f,sep="\t",names=['np','aa_change'])
    df = pd.DataFrame(d)
    print(df)

    return df
    
def read_sazonovs_file():
    d = pd.read_csv("/storage/home/hcoda1/1/castore3/p-ggibson3-0/sazonov/sazonovs_sigvar_aminoacid_ep_map.txt",sep="\t")
    df = pd.DataFrame(d)
    print(df)

    return df


def annotate_ep_scores(protein_mod_df,sazonovs_map_df,out_file_name):

    d = pd.read_csv("/storage/home/hcoda1/1/castore3/p-ggibson3-0/evolutionary_probability/skumar_ep_matrix/all_protein_eps.csv",sep=",")
    ep_df = pd.DataFrame(d)
    
    np_pos_ep_d = []
    for i,j in zip(protein_mod_df['np'], protein_mod_df['aa_change']):
        np_i = str(i).strip(" ")
        ref_aa = str(j)[0].lower()
        new_aa = str(j)[-1].lower()
        pos = int(str(j)[1:-1].strip(" "))

        if (ref_aa == "x") or (new_aa == "x"):
            ep_i = 'NA'
        else:
            prot_mod = np_i + ":" + ref_aa + str(pos) + new_aa
            
            ep_df_i = ep_df[(ep_df['protcore'] == np_i) & (ep_df['aa_pos'] == pos)]

            if len(ep_df_i) == 0:
                ep_i = 'NA'
            else:
                ep_i = ep_df_i[new_aa].iloc[0]
        
        print("EP-score for",np_i,j,"is:",ep_i)
        np_pos_ep_d.append(
            {
                'NP': np_i,
                'AminoAcidSub': j,
                'ep': ep_i
            }
        )
    protein_mod_ep_df = pd.DataFrame(np_pos_ep_d)

    print(protein_mod_ep_df)

    final_df = pd.merge(sazonovs_map_df,protein_mod_ep_df,on="AminoAcidSub")

    final_df['Pos'] = final_df['Pos'].str.replace(",","")
    print(final_df)

    final_df.to_csv(out_file_name,sep="\t",index=False,header=True)

    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help = "Please input directory to your individual chr directories")
    args = parser.parse_args()

    out_file_name = "EPscore_" + os.path.basename(str(args.f))

    inputVar = read_input_var_file(args.f)
    sazonovsMap = read_sazonovs_file()

    annotateEP = annotate_ep_scores(inputVar,sazonovsMap,out_file_name)
 
if __name__ == "__main__":
    main()
