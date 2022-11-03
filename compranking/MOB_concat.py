#!/usr/bin/env python
# title             :plasMOB_concat.py -> plasMOB_concat.ipynb
# description       :treat plascad result
# author            :Gaoyang Luo
# date              :20221029
# version           :1.0
# usage             :import plasMOB_concat
# required packages :pandas, 
# notification: enjoy yourself
#==============================================================================

import pandas as pd
import os
import path


def plasMOB_concat(input_mob_conj,input_mob_unconj,df_AMR_annotate_contig):
    """
    A plasmid containing relaxase, T4CP, and T4SSs is considered a conjugative plasmid, 
    whereas one encoding only a relaxase is called a mobilizable plasmid.
    Conjugative:T4CP, and T4SSs Mobilizable:relaxase
    
    """
    #load mob ref and AMR annotate_contig
    
    df_mob_conj=pd.read_csv(input_mob_conj,sep="\t",header=None)
    df_mob_unconj=pd.read_csv(input_mob_unconj,sep="\t",header=None)
    df_mob_conj["MOB"]="mob_conj"
    df_mob_unconj["MOB"]="mob_unconj"
    df_mob_conj=pd.concat([df_mob_conj,df_mob_unconj],ignore_index=True)
    
    #MOB merge summary table
    df_AMR_annotate_MOB_contig=pd.merge(df_AMR_annotate_contig,df_mob_conj,left_on="Contig",right_on=0,how="left")
    df_AMR_annotate_MOB_contig=df_AMR_annotate_MOB_contig.drop([0],axis=1,inplace=False)
    df_AMR_annotate_MOB_contig=df_AMR_annotate_MOB_contig.fillna("-")
    
    #according to MOB to filter phage again
    a=df_AMR_annotate_MOB_contig[df_AMR_annotate_MOB_contig.CompRanking_MGE_prediction == "phage"]
    b=a[a.MOB == "mob_conj"]
    index_list=[]
    for i, name in b.iterrows():
        index_list.append(i)

    for i in index_list:
        df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="ambiguous (conj/phage)"
    
    
    #according to MOB to filter ambiguous (plasmid/phage), if MOB == conj, ambiguous = plasmid
    a=df_AMR_annotate_MOB_contig[df_AMR_annotate_MOB_contig.CompRanking_MGE_prediction == "ambiguous (plasmid/phage)"]
    b=a[a.MOB == "mob_conj"]
    index_list=[]
    for i, name in b.iterrows():
        index_list.append(i)

    for i in index_list:
        df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="plasmid"
    
    #according to MobileOG to filter ambiguous (plasmid/phage) again to generate final phage
    a=df_AMR_annotate_MOB_contig[df_AMR_annotate_MOB_contig.CompRanking_MGE_prediction == "ambiguous (plasmid/phage)"]
    b=a[a.Taxonomy == "phage"]
    index_list=[]
    for i, name in b.iterrows():
        index_list.append(i)

    for i in index_list:
        df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="phage"  
        
        
    ####save####
    # df_AMR_annotate_contig.to_csv("/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/CompRanking_ERR1191817_AMR_result.csv",sep="\t",index=0)
    # AMRfile=df_AMR_annotate_MOB_contig.to_csv(output + "/CompRanking_" + filebase + "_AMR_MOB_prediction.tsv", sep="\t", index=0)
    
    return df_AMR_annotate_MOB_contig  
    
    
    
    
    
    
    
    if __name__ == "__main__":
        input_mob_conj="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/MGE/plascad/ERR1191817.contigs_5M_contigs_Conj_plasmids_id_out"
        input_mob_unconj="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/MGE/plascad/ERR1191817.contigs_5M_contigs_mob_unconj_plasmids_id_out"
        df_AMR_annotate_contig="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_result/CompRanking_ERR1191817.contigs_AMR_prediction.tsv"
        
        
        
        