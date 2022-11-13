#!/usr/bin/env python
# title             :summary.py
# description       :combine AMR and VF and PATH results
# author            :Gaoyang Luo
# date              :20221113
# version           :1.0
# usage             :import summary
# required packages :pandas,
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
# from compranking import path

def sum_all(df_AMR_MOB, df_VF_PATH):
    """
    first step: sum all the result
    second step: give special risk rank of ARGs by CompRanking criterion
    """
    #first step
    df_sum=pd.merge(df_AMR_MOB,df_VF_PATH,left_on="ORF_ID",right_on="ORF",how="left")
    df_sum=df_sum.drop(["Contig_y","ORF","VF_Host_Bacteria","coverage","SNPs","mobileOG_ID","VF_Name"],axis=1,inplace=False)
    df_sum.columns=["Contig","ORF_ID","ARG_prediction","ARG_class","ARG_Database","ARG_rank","MGE_prediction","MGE_Name","MGE_Taxonomy","MGE_Category","MGE_Database","Mobility","Viculence_factor","VF_Category","VF_ESKAPE","Pathogenicity"]
    #second step
    cpr_rank=[]
    for i, name in df_sum.iterrows():
        if name["ARG_prediction"] == "-":
            cpr_rank.append("IV")
        elif name["ARG_prediction"]!= "-" and name["MGE_Name"] == "-":
            cpr_rank.append("III")
        elif name["ARG_prediction"]!= "-" and name["MGE_Name"] != "-" and name["Pathogenicity"] == "-":
            cpr_rank.append("II")
        elif name["ARG_prediction"]!= "-" and name["MGE_Name"] != "-" and name["Pathogenicity"] != "-":
            cpr_rank.append("I")
        else:
            cpr_rank.append("-")
    
    return df_sum
    