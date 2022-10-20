#!/usr/bin/env python
# title             :ARG_ranker.py
# description       :function ARG ranking from SARG
# author            :Gaoyang Luo
# date              :202209020
# version           :1.0
# usage             :import AMRcombined
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================

import pandas as pd


def arg_rank(self, input_sarg, ):
    #coveragy cal
    df_sarg=pd.read_csv(input_sarg, sep="\t",header=None)
    df_sarg.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']

    #open df_sarg
    df_SARG_len=pd.read_csv("/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.db.fasta.length", sep="\t", header=None )
    df_SARG_len=df_SARG_len[[1,2]]

    #merge len and df_sarg
    df_sarg_len=pd.merge(df_sarg,df_SARG_len,left_on="sub_id",right_on=1,how="left")
    df_sarg_len["coverage"]="-"

    #cover cal
    for index, name in df_sarg_len.iterrows():
        df_sarg_len["coverage"][index]=df_sarg_len["alignLen"][index]/df_sarg_len[2][index]

    # filter out contigs identity under 
    df_sarg_iden95 = df_sarg_len[df_sarg_len.identity > 90]
    # filter out contigs alignment length under 25
    df_sarg_alen = df_sarg_iden95[df_sarg_iden95.coverage>0.8]
    df_sarg_tmp=df_sarg_alen.iloc[:,[0,1]]
    df_sarg_tmp.columns=["Contig_ID", "query"]
    
    
    df_sarg_structure=pd.read_csv(input_sarg_structure, sep="\t", header=0)

        
#arg rank processing
input_sarg="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/AMR/ARGranking/ERR1191817.contigs_5M_contigs_SARG_Protein_diamond.txt"
input_sarg_structure="/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.structure.txt"



    
    