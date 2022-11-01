#!/usr/bin/env python
# title             :ARG_ranker.py -> SARG_processing.ipynb
# description       :function ARG ranking from SARG
# author            :Gaoyang Luo
# date              :202209020
# version           :1.0
# usage             :import AMRcombined
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================

import pandas as pd
import os
import path


def arg_rank(input_sarg, input_sarg_length,input_sarg_structure, input_argrank,filebase, output):
    #coveragy cal
    df_sarg=pd.read_csv(input_sarg, sep="\t",header=None)
    df_sarg.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']

    #open df_sarg
    df_SARG_len=pd.read_csv(input_sarg_length, sep="\t", header=None )
    df_SARG_len=df_SARG_len[[1,2]]

    #merge len and df_sarg
    df_sarg_len=pd.merge(df_sarg,df_SARG_len,left_on="sub_id",right_on=1,how="left")
    df_sarg_len["coverage"]="-"

    #cover cal
    for index, name in df_sarg_len.iterrows():
        df_sarg_len["coverage"][index]=df_sarg_len["alignLen"][index]/df_sarg_len[2][index]

    # filter out contigs identity under 60
    df_sarg_iden80 = df_sarg_len[df_sarg_len.identity > 60]
    # filter out contigs coverage more than 0.7
    df_sarg_cov = df_sarg_iden80[df_sarg_iden80.coverage>0.7]
    df_sarg_tmp=df_sarg_cov.iloc[:,[0,1]]
    df_sarg_tmp.columns=["Contig_ID", "query"]
    
    #load SARG structure file
    df_sarg_structure=pd.read_csv(input_sarg_structure, sep="\t", header=0)
    
    #merge structure and filtered file
    df_sarg_type=pd.merge(df_sarg_tmp,df_sarg_structure,left_on="query",right_on="Reference",how="left")
    df_sarg_type1=df_sarg_type.drop(['Reference'], axis=1, inplace=False)
    
    #standardizing format
    df_sarg_type1["class"]=df_sarg_type1["Genotype"].str.split("__",expand=True)[1]
    df_sarg_type2=df_sarg_type1.drop(['Genotype'], axis=1, inplace=False)
    df_sarg_type2=df_sarg_type2[["Contig_ID","query","class","Phenotype"]]
    
    #processing ARG ranking index
    #open arg rank index file
    df_argrank=pd.read_csv(input_argrank, sep="\t", header=0)
    
    #chose columns
    df_arg1=df_argrank[["ARG","Rank"]]
    
    #add rank
    df_argrank2=pd.merge(df_sarg_type2,df_arg1,left_on="query",right_on="ARG",how="inner")
    df_argrank2=df_argrank2.drop(["ARG"],axis=1,inplace=False)
    df_argrank2.columns=["Contig_ID","query","class","Phenotype","SARG_Rank"]
    
    #save
    df_argrank2.to_csv(output + "/" + filebase +"_SARGrank_Protein60_Result.tsv", sep="\t")
    

if __name__=="__main__":   
    import pandas as pd
    import os
    import path    
    
    #globle settings
    #default settings
    ##flexible settings
    input_dir="/lomi_home/gaoyang/software/CompRanking/test" #-i manual setting
    project_prefix="CompRanking" #-p manual setting
    ##fixed settings
    input_argrank="/lomi_home/gaoyang/software/CompRanking/databases/SARG/ARG_rank.txt"
    input_sarg_length="../databases/SARG/SARG.db.fasta.length"
    input_sarg_structure="../databases/SARG/SARG.structure.txt" #"/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.structure.txt"
    output=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking") #fixed path
    file_abs_path=path.file_abs_path_list_generation(input_dir) #fixed path
    file_name_base = path.file_base_acquire(file_abs_path) #fixed path
    
    #arg rank processing
    # input_sarg="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/ARGranking/ERR1191817.contigs_5M_contigs_SARG_Protein_diamond.txt"
    for i in file_name_base:
        input_sarg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking", i+"_5M_contigs_SARG_Protein_diamond.txt")
        arg_rank(input_sarg, input_sarg_length,input_sarg_structure, input_argrank,i, output)
    
    
    



