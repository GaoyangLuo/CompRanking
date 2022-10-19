#!/usr/bin/env python
# title             :Virulence_processing.py
# description       :process all the outputs
# author            :Gaoyang Luo
# date              :202201019
# version           :1.0
# usage             :from Virulence_processing import *
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import path

def VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,output):
    #load index
    df_contig=pd.read_csv(input_contig,sep="\t",header=None)
    #load VFDB output
    df_ERR_VFDB_output=pd.read_csv(input_ERR_VFDB_output,sep="\t",header=None)
    df_ERR_VFDB_output.columns=['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
    df_ERR_VFDB_output60= df_ERR_VFDB_output[df_ERR_VFDB_output.identity > 60]
    df_ERR_VFDB_output_algn=df_ERR_VFDB_output60[df_ERR_VFDB_output60.alignLen>40]
    
    #merge with Contig ID
    df_VFDB_merge=pd.merge(df_contig,df_ERR_VFDB_output_algn,left_on=1,right_on="id",how="inner")
    df_VFDB_merge=df_VFDB_merge[[0,"id","sub_id"]] #filter 3
    
    #load CompRankingVirluenceSummary
    df_cpr_VF_sum=pd.read_csv(input_cpr_VF_sum, index_col=0,sep="\t",header=0)
    
    #merge CompRankingVirluenceSummary and VFDB output
    df_VFDB_output_annote=pd.merge(df_VFDB_merge,df_cpr_VF_sum,left_on="sub_id",right_on="fasID",how="left")
    df_VFDB_output_annote=df_VFDB_output_annote[[0,"id","sub_id","VFID","Virulence_factor","VF_Name","Bacteria","VFcategory","ESKAPE"]]
    df_VFDB_output_annote.columns=["Contigs","id","sub_id","VFID","Virulence_factor","VF_Name","Bacteria","VFcategory","ESKAPE"]
    
    #save
    
    VFfile=df_VFDB_output_annote.to_csv(output + "/CompRanking_Virulence_VFDB_output.csv",sep="\t",index=None)
    
    return VFfile
    
        
    

if __name__ == "__main__":
    input_contig="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.index"
    input_ERR_VFDB_output="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/Virulence/ERR1191817.contigs_5M_contigs_VFDB_setA1e-5.out"
    input_cpr_VF_sum="/lomi_home/gaoyang/software/CompRanking/databases/CompRanking_VirulenceDB/CompRanking_HMM_Virulence_Summary.csv"
    output="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_result"
    
    VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,output)
    
    