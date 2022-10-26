#!/usr/bin/env python
# title             :Virulence_processing.py -> VFDB&PATH_processing.ipynb
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

def VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,output,filebase):
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
    
    #PATH-like sequences processing
    #load patric output
    df_patric=pd.read_csv(input_patric, sep="\t", header=None)
    df_patric.columns=['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
    #iden60 and alen150
    df_patric_iden60=df_patric[df_patric.identity > 80]
    df_patric_den60_alen150=df_patric_iden60[df_patric_iden60.alignLen > 150]
    df_PATH_contig=df_patric_den60_alen150[["id"]].drop_duplicates()
    df_PATH_contig["Pathogenicity"]="Pathogenic"
    #first merge VFs with index
    df_VFs_contig=pd.merge(df_contig,df_VFDB_output_annote,left_on=1,right_on="id",how="left")
    df_VFs_PATH_contig=pd.merge(df_VFs_contig,df_PATH_contig,left_on=0,right_on="id",how="left")
    #drop
    df_VFs_PATH_contig=df_VFs_PATH_contig.drop(["Contigs","id_x","sub_id","VFID","id_y"],axis=1,inplace=False)
    #fill pathogenicity
    df_VFs_PATH_contig=df_VFs_PATH_contig.fillna("-")
    for index ,name in df_VFs_PATH_contig.iterrows():
        if df_VFs_PATH_contig["Virulence_factor"][index] != "-":
            df_VFs_PATH_contig["Pathogenicity"][index]="Pathogenic"
    #rename column
    df_VFs_PATH_contig.columns=["Contig","ORF","Virulence_factor","VF_Name","VF_Host_Bacteria","VFcategory","ESKAPE","Pathogenicity"]
    
    #save
    VFfile=df_VFDB_output_annote.to_csv(output + "/CompRanking_"+filebase+"_Virulence_VFDB_output.tsv",sep="\t",index=None)
    
    return VFfile
    
        
    

if __name__ == "__main__":
    import pandas as pd
    import re
    import glob
    import os
    import path
    
    input_contig="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.index"
    input_ERR_VFDB_output="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/Virulence/ERR1191817.contigs_5M_contigs_VFDB_setA1e-5.out"
    
    input_patric="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/Virulence/PATRIC/ERR1191817.contigs_5M_contigs_PATRIC.out"
    
    #gloab settings
    #flexible setting
    input_dir="/lomi_home/gaoyang/software/CompRanking/test" 
    #fixed setting
    input_cpr_VF_sum="../databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" 
    output=os.path.join(input_dir,"CompRanking/CompRanking_result")
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    
    
    # print(file_abs_path)
    for i in file_name_base:
        VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,output,i)
    
    