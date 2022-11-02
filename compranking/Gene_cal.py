#!/usr/bin/env python
# title             :Gene_cal.py -> Gene_cal.ipynb
# description       :calculate relative abundance of genes
# author            :Gaoyang Luo
# date              :202201101
# version           :1.0
# usage             :import Gene_cal
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import path
from Bio import SeqIO


def get_DB_DeepARG_len(input_deeparg_length):
    #load_Deeparg_structure
    df_deepARG_strucure=pd.read_csv(input_deeparg_length,sep="\t",header=None)
    df_deepARG_strucure
    df_deepARG_strucure["id"]=df_deepARG_strucure[0].str.split("|",expand=True)[0]
    #write dic
    DB_deepARG_length={}
    for i, name in df_deepARG_strucure.iterrows():
        DB_deepARG_length.setdefault(str(name["id"]), name[1])
    return DB_deepARG_length
    
    
def get_DB_SARG_len(input_sarg_structure):
    #load SARG_structure
    df_sarg_structure=pd.read_csv(input_sarg_structure, sep="\t", header=None)
    df_sarg_structure=df_sarg_structure[[1,2]]
    #write dic
    DB_SARG_length={}
    for i, name in df_sarg_structure.iterrows():
        DB_SARG_length.setdefault(str(name[1]), name[2])
    return DB_SARG_length     

def get_MobilOG_len(input_mobileOG_structure):
    #load mobileOG_structure
    df_mobileOG_structure=pd.read_csv(input_mobileOG_structure,sep="\t", header=0)
    #write dic
    DB_MobileOG_length={}
    for i, name in df_mobileOG_structure.iterrows():
        DB_MobileOG_length.setdefault(str(name["mobileOG_ID"]), name["length"])
    return DB_MobileOG_length
    
    
def RB_gene_sum(DB_deepARG_length,DB_SARG_length, DB_MobileOG_length, 
                input_AMR_sum,input_kk2,input_deeparg_sure,input_rgi,input_SARG,filebase):
    #load final output
    df_AMR_sum=pd.read_csv(input_AMR_sum,sep="\t",header=0)
    df_AMR_hit=df_AMR_sum[df_AMR_sum.ARG_prediction != "-"]
    df_AMR_hit1=df_AMR_hit[["ORF_ID","Database"]]
    df_AMR_hit1["db_final"]=df_AMR_hit1["Database"].str.split("/", expand=True)[0]
    
    #record hit database and orf_id
    Record_db_orf={}
    for i, name in df_AMR_hit1.iterrows():
        Record_db_orf.setdefault(str(name["ORF_ID"]), str(name["db_final"]))
    
    #calculate 16s copies
    #load kk2
    kraken=input_kk2
    copy_16S = 1
    gene_length = 1550 # 16S

    # metagenomes
    for lines in open(kraken,'r'):
        content = lines.split('\n')[0].split('\t')
        if 'Bacteria' in lines:
            copy_16S = float(content[1]) # pair end
            break
        
    #load deeparg
    df_deeparg_sure=pd.read_csv(input_deeparg_sure, sep="\t")
    #load rgi
    df_RGI=pd.read_csv(input_rgi, sep="\t")
    df_RGI=df_RGI.fillna("-")
    #load sarg
    df_SARG=pd.read_csv(input_SARG, sep="\t", header=0,index_col=0)
    df_SARG["class"]=df_SARG["class"].str.upper()
    df_SARG.columns=["ORF","query","class","Phenotype","ARG_rank"]

    #get each ARG output dic
    deepARG_output_dic={}
    SARG_output_dic={}
    for i, name in df_deeparg_sure.iterrows():
        deepARG_output_dic.setdefault(str(name["read_id"]), str(name["best-hit"].split("|")[0]))
    for i, name in df_SARG.iterrows():
        SARG_output_dic.setdefault(str(name["ORF"]), str(name["query"]))

    #generate DB_xxx_length_res
    DB_deepARG_length_res={}
    DB_SARG_length_res={}
    DB_CARD_length_res={}
    #deeparg_output_ARG_cal
    for i in deepARG_output_dic:
        DB_deepARG_length_res.setdefault(i,DB_deepARG_length[deepARG_output_dic[i]])
    #SARG_output_ARG_cal
    for i in SARG_output_dic:
        DB_SARG_length_res.setdefault(i,DB_SARG_length[SARG_output_dic[i]])
    #RGI_output_ARG_cal
    for i, name in df_RGI.iterrows():
        DB_CARD_length_res.setdefault(str(name["ORF_ID"]), len(str(name["CARD_Protein_Sequence"])))
    
    #cal ARGs relative abundance    
    abundance_arg=0
    for orf in Record_db_orf:
        find_db=''
        if Record_db_orf[orf]:
            find_db=Record_db_orf[orf]
            if find_db=="DeepARG":
                abundance_arg += (gene_length/copy_16S)*(1/DB_deepARG_length_res[orf])
            elif find_db=="RGI":
                abundance_arg += (gene_length/copy_16S)*(1/DB_CARD_length_res[orf])
            elif find_db=="SARG":
                abundance_arg += (gene_length/copy_16S)*(1/DB_SARG_length_res[orf])
            else:
                continue
    print(abundance_arg)
    
    abundance_MGE=0
    
    
    result=[abundance_arg,abundance_MGE]
    
    return result
    
    
if __name__ == "__main__":
    #global settings
    #gloab settings
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    output=os.path.join(input_dir,"CompRanking/CompRanking_result")
    project_prefix="CompRanking"
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    
    #input_db
    input_deeparg_length="/lomi_home/gaoyang/software/CompRanking/databases/deepargdata1.0.2/database/v2/features.gene.length"
    input_sarg_structure="/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.db.fasta.length"
    input_mobileOG_structure="/lomi_home/gaoyang/software/CompRanking/databases/MobileOG-db/MobileOG-db_structure.tsv"
    
    #input_result
    input_AMR_sum="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_result/CompRanking_ERR1191817.contigs_AMR_prediction.tsv"
    input_kk2="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_report_kk2_mpaStyle.txt"
    input_deeparg_sure="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/DeepARG/ERR1191817.contigs_5M_contigs_DeepARG.out.mapping.ARG"
    input_rgi="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/RGI/ERR1191817.contigs_5M_contigs.RGI.out.txt"
    input_SARG="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/ARGranking/ERR1191817.contigs_SARGrank_Protein60_Result.tsv"
    
    #calculate relative abundance of functional genes
    for i in file_name_base:
        #load ARGs result
        input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
        input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein60_Result.tsv")
        input_deeparg_sure=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG", i+"_5M_contigs_DeepARG.out.mapping.ARG")
        input_kk2=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs", i+"_report_kk2_mpaStyle.txt")
        input_AMR_sum=os.path.join(input_dir,project_prefix,"CompRanking_result","CompRanking_"+i+"_AMR_prediction.tsv")
        
        #load reference length
        DB_deepARG_length = get_DB_DeepARG_len(input_deeparg_length)
        DB_SARG_length =get_DB_SARG_len(input_sarg_structure)
        DB_MobileOG_length=get_MobilOG_len(input_mobileOG_structure)
        
        
        result=RB_gene_sum(DB_deepARG_length,DB_SARG_length, DB_MobileOG_length, 
                input_AMR_sum,input_kk2,input_deeparg_sure,input_rgi,input_SARG,i)
        output="\t".join(map(str, result))
    
        with open(os.path.join(input_dir,"CompRanking/CompRanking_result")+"Gene_Abundance_Sum.txt", "a") as f:
            f.write("\n" + i + "\t" + output)
        
        
    