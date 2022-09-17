#!/usr/bin/env python
# title             :hmm_processing.py
# description       :process hmm output for virulence prediction
# author            :Gaoyang Luo
# date              :202209012
# version           :1.0
# usage             :import hmm_processing
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import cpr_run as cpr

#acquire abs path of sample
# def file_abs_path_list_generation(input_dir):
#     input=input_dir+"/*fa"
#     file_abs_path_list=glob.glob(input)
#     return file_abs_path_list

# #acquire sample name base
# def file_base_acquire(file_abs_path): #get names of all the sample files
#     base_list=[]
#     #Acquire file_base
#     for i in file_abs_path:
#         fileBase_1=i.split("/")[-1]
#         fileBase_2=fileBase_1.strip(".fa")
#         base_list.append(fileBase_2)
#     return base_list

def change_tab_hmmscan(input_hmmscan, output_hmm_csv):
    #input_hmmscan is the output file of hmm
    #读入hmmscan原始文件
    #sed '/^#/ d' ERR1191817.contigs_5M_contigs.hmmscan | sed 's/ \+/\\t/g' > 111.hmm.txt
    df=pd.read_csv(input_hmmscan, 
                    sep=r'\\t', header=None, engine='python')
    df_tmp=df.iloc[:,[0,2,4,5]]
    df_tmp.columns=["ID","query","evalue","score"] #rename column
    empty=[]
    for i in df_tmp["query"]:
        i=re.sub("_.*","",i)
        empty.append(i)
    df_tmp_copy=df_tmp.copy() #copy to avoid error warning
    df_tmp_copy["query"]=empty
    #save result
    output_hmm=df_tmp_copy.to_csv(output_hmm_csv, sep="\t")
    
    return output_hmm

def VF_predition(input_hmmcsv, positive_domains,name_VF_output):
    df_hmm=pd.read_csv(input_hmmcsv,  #"./test/CompRanking/Virulence/hmm_result/ERR1191817.contigs_5M_contigs_hmm.csv"
                  sep='\t', header=0, engine="python")
    positive_query=pd.read_csv(positive_domains,  #"./databases/models_and_domains/positive_domains.tsv"
                  sep='\t', header=None, engine="python")
    positive_query_1=positive_query.copy()
    positive_query_1["class"]="pathogenic"
    positive_query_1.columns=["query_positive","class"]
    positive_merge_tmp=pd.merge(df_hmm,positive_query_1,left_on="query", right_on="query_positive",how="inner")
    Virulence_non_redundant=positive_merge_tmp.loc[:,["ID","class"]]
    Virulence_non_redundant_unique=Virulence_non_redundant.drop_duplicates()
    #"./test/CompRanking/Virulence/ERR1191817.contigs_5M_contigs_VF_Prediction.csv"
    output_VF_prediction=Virulence_non_redundant_unique.to_csv(name_VF_output, sep="\t")
    
    return output_VF_prediction

if __name__ == '__main__':
    #获取输入样本的文件名base
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    file_abs_path=cpr.file_abs_path_list_generation(input_dir)
    file_name_base = cpr.file_base_acquire(file_abs_path)
    prefix="CompRanking"
    
    #generate hmm.csv
    # input="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/Virulence/111.hmm.txt"
    suffix="_5M_contigs"
    for i in file_name_base:
        hmm_file=os.path.join(input_dir,prefix,"Virulence",i+suffix+"_tmp_hmm.txt")
        output_file=os.path.join(input_dir,prefix,"Virulence","hmm_result",i+suffix+"_hmm.csv")
        change_tab_hmmscan(input_hmmscan=hmm_file,output_hmm_csv=output_file)
    
    #hmmcsv2predcsv
    positive_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/models_and_domains/positive_domains.tsv") 
    for i in file_name_base:
        input=os.path.join(input_dir,prefix,"Virulence","hmm_result",i+suffix+"_hmm.csv")
        output=os.path.join(input_dir,prefix,"Virulence","hmm_result",i+suffix+"_VF_Prediction.csv")
        VF_predition(input_hmmcsv=input,positive_domains=positive_path,name_VF_output=output)
        
        
    
    





