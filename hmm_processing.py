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

#acquire abs path of sample
def file_abs_path_list_generation(input_dir):
    input=input_dir+"/*fa"
    file_abs_path_list=glob.glob(input)
    return file_abs_path_list

#acquire sample name base
def file_base_acquire(file_abs_path): #get names of all the sample files
    base_list=[]
    #Acquire file_base
    for i in file_abs_path:
        fileBase_1=i.split("/")[-1]
        fileBase_2=fileBase_1.strip(".fa")
        base_list.append(fileBase_2)
    return base_list

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



if __name__ == '__main__':
    import pandas as pd
    import re
    import glob
    import os
    
    #获取输入样本的文件名base
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    file_abs_path=file_abs_path_list_generation(input_dir)
    file_name_base = file_base_acquire(file_abs_path)
    prefix="CompRanking"
    
    # input="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/Virulence/111.hmm.txt"
    output=input_dir +"/Virulence/"
    suffix="_5M_contigs"
    for i in file_name_base:
        hmm_file=os.path.join(input_dir,prefix,"Virulence",i+suffix+"_tmp_hmm.txt")
        output_file=os.path.join(input_dir,prefix,"Virulence","hmm_result",i+suffix+"_hmm.csv")
        change_tab_hmmscan(input_hmmscan=hmm_file,output_hmm_csv=output_file)
    
    
    





