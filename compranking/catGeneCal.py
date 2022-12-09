#!/usr/bin/env python
# title             :catGeneCal.py 
# description       :calculate relative abundance of genes
# author            :Gaoyang Luo
# date              :202201101
# version           :1.0
# usage             :python catGeneCal.py  -i <input_dir> 
#                                          -f <file_dir> 
#                                          -p <project_prefix> 
#                             
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
import path
import optparse

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")
parser.add_option("-f", "--fileDir", action = "store", type = "string", dest = "subfile_dir",
				 help = "directory of tmp calculation results")


(options, args) = parser.parse_args()
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
subfile_dir=options.file_dir
output=subfile_dir
file_abs_path=path.file_abs_path_list_generation(input_dir)
file_name_base = path.file_base_acquire(file_abs_path)

#########concat all the abu result############
#concat ARG result
#concat 16S
name_list_16S=[]
for i in file_name_base:
    name_list_16S.append(i+"_ARG_16sAbu_tmp.txt")
init=0
df_main=pd.read_csv(os.path.join(subfile_dir,name_list_16S[0]),sep="\t", header=None)
df_main.columns=["type",name_list_16S[0]]
for i,name in enumerate(name_list_16S):
    if i < len(name_list_16S)-1:
        init+=1
        if name_list_16S[init]:
            df_2=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_16S[init]),sep="\t", header=None)
            df_2.columns=["type",name_list_16S[init]]
            df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
#save 16s subtype abu
df_main.to_csv(os.path.join(
            subfile_dir,
                    project_prefix+"_Abundance_ARGs_subtypes_16S.txt"),sep="\t",index=None)
#cal rpkm
name_list_rpkm=[]
for i in file_name_base:
    name_list_rpkm.append(i+"_ARG_rpkmAbu_tmp.txt")
init=0
df_main=pd.read_csv(os.path.join(subfile_dir,name_list_rpkm[0]),
                    sep="\t", header=None)
df_main.columns=["type",name_list_rpkm[0]]
for i,name in enumerate(name_list_rpkm):
    if i < len(name_list_16S)-1:
        init+=1
        if name_list_rpkm[init]:
            df_2=pd.read_csv(os.path.join(subfile_dir,name_list_rpkm[init]),
                                sep="\t", header=None)
            df_2.columns=["type",name_list_rpkm[init]]
            df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
#save rpkm subtype abu
df_main.to_csv(os.path.join(
                subfile_dir,
                        project_prefix+"_Abundance_ARGs_subtypes_rpkm.txt"),sep="\t",index=None)
#concat MGE result
#concat 16S
name_list_16S=[]
for i in file_name_base:
    name_list_16S.append(i+"_MGE_16sAbu_tmp.txt")
init=0
df_main=pd.read_csv(os.path.join(subfile_dir,name_list_16S[0]),sep="\t", header=None)
df_main.columns=["type",name_list_16S[0]]
for i,name in enumerate(name_list_16S):
    if i < len(name_list_16S)-1:
        init+=1
        if name_list_16S[init]:
            df_2=pd.read_csv(os.path.join(subfile_dir,name_list_16S[init]),sep="\t", header=None)
            df_2.columns=["type",name_list_16S[init]]
            df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
#save 16s subtype abu
df_main.to_csv(os.path.join(
            subfile_dir,
                    project_prefix+"_Abundance_MGEs_subtypes_16S.txt"),sep="\t",index=None)
#cal rpkm
name_list_rpkm=[]
for i in file_name_base:
    name_list_rpkm.append(i+"_MGE_rpkmAbu_tmp.txt")
init=0
df_main=pd.read_csv(os.path.join(subfile_dir,name_list_rpkm[0]),
                    sep="\t", header=None)
df_main.columns=["type",name_list_rpkm[0]]
for i,name in enumerate(name_list_rpkm):
    if i < len(name_list_16S)-1:
        init+=1
        if name_list_rpkm[init]:
            df_2=pd.read_csv(os.path.join(subfile_dir,name_list_rpkm[init]),
                                sep="\t", header=None)
            df_2.columns=["type",name_list_rpkm[init]]
            df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
#save rpkm subtype abu
df_main.to_csv(os.path.join(
                subfile_dir,
                        project_prefix+"_Abundance_MGEs_subtypes_rpkm.txt"),sep="\t",index=None)
