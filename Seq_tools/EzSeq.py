#!/usr/bin/env python
# title             :Seq_extract.py -> seq_extract.ipynb
# description       :extract seq of interest with seq id input and fasta files
# author            :Gaoyang Luo
# date              :20221113
# version           :1.0
# usage             :import EzSeq
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
import re
from compranking import path
from Bio import SeqIO
import optparse

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director containing raw fasta files, the same as the input of cpr input directory")
parser.add_option("-c", "--class_type", action = "store", type = "string", dest = "class_type",
				 help = "class:ARG or PATH")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")

(options, args) = parser.parse_args()

####################################preseting#################################### 
#path configeration
input_dir=options.input_dir
class_type=options.class_type
project_prefix=options.project_prefix
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name

if (options.class_type is None):
    raise TypeError("there must a type chosen...") #choose something you want to extract

def seq_id_extract_ARG_host(sumTable):
    """
    extract sequences of interest from fasta files
    First, you need to make sure what kind of seqs you want to extract
    according to the summary table of CompRanking output.
    This function will give you an output type of dataframe 
    and a fasta file containing sequences of interst
    """
    #load input
    seq_id_ARG_host=[]
    df_summary=pd.read_csv(sumTable, sep="\t",header=0)
    ARG_filter=df_summary[df_summary.ARG_prediction	!= "-"]
    for i ,name in ARG_filter.iterrows():
        seq_id_ARG_host.append(name["Contig"])
    return seq_id_ARG_host
    
    
def seq_id_extract_PATH_host(sumTable):
    """
    extract sequences of interest from fasta files
    First, you need to make sure what kind of seqs you want to extract
    according to the summary table of CompRanking output.
    This function will give you an output type of dataframe 
    and a fasta file containing sequences of interst
    """
    seq_id_PATH_host=[]
    df_summary=pd.read_csv(sumTable, sep="\t",header=0)
    PATH_filter=df_summary[df_summary.Pathogenicity =="Pathogenic"]
    for i ,name in PATH_filter.iterrows():
        seq_id_PATH_host.append(name["Contig"])
    return seq_id_PATH_host

def subSeqFileGeneration(FASTA_file):
	dicfq = {}
	for i in FASTA_file:
		if re.match(">",i):
			dicfq[i]=""
			flag = i
		else:
			dicfq[flag] = dicfq[flag]+i
	return dicfq
    
    
if __name__ == "__main__":
    #read file names
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    
    if class_type == "ARG":
        try:
            for name in file_name_base:
                inputFASTA=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",name+"_5M_contigs.fa")
                inputTable=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs","CompRanking_"+name+"_Summary.tsv")
                print("You must want to extract ARG host sequences...")
                ARG_host_list=seq_id_extract_ARG_host(inputTable)
                #generate fasta file
                # #写入“ContigsID.txt”
                # filename=file_name_base + "_ARGHost_ContigsID.txt"
                # with open (filename, "w") as f1:
                #     for i in ARG_host_list:
                #         f1.write(i + "\n")
                
                #write into new fasta file
                fasta_file = open(inputFASTA,"r")
                f_list = ARG_host_list
                dickf = subSeqFileGeneration(fasta_file)
                for i in f_list:
                    i = i.strip()
                    for j in dickf.keys():
                        if re.match(">"+i,j):
                            print(j,end="")
                            print(dickf[j])
                fasta_file.close()            
        except:
            raise TypeError("Something wrong...")
        
    if class_type == "PATH":
        try:
            for name in file_name_base:
                inputFASTA=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",name+"_5M_contigs.fa")
                inputTable=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs","CompRanking_"+name+"_Summary.tsv")
                print("You must want to extract ARG host sequences...")
                PATH_host_list=seq_id_extract_PATH_host(inputTable)
                #generate fasta file
                # #写入“ContigsID.txt”
                # filename=file_name_base + "_ARGHost_ContigsID.txt"
                # with open (filename, "w") as f1:
                #     for i in ARG_host_list:
                #         f1.write(i + "\n")
                
                #write into new fasta file
                fasta_file = open(inputFASTA,"r")
                f_list = PATH_host_list
                dickf = subSeqFileGeneration(fasta_file)
                for i in f_list:
                    i = i.strip()
                    for j in dickf.keys():
                        if re.match(">"+i,j):
                            print(j,end="")
                            print(dickf[j])
                fasta_file.close()            
        except:
            raise TypeError("Something wrong...")
            
            
#python EzSeq.py -i /lomi_home/gaoyang/software/CompRanking/test
    

    

    
        