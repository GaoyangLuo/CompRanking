#!/usr/bin/env python
# title             :Seq_extract.py -> seq_extract.ipynb
# description       :extract seq of interest with seq id input and fasta files
# author            :Gaoyang Luo
# date              :20221113
# version           :1.0
# usage             :EzSeq -i <input_director> - <what type of seq u wannna extract>
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
import re
import optparse
import sys
import datetime
from Bio import SeqIO
sys.path.append("..")
from compranking import path

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

def writeSeq2Dic(FASTA_file):
    seq_store={}
    for record in SeqIO.parse(FASTA_file, 'fasta'):
        seq_store.setdefault(">" + str(record.id),
                                    str(record.seq))
    return seq_store

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
            start = datetime.datetime.now() #time start
            
            for name in file_name_base:
                #read summary table
                inputFASTA=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",name+"_5M_contigs.fa")
                inputTable=os.path.join(input_dir,project_prefix,"CompRanking_result","CompRanking_"+name+"_Summary.tsv")
                print("You must want to extract PATH host sequences...")
                ARG_host_list=seq_id_extract_ARG_host(inputTable)
                
                #write into new fasta file
                fasta_file = open(inputFASTA,"r")
                f_list = ARG_host_list
                f_list=set(f_list)
                dicfq = subSeqFileGeneration(fasta_file)
                seq_store = writeSeq2Dic(inputFASTA)
                with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                                       name+"_ARG_tmp_host.txt"),"w") as file:
                    for i in f_list:
                        i = i.strip()
                        i = i.split(" ")[0]
                        # for j in dicfq.keys():
                        #     if re.match(">"+i,j):
                        #         print(j,end="")
                        #         print(dicfq[j])
                        #         file.write(">"+i + '\n'  + str(
                        #     dicfq[j]))
                        if (">"+i) in seq_store.keys():
                            print(">"+ i )
                            print(seq_store[">"+i])
                            file.write(">"+ i + '\n'  + str(
                                seq_store[">"+i]) + '\n')
                        
                fasta_file.close()
                file.close()         
                
            end = datetime.datetime.now() #time end
            print("ARG extract cost time: {}".format(end-start))   
        except:
            raise TypeError("Something wrong...")
        
    if class_type == "PATH":
        try:
            start = datetime.datetime.now() #time start
            for name in file_name_base:
                #read summary table
                inputFASTA=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",name+"_5M_contigs.fa")
                inputTable=os.path.join(input_dir,project_prefix,"CompRanking_result","CompRanking_"+name+"_Summary.tsv")
                print("You must want to extract Pathogen host sequences...")
                PATH_host_list=seq_id_extract_PATH_host(inputTable)
                                
                #write into new fasta file
                fasta_file = open(inputFASTA,"r")
                f_list = PATH_host_list
                f_list = set(f_list)
                dicfq = subSeqFileGeneration(fasta_file)
                seq_store = writeSeq2Dic(inputFASTA)
                with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                                       name+"_PATH_tmp_host.txt"),"w") as file:
                    for i in f_list:
                        i = i.strip()
                        i = i.split(" ")[0]
                        if (">"+i) in seq_store.keys():
                            print(">"+ i )
                            print(seq_store[">"+i])
                            file.write(">"+ i + '\n'  + str(
                                seq_store[">"+i]) + '\n')
                        
                fasta_file.close()
                file.close()         
                
            end = datetime.datetime.now() #time end
            print("PATH extract cost time: {}".format(end-start))       
        except:
            raise TypeError("Something wrong...")
            
            
#python EzSeq.py -i /lomi_home/gaoyang/software/CompRanking/test
    

    

    
        