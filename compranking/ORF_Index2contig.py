#!/usr/bin/env python
# title             :ORF_Index2contig.py -> index.ipynb
# description       :To link orf index to contig index
# author            :Gaoyang Luo
# date              :20231202
# version           :1.0
# usage             :python ORF_Index2contig.py -i <input_dir> -p <project_prefix>
# required packages :pandas, numpy, multiprocessing
# notification: enjoy yourself
#=============================================================================

import pandas as pd
import numpy as np
import os, sys
import optparse
import multiprocessing
import path

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")

(options, args) = parser.parse_args()
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
    
def orf2contig(sample_list):
    sample_name=sample_list
    file_path=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs/"+ sample_name +"_5M_contigs.fna2faa.index")
    df_contig=pd.read_csv(file_path, 
                    sep='\t', header=0)
    df_contig.columns=["contig", "orf"]
    array=np.array(df_contig)
    array=array.tolist()
    contig_orf_dic={}
    for i in array:
        tmp=i[0].split("_")[0:-1]
        contig_orf_dic.setdefault(str(i[1]),str(tmp[0]+"_"+tmp[1]))
    
    df_contig_orf = pd.DataFrame(pd.Series(contig_orf_dic))
    df_contig_orf.reset_index(inplace=True)

    df_contig_orf.columns=["orf", "contigs"]
    df_contig_orf=df_contig_orf[["contigs","orf"]]
    
    df_contig_orf.to_csv(os.path.join(
                        input_dir,project_prefix,
                            "CompRanking_intermediate/preprocessing/5M_contigs",
                                sample_name +"_5M_contigs.index"),
                                    sep="\t",header=False, index=None)


def multi_info_sum():
    openthreads = len(sample_list) 
    exfiles = []
    for i in range(openthreads):
        worker = multiprocessing.Process(target=orf2contig,args=([sample_list[i]]))
        worker.start()
        print("Now processing:{}".format(sample_list[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()
    
if __name__ == "__main__":   
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    sample_list= path.file_base_acquire(file_abs_path) #sample name without suffix .fa
    
    multi_info_sum()
    
    # for i in *gff; do base=${i%%.gf*}; sed -i '/^#/d' ${i}; cut -f 1,9 ${i} |cut -d';' -f1| sed 's/ID=//' > ${base}.fna2faa.index; done
    # for i in *index; do base=${i/.index/}; mv ${i} ${base}.fna2faa.index; done