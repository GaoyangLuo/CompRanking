#!/usr/bin/env python
# title             :baseInfoExtra_Celleq.py -> baseInfoExtra.ipynb
# description       :This script is for benchmarking using average genome size as normalization bases to gener
#                    generate risk score
# author            :Gaoyang Luo
# date              :202301210
# version           :1.0
# usage             :python baseInfoExtra_Celleq.py -i <input_dir> -p <project_prefix>
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
import path
import math
import optparse
import multiprocessing

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_dir",
				 help = "director tontained output files")     


(options, args) = parser.parse_args()
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
output=options.output_dir
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.output_dir is None):
    output = os.path.join(input_dir,project_prefix,"CompRanking_result") #default output directory

    
# #假设输入文件为示例文件，放在for循环的开头第一层 for i = samle_name
# file_abs_path=path.file_abs_path_list_generation(input_dir)
# sample_list= path.file_base_acquire(file_abs_path) #sample name without suffix .fa
# write_file=output + "/CompRanking_"+ project_prefix + "_Contigs_Risk_Summary.txt"
# with open(write_file, "w") as f1:
#     f1.write("sample_name/index\tnContigs\tnARGs_contigs\tnMGEs_contig\tnMGEs_plasmid_contig\tnMGEs_phage_contigs\tnPAT_contigs\tnARGs_MGEs_contig\tnARGs_MGEs_plasmid_contigs\tnARGs_MGEs_phage_contigs\tnARGs_MGEs_PAT_contigs\tfARG\tfMGE\tfMGE_plasmid\tfMGE_phage\tfPAT\tfARG_MGE\tfARG_MGE_plasmid\tfARG_MGE_phage\tfARG_MGE_PAT\tscore_pathogenic\tscore_phage\tscore_plasmid")

def info_sum(sample_list, genome_equivalents, copy_cell):
# for i in sample_list: #获取sample_name
    sample_name=sample_list
    print(sample_name)
    # 获取结果表
    file_path=os.path.join(input_dir,project_prefix,"CompRanking_result/CompRanking_"+ sample_name +"_Summary.tsv")
    df=pd.read_csv(file_path, 
                    sep='\t', header=0)
    #ncontigs_counts
    nContigs_num=len(df.Contig.unique()) # numbers of contigs
    #get genome equivalent
    nContigs=float(copy_cell)* float(genome_equivalents)* 0.01
    nContigs=float(copy_cell)
    print("Cell equivalent is {}".format(nContigs))
    #nARGs_contigs_counts
    #RankIII_Risk
    #ARGs_List&&ARGs_typeCount
    ARGs_List=[]
    arg_filter=df[df.ARG_prediction!="-"]
    nARGs_contigs=len(arg_filter.Contig.unique()) # numbers of contigs containing predicted ARGs
    #nMGEs_counts
    MGEs_filter=df[df.MGE_prediction!="unclassified" ]
    nMGEs_contigs=len(MGEs_filter.Contig.unique()) # numbers of contigs containing predicted MGEs
    #nMGEs_plasmids_counts
    MGEs_plasmid_filter=df[df.MGE_prediction=="plasmid" ]
    nMGEs_plasmid_contigs=len(MGEs_plasmid_filter.Contig.unique())
    #nMGEs_phage_counts
    MGEs_phage_filter=df[df.MGE_prediction=="phage" ]
    nMGEs_phage_contigs=len(MGEs_phage_filter.Contig.unique())
    #nPAT(pathogenic)_counts
    PAT_filter=df[df.Pathogenicity=="-"]
    nPAT_pathogenic_contigs=len(PAT_filter.Contig.unique())
    #nARGs_MGEs_counts
    #RankII_Risk
    # nARGs_MGEs_filter=arg_filter[arg_filter.Pathogenicity=="-"]
    nARGs_MGEs_filter=arg_filter[arg_filter.MGE_prediction!="unclassified"]
    nARGs_MGEs_contigs=len(nARGs_MGEs_filter.Contig.unique())
    #nARGs_MGEs_plasmid_counts
    nARGs_MGEs_plasmid_filter=arg_filter[arg_filter.MGE_prediction=="plasmid"]
    nARGs_MGEs_plasmid_contigs=len(nARGs_MGEs_plasmid_filter.Contig.unique())
    #nARGs_MGEs_phage_counts
    nARGs_MGEs_phage_filter=arg_filter[arg_filter.MGE_prediction=="phage"]
    nARGs_MGEs_phage_contigs=len(nARGs_MGEs_phage_filter.Contig.unique())
    #nARGs_MGEs_PAT(pathogenic)
    #RankI_Risk
    nARGs_MGEs_PAT_filter=arg_filter[arg_filter.MGE_prediction!="unclassified"]
    nARGs_MGEs_PAT_filter=nARGs_MGEs_PAT_filter[nARGs_MGEs_PAT_filter.Pathogenicity=="Pathogenic"]
    # nARGs_MGEs_PAT_filter=nARGs_MGEs_filter[nARGs_MGEs_filter.Pathogenicity=="Pathogenic"]
    nARGs_MGEs_PAT_contigs=len(nARGs_MGEs_PAT_filter.Contig.unique())
    #risk calculation
    # nContigs
    # nARGs_contigs
    # nMGEs_contigs
    # nPAT_contgis
    # nARGs_MGEs_contigs
    # nARGs_MGEs_PAT_contigs
    fMGE=float(nMGEs_contigs)/nContigs
    fMGE_plasmid=float(nMGEs_plasmid_contigs)/nContigs
    fMGE_phage=float(nMGEs_phage_contigs)/nContigs
    fPAT_pathogenic=float(nPAT_pathogenic_contigs)/nContigs
    fARG = float(nARGs_contigs)/nContigs
    fARG_MGE = float(nARGs_MGEs_contigs)/nContigs
    fARG_MGE_plasmid=float(nARGs_MGEs_plasmid_contigs)/nContigs
    fARG_MGE_phage=float(nARGs_MGEs_phage_contigs)/nContigs
    fARG_MGE_PAT= float(nARGs_MGEs_PAT_contigs)/nContigs
    #caclulate distance between sample of interest and theriotic point
    # distance_all = math.sqrt((0.01 - fARG)**2 + (0.01 - fARG_MGE)**2 + (0.01 - fARG_MGE_allPAT)**2)
    experience_para=float(0.015)
    distance_pathogenic = math.sqrt((0.01 - fARG)**2 + (experience_para - fARG_MGE)**2 + (0.01 - fARG_MGE_PAT)**2)
    distance_phage=math.sqrt((0.01 - fARG)**2 + (experience_para - fARG_MGE_phage)**2 + (0.01 - fARG_MGE_PAT)**2)
    distance_plasmid= math.sqrt((0.01 - fARG)**2 + (experience_para - fARG_MGE_plasmid)**2 + (0.01 - fARG_MGE_PAT)**2)
    #calculate risk score
    # score_all = 1.0 / ( (2 + math.log10(distance_all))**2 )
    score_pathogenic= 1.0 / ( (2 + math.log10(distance_pathogenic))**2 )
    score_phage= 1.0 / ( (2 + math.log10(distance_phage))**2 )
    score_plasmid= 1.0 / ( (2 + math.log10(distance_plasmid))**2 )

    result=[nContigs,nContigs_num, nARGs_contigs, nMGEs_contigs,nMGEs_plasmid_contigs,nMGEs_phage_contigs, nPAT_pathogenic_contigs, nARGs_MGEs_contigs,nARGs_MGEs_plasmid_contigs, nARGs_MGEs_phage_contigs,nARGs_MGEs_PAT_contigs, fARG, fMGE,fMGE_plasmid,fMGE_phage, fPAT_pathogenic, fARG_MGE, fARG_MGE_plasmid,fARG_MGE_phage,fARG_MGE_PAT,score_pathogenic,score_phage,score_plasmid]
    output="\t".join(map(str, result))

    with open(write_file, "a") as f:
        f.write("\n"+ sample_name + "\t" + output)

# def get_genome_equivalent(input_AGS, prefix_list):
#     genome_equivalents_dic={}
#     for index, j in enumerate(input_AGS): 
#         for lines in open(j,'r'):
#                             if lines.startswith('genome_equivalents'):
#                                 lines_set = lines.split('\n')[0].split('\t')
#                                 genome_equivalents = float(lines_set[1])
#                                 genome_equivalents_dic.setdefault(prefix_list[index],float(genome_equivalents))
#                                 print("The Average Genome Equivalent of file {} is {}".format(prefix_list[index],genome_equivalents))
#     return genome_equivalents_dic

def get_genome_equivalent(input_AGS, prefix_list):
    genome_equivalents_dic={}
    for index, j in enumerate(input_AGS): 
        for lines in open(j,'r'):
                            if lines.startswith('genome_equivalents'):
                                lines_set = lines.split('\n')[0].split('\t')
                                genome_equivalents = float(lines_set[1])
                                genome_equivalents_dic.setdefault(prefix_list[index],float(genome_equivalents))
                                print("The Average Genome size of file {} is {}".format(prefix_list[index],genome_equivalents))
    return genome_equivalents_dic

def get_cell_number(input_kk2):
    
    kraken=input_kk2
    for lines in open(kraken,'r'):
        content = lines.split('\n')[0].split('\t')
        if 'Bacteria' in lines:
            copy_cell = float(content[1]) # pair end
            break
    return copy_cell


def multi_info_sum():
    openthreads = len(sample_list) 
    exfiles = []
    for i in range(openthreads):
        worker = multiprocessing.Process(target=info_sum,args=([sample_list[i],genome_equivalents_dic[sample_list[i]],copy_cell_dic[sample_list[i]]]))
        worker.start()
        print("Now processing:{}".format(sample_list[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()  



if __name__ == "__main__":
    
    #假设输入文件为示例文件，放在for循环的开头第一层 for i = samle_name
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    sample_list= path.file_base_acquire(file_abs_path) #sample name without suffix .fa
    write_file=output + "/CompRanking_"+ project_prefix + "_Contigs_Risk_Summary_Celleq.txt"
    with open(write_file, "w") as f1:
        f1.write("sample_name/index\tnCelleq\tnContigs_num\tnARGs_contigs\tnMGEs_contig\tnMGEs_plasmid_contig\tnMGEs_phage_contigs\tnPAT_contigs\tnARGs_MGEs_contig\tnARGs_MGEs_plasmid_contigs\tnARGs_MGEs_phage_contigs\tnARGs_MGEs_PAT_contigs\tfARG\tfMGE\tfMGE_plasmid\tfMGE_phage\tfPAT\tfARG_MGE\tfARG_MGE_plasmid\tfARG_MGE_phage\tfARG_MGE_PAT\tscore_pathogenic\tscore_phage\tscore_plasmid")
    
    #load AGS
    input_AGS_dir=os.path.join(input_dir,project_prefix,
                            "CompRanking_intermediate/preprocessing/5M_contigs/AGS")
    file_list,prefix=path.getPrefix(input_AGS_dir)
    print(file_list,prefix)
    genome_equivalents_dic=get_genome_equivalent(file_list,prefix)
    print(genome_equivalents_dic)
    
    #load cell number
    copy_cell_dic={}
    for i in file_name_base:
        input_kk2=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_report_kk2_mpaStyle.txt")
        copy_cell=get_cell_number(input_kk2)
        copy_cell_dic.setdefault(i,float(copy_cell))
        
    
    
    
    #run
    multi_info_sum()
    