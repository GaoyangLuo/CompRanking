#!/usr/bin/env python
# title             :baseInfoExtra.py -> baseInfoExtra.ipynb
# description       :This script is for benchmarking using single copy gene (scg) copies 
#                    as normalization bases to generate risk score
# author            :Gaoyang Luo
# date              :202201119
# version           :1.0
# usage             :python baseInfoExtra.py -i <input_dir> -p <project_prefix>
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import numpy as np
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

def info_sum(sample_list, copy_scg):
# for i in sample_list: #获取sample_name
    sample_name=sample_list
    print(sample_name)
    # 获取结果表
    file_path=os.path.join(input_dir,project_prefix,"CompRanking_result/CompRanking_"+ sample_name +"_Summary.tsv")
    df=pd.read_csv(file_path, 
                    sep='\t', header=0)
    #ncontigs_counts
    nContigs_num=len(df.Contig.unique()) # numbers of scg copies
    #get genome equivalent
    # nContigs=float(copy_16S)/1550
    nContigs=copy_scg #scg as bases
    print("scg_copy base is {}".format(nContigs))
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
    hazard_point=float(0.01)
    distance_pathogenic = math.sqrt((hazard_point - fARG)**2 + (hazard_point - fARG_MGE)**2 + (hazard_point - fARG_MGE_PAT)**2)
    distance_phage=math.sqrt((hazard_point - fARG)**2 + (hazard_point - fARG_MGE_phage)**2 + (hazard_point - fARG_MGE_PAT)**2)
    distance_plasmid= math.sqrt((hazard_point - fARG)**2 + (hazard_point - fARG_MGE_plasmid)**2 + (hazard_point - fARG_MGE_PAT)**2)
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

# def get_genome_equivalent(input_AGS, prefix_list):
#     genome_equivalents_dic={}
#     for index, j in enumerate(input_AGS): 
#         for lines in open(j,'r'):
#                             if lines.startswith('average_genome_size'):
#                                 lines_set = lines.split('\n')[0].split('\t')
#                                 genome_equivalents = float(lines_set[1])
#                                 genome_equivalents_dic.setdefault(prefix_list[index],float(genome_equivalents))
#                                 print("The Average Genome size of file {} is {}".format(prefix_list[index],genome_equivalents))
#     return genome_equivalents_dic

def multi_info_sum():
    openthreads = len(sample_list) 
    exfiles = []
    for i in range(openthreads):
        input_scg=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/cov", 
                                    sample_list[i]+"_5M_contigs_scg_Protein_dimond.txt") 
        input_rpkm=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/cov",
                                    sample_list[i]+"_5M_contigs_gene.rpkm")
        input_indexFile=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    sample_list[i]+"_5M_contigs.fna2faa.index")
        # metagenomes
        df_scg=pd.read_csv(input_scg, sep="\t",header=None)
        df_scg.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        # filter out contigs identity under 30
        df_scg_iden50 = df_scg[df_scg.identity > 30]
        # filter out contigs length larger than 0
        df_scg_len = df_scg_iden50[df_scg_iden50.alignLen > 0]
        # filter out bit larger than 50
        df_scg_bit50 = df_scg_len[df_scg_len.bit > 50]
        
        #load rpkm alinged reads number and mapped reads number
        # input_rpkm="/lomi_home/gaoyang/software/CompRanking/tmp_DSR/DSR/CompRanking_intermediate/preprocessing/5M_contigs/cov/S0PCL_clean.sorted_filtered.rpkm"
        #get reads numbers
        for lines in open(input_rpkm,'r'):
            #content = lines.split('\n')[0].split('\t')
            content = lines.split('\t')
            if '#Mapped' in lines:
                mapped_reads = float(content[1]) # pair end total num of mapped reads
                break
        #load rpkm table
        df_rpkm=pd.read_csv(input_rpkm,sep="\t",header=4)
        array=np.array(df_rpkm)
        array=array.tolist()
        #write aligned reads into dic
        """
        refer rpkm_dic by using contig_orf
        """
        rpkm_dic={}
        for j in array:
            rpkm_dic.setdefault(str(j[0]),float(j[4]))
        #calculate mrna length
        gene_Length_dic={}
        for j in array:
            gene_Length_dic.setdefault(str(j[0]),float(j[1]))
        #change orf id to contigs id
        df_index=pd.read_csv(input_indexFile,sep="\t",header=None)
        array_indexFile=np.array(df_index)
        array_indexFile=array_indexFile.tolist()
        index_dic={}
        for j in array_indexFile:
            index_dic.setdefault(str(j[1]),str(j[0]))
        #calculate scg base for cell number normalization 
        array_scg=np.array(df_scg_bit50)
        array_scg=array_scg.tolist()
        scg_sum=0
        for j in array_scg:
            print(j[0])
            contig_orf=index_dic[j[0]]
            scg_sum+=(rpkm_dic[contig_orf] / gene_Length_dic[contig_orf])
        num_scg=float(scg_sum / 16 * 10000)
        # print("The sample ID {0} and number of scg {1} are under calculated...".format(sample_list[i],num_scg))
        
        worker = multiprocessing.Process(target=info_sum,args=([sample_list[i],num_scg]))
        worker.start()
        print("Now processing:{}".format(sample_list[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()  



if __name__ == "__main__":
    
    #假设输入文件为示例文件，放在for循环的开头第一层 for i = samle_name
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    sample_list= path.file_base_acquire(file_abs_path) #sample name without suffix .fa
    write_file=output + "/CompRanking_"+ project_prefix + "_Contigs_Risk_Summary_fna2faa_scg.txt"
    with open(write_file, "w") as f1:
        f1.write("sample_name/index\tnSCG_Base\tnContigs_num\tnARGs_contigs\tnMGEs_contig\tnMGEs_plasmid_contig\tnMGEs_phage_contigs\tnPAT_contigs\tnARGs_MGEs_contig\tnARGs_MGEs_plasmid_contigs\tnARGs_MGEs_phage_contigs\tnARGs_MGEs_PAT_contigs\tfARG\tfMGE\tfMGE_plasmid\tfMGE_phage\tfPAT\tfARG_MGE\tfARG_MGE_plasmid\tfARG_MGE_phage\tfARG_MGE_PAT\tscore_pathogenic\tscore_phage\tscore_plasmid")
    
    # #load AGS
    # input_AGS_dir=os.path.join(input_dir,project_prefix,
    #                         "CompRanking_intermediate/preprocessing/5M_contigs/AGS")
    # file_list,prefix=path.getPrefix(input_AGS_dir)
    # print(file_list,prefix)
    # genome_equivalents_dic=get_genome_equivalent(file_list,prefix)
    # print(genome_equivalents_dic)
    
    
    
    #run
    multi_info_sum()
    