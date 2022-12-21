#!/usr/bin/env python
# title             :baseInfoExtra.py -> risk_cal.ipynb
# description       :calculate relative abundance of genes
# author            :Gaoyang Luo
# date              :202201119
# version           :1.0
# usage             :python risk_cal.py -i <input_dir>
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
import path
import math
import optparse

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
threads=options.threads
database=options.database
output=os.path.join(input_dir,"CompRanking/CompRanking_result")
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.output_dir is None):
    output = os.path.join(input_dir,"CompRanking/CompRanking_result") #default output directory

    
#假设输入文件为示例文件，放在for循环的开头第一层 for i = samle_name
file_abs_path=path.file_abs_path_list_generation(input_dir)
sample_list= path.file_base_acquire(file_abs_path)
write_file="CompRanking_"+ project_prefix + "_Contigs_Risk_Summary.txt"
with open(write_file, "w") as f1:
    f1.write("sample_name/index\tnContigs\tnARGs_contigs\tnMGEs_contig\tnMGEs_plasmid_contig\tnMGEs_phage_contigs\tnPAT_pathogenic_contigs\tnPAT_all_contigs\tnARGs_MGEs_contig\tnARGs_MGEs_plasmid_contigs\tnARGs_MGEs_phage_contigs\tnARGs_MGEs_PAT_contigs\tfARG\tfMGE\tfMGE_plasmid\tfMGE_phage\tfPAT_pathogenic\tfARG_MGE\tfARG_MGE_plasmid\tfARG_MGE_phage\tfARG_MGE_PAT_pathogenic\tscore_pathogenic\tscore_phage\tscore_plasmid")

for i in sample_list: #获取sample_name
    sample_name=i
    print(sample_name)
    # 获取结果表
    file_path=os.path.join(input_dir,project_prefix,"CompRanking_"+ sample_name +"_Summary.tsv")
    df=pd.read_csv(file_path, 
                    sep='\t', header=0)
    df=df.drop(["Unnamed: 0"],axis=1)
    #ncontigs_counts
    nContigs=len(df.Contig_ID.unique()) # numbers of contigs
    #nARGs_contigs_counts
    #RankI_Risk
    #ARGs_List&&ARGs_typeCount
    ARGs_List=[]
    arg_filter=df[df.ARG!="-"]
    nARGs_contigs=len(arg_filter.Contig.unique()) # numbers of contigs containing predicted ARGs
    #nMGEs_counts
    MGEs_filter=df[df.MGE_prediction!="unclassified" ]
    MGEs_filter=MGEs_filter[MGEs_filter.MGE_prediction!="chromosome"]
    MGEs_filter=MGEs_filter[MGEs_filter.MGE_prediction!="ambiguous (phage/chromosome)"]
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
    #nPAT_all_counts
    allPAT_filter=df[df.Pathogenicity!="-"]
    nPAT_all_contigs=len(allPAT_filter.Contig.unique())
    #nARGs_MGEs_counts
    #RankII_Risk
    nARGs_MGEs_filter=arg_filter[arg_filter.Pathogenicity!="-"]
    nARGs_MGEs_filter=nARGs_MGEs_filter[nARGs_MGEs_filter.peram_MGE_prediction!="chromosome"]
    nARGs_MGEs_filter=nARGs_MGEs_filter[nARGs_MGEs_filter.peram_MGE_prediction!="ambiguous (phage/chromosome)"]
    nARGs_MGEs_contigs=len(nARGs_MGEs_filter.Contig.unique())
    #nARGs_MGEs_plasmid_counts
    nARGs_MGEs_plasmid_filter=arg_filter[arg_filter.peram_MGE_prediction=="plasmid"]
    nARGs_MGEs_plasmid_contigs=len(nARGs_MGEs_plasmid_filter.Contig.unique())
    #nARGs_MGEs_phage_counts
    nARGs_MGEs_phage_filter=arg_filter[arg_filter.peram_MGE_prediction=="phage"]
    nARGs_MGEs_phage_contigs=len(nARGs_MGEs_phage_filter.Contig.unique())
    #nARGs_MGEs_allPAT_contigs
    nARGs_MGEs_allPAT_filter=nARGs_MGEs_filter[nARGs_MGEs_filter.Pathogenicity !="-"]
    nARGs_MGEs_allPAT_contigs=len(nARGs_MGEs_allPAT_filter.Contig.unique())
    #nARGs_MGEs_PAT(pathogenic)
    #RankIII_Risk
    nARGs_MGEs_PAT_filter=nARGs_MGEs_filter[nARGs_MGEs_filter.Pathogenicity=="pathogenic"]
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
    distance_pathogenic = math.sqrt((0.01 - fARG)**2 + (0.01 - fARG_MGE)**2 + (0.01 - fARG_MGE_PAT)**2)
    distance_phage=math.sqrt((0.01 - fARG)**2 + (0.01 - fARG_MGE_phage)**2 + (0.01 - fARG_MGE_PAT)**2)
    distance_plasmid= math.sqrt((0.01 - fARG)**2 + (0.01 - fARG_MGE_plasmid)**2 + (0.01 - fARG_MGE_PAT)**2)
    #calculate risk score
    # score_all = 1.0 / ( (2 + math.log10(distance_all))**2 )
    score_pathogenic= 1.0 / ( (2 + math.log10(distance_pathogenic))**2 )
    score_phage= 1.0 / ( (2 + math.log10(distance_phage))**2 )
    score_plasmid= 1.0 / ( (2 + math.log10(distance_plasmid))**2 )

    result=[nContigs, nARGs_contigs, nMGEs_contigs,nMGEs_plasmid_contigs,nMGEs_phage_contigs, nPAT_pathogenic_contigs,nPAT_all_contigs, nARGs_MGEs_contigs,nARGs_MGEs_plasmid_contigs, nARGs_MGEs_phage_contigs,nARGs_MGEs_PAT_contigs, fARG, fMGE,fMGE_plasmid,fMGE_phage, fPAT_pathogenic, fARG_MGE, fARG_MGE_plasmid,fARG_MGE_phage,fARG_MGE_PAT,score_pathogenic,score_phage,score_plasmid]
    output="\t".join(map(str, result))

    with open(write_file, "a") as f:
        f.write("\n"+ sample_name + "\t" + output)