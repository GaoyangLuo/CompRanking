#!/usr/env/bin python
#coding :UTF-8
# title             :cpr_run.py
# description       :run the whole process
# author            :Gaoyang Luo
# date              :20220901
# version           :1.0
# usage             :python cpr_run.py -i <input_dir> 
#                                      -n <input_file_name>
#                                      -p <project_name default:[CompRanking]> 
#                                     
#example            :python cpr_run.py -i /lomi_home/gaoyang/software/CompRanking/test -p CompRanking -r 0 -t 2
# required packages :optparse, subprocess 
# notification: makesure that put PathoFact table in to ./<Project_name>_PathoFact_result 
# and put seeker result into ./<Project_name>_seeker_result
#==============================================================================
import subprocess
import os
import glob
import sys
import optparse
import datetime
import multiprocessing
import threading
import yaml
sys.path.append("..")
from compranking import ARG_ranker
from compranking import path, summary_all, MOB_concat,Virulence_processing
from compranking.AMR_combination import AMRCombined

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-n", "--nameoffile", action = "store", type = "string", dest = "input_file", 
                  help = "name of input file")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")
(options, args) = parser.parse_args()

####################################preseting#################################### 
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
input_file=options.input_file

###################Combine AMR########################
#gloab settings
input_cpr_VF_sum="../databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" #fixed
input_cpr_VF_sum=os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" )#fixed
output=os.path.join(input_dir,project_prefix,"CompRanking_result")

AMR_combine=AMRCombined()

start_sum = datetime.datetime.now() #time start 

i=input_file

start_all = datetime.datetime.now()
#if file large it will be very slow
#set ARG&MGE input
input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs_DeepARG.out.mapping.ARG")
input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein60_Result.tsv")
input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
input_dvf=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DVF",i+"_5M_contigs.fa_gt500bp_dvfpred.txt")
input_plasflow=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Plasflow",i+"_5M_contigs_plasflow_predictions.tsv")
seeker_table=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Seeker","seeker_"+i+"_5M_contigs_output.txt")
input_mobileOG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/MobileOG",i+"_5M_contigs_mobileOG_diamond.txt")
input_mob_conj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_Conj_plasmids_id_out")
input_mob_unconj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_mob_unconj_plasmids_id_out")
input_def=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DEF",i+"_5M_contigs.fa_pred_one-hot_hybrid.txt")

#set VF input
input_contig=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
input_ERR_VFDB_output=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/VFDB",i+"_5M_contigs_VFDB.out")
input_patric=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/PATRIC",i+"_5M_contigs_PATRIC.out")

#generate sum table without mob ref
df_AMR_annotate_contig=AMR_combine.AMR_combined(input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_def,seeker_table,input_mobileOG)
#generate sum table with mob ref
df_AMR_annotate_MOB_contig=MOB_concat.plasMOB_concat(input_mob_conj,input_mob_unconj,df_AMR_annotate_contig)
df_AMR_annotate_MOB_refFilter_contig=AMR_combine.refFilter(df_AMR_annotate_MOB_contig)
#save ARG to tsv
df_AMR_annotate_MOB_refFilter_contig.to_csv(output + "/CompRanking_" + i + "_AMR_MOB_prediction.tsv", sep="\t", index=None)

#vf processing
df_VFs_PATH_contig=Virulence_processing.VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,input_patric)
#save vf to tsv
df_VFs_PATH_contig.to_csv(output + "/CompRanking_"+ i +"_Virulence_Pathogenic_prediction.tsv",sep="\t",index=None)

#summary
df_sum=summary_all.sum_all(df_AMR_annotate_MOB_refFilter_contig,df_VFs_PATH_contig)
df_sum.to_csv(output + "/CompRanking_" + i + "_Summary.tsv", sep="\t", index=None)
    
end_sum = datetime.datetime.now() #time end
end_all = datetime.datetime.now()

print("Summary ouput takes time {}: ".format(end_sum-start_sum))
print("All process takes time {}: ".format(end_all-start_all))
    