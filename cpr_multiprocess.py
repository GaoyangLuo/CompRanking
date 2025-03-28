#!/usr/env/bin python
#coding :UTF-8
# title             :cpr_run.py
# description       :run the whole process
# author            :Gaoyang Luo
# date              :20220901
# version           :1.0
# usage             :python cpr_run.py -i <input_dir> 
#                                      -o <output_dir default: [input_dir]> 
#                                      -t <threads defalut: 24> 
#                                      -p <project_name default:[CompRanking]> 
#                                      -c <config_file default:[config.yaml]>
#                                      -r <default: 1 (restart), 0 means don't restat>
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
import time
sys.path.append("..")
from compranking import ARG_ranker
from compranking import path, summary_all, MOB_concat,Virulence_processing
from compranking.AMR_combination import AMRCombined

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_dir",
				 help = "director tontained output files")
parser.add_option("-t", "--threads", action = "store", type = "string", dest = "threads",
				 help = "how many cpus you want use")      
parser.add_option("-c", "--config_file", action = "store", type = "string", dest = "config_file", 
                  help = "file contains basic configeration information")           
parser.add_option("-r", "--restart", action = "store", type = "string", dest = "restart",
									default='1', help = "restart all the processs")
parser.add_option("-s", "--splitLength", action = "store", type = "string", dest = "splitLength",
									default='10000', help = "split length of the input")
parser.add_option("-l", "--filterLength", action = "store", type = "string", dest = "filterLength",
									default='500', help = "split length of the input")

(options, args) = parser.parse_args()

####################################preseting#################################### 
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
output_dir=options.output_dir
threads=options.threads
config_path=options.config_file
restart=options.restart
splitLength=options.splitLength
filterLength=options.filterLength
#scrip path
PREPROCESSING="./scripts/preprocessing_run.sh"
AMR1_PREDICTION="./scripts/RGI_run.sh"
AMR2_PREDICTION="./scripts/DeepARG_run.sh"
AMR3_PREDICTION="./scripts/SARG_run.sh"
MGE1_PREDICTION="./scripts/Plasflow_run.sh"
MGE2_PREDICTION="./scripts/DVF_Seeker_run.sh"
MGE3_PREDICTION="./scripts/mobileOG_run.sh"
MGE4_PREDICTION="./scripts/seeker_separate.sh"
MGE5_PREDICTION="./scripts/def.sh"
VIRULENCE_PREDICTION="./scripts/Virulence_run.sh"
PLASCAD_PREDICTION="./scripts/plascad_run.sh"
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.threads is None):
    threads = "32" #default threads
if (options.config_file is None):
    config_path = "./test_yaml.yaml"#default config_file path
if options.restart == "1":
    os.system("rm checkdone/" + project_prefix + "*done")
if (options.splitLength is None):
    options.splitLength == "10000"
    

#===============================================================================
####################################Function####################################

#MULTIPROCESSING DEFINITION
##AMR PREDICTION
def ARG1_prediction(): #rgi
    subprocess.call(["bash", AMR1_PREDICTION, 
                     "-i", rgi_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])
def ARG2_prediction():#deeparg
    subprocess.call(["bash", AMR2_PREDICTION, 
                     "-i", deeparg_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])   
def ARG3_prediction():#sarg
    subprocess.call(["bash", AMR3_PREDICTION, 
                    "-i", sarg_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])     
##MGE prediction
def MGE1_prediction(): #plasflow
    subprocess.call(["bash", MGE1_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str])
def MGE2_prediction(): #dvf
    subprocess.call(["bash", MGE2_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str, "-l", filterLength])  
def MGE3_prediction(): #mobileOG
    subprocess.call(["bash", MGE3_PREDICTION, 
                     "-i", mobileog_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])
def MGE4_prediction(): #seeker
    subprocess.call(["bash", MGE4_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str]) 
def MGE5_prediction(): #def
    subprocess.call(["bash", MGE5_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str])     
##Virlence prediction
def VIR_prediction(): #VF&Pathogen
    subprocess.call(["bash", VIRULENCE_PREDICTION, 
                     "-i", VF_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])
#plascad
def plascad_prediction(): #plascad   
    subprocess.call(["bash", PLASCAD_PREDICTION, 
                     "-i", plascad_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])     

##RESULT RPOCESS
#combine results
def AMR_process(file_name_base):
    i=file_name_base
    #set ARG&MGE input
    # input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
    input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.fna2faa.RGI.out.txt")
    # input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs_DeepARG.out.mapping.ARG")
    input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs.fna2faa_DeepARG.out.mapping.ARG")
    input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein60_Result.tsv")
    # input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
    input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.fna2faa.index")
    input_dvf=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DVF",i+"_5M_contigs.fa_gt"+str(filterLength)+"bp_dvfpred.txt")
    # input_plasflow=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Plasflow",i+"_5M_contigs_plasflow_predictions.tsv")
    seeker_table=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Seeker","seeker_"+i+"_5M_contigs_output.txt")
    input_mobileOG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/MobileOG",i+"_5M_contigs_mobileOG_diamond.txt")
    input_mob_conj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_Conj_plasmids_id_out")
    input_mob_unconj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_mob_unconj_plasmids_id_out")
    input_def=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DEF",i+"_5M_contigs.fa_pred_onehot_hybrid.tsv")
    
    #set VF input
    # input_contig=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
    input_contig=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.fna2faa.index")
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

#threads combine
def AMR_combine_main():
    openthreads = len(file_name_base) 
    exfiles = []
    for i in range(openthreads):
        # worker = threading.Thread(target=AMR_process,args=([file_name_base[i]]))
        worker = multiprocessing.Process(target=AMR_process,args=([file_name_base[i]]))
        worker.start()
        print("正在处理:{}".format(file_name_base[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()
#===============================================================================
####################################Get Started#################################
if __name__ == '__main__': 
    #### Step 0 Presetting ####
    yaml_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),config_path)
    faaFile_input=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs")
    rgi_input=faaFile_input
    VF_input=faaFile_input
    mobileog_input=faaFile_input
    deeparg_input=faaFile_input
    sarg_input=faaFile_input
    plascad_input= input_dir #os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/ori_file")
    output_prefix=os.path.join(input_dir,project_prefix,"CompRanking_result")
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    print(yaml_path)

    #multiprocessing
    AMR_PRED1 = multiprocessing.Process(target=ARG1_prediction) #rgi
    AMR_PRED2 = multiprocessing.Process(target=ARG2_prediction) #deeparg
    AMR_PRED3 = multiprocessing.Process(target=ARG3_prediction) #sarg
    # MGE1_PRED = multiprocessing.Process(target=MGE1_prediction) #plasflow
    MGE2_PRED = multiprocessing.Process(target=MGE2_prediction) #dvf
    MGE3_PRED = multiprocessing.Process(target=MGE3_prediction) #mobileOG
    MGE4_PRED = multiprocessing.Process(target=MGE4_prediction) #seeker
    MGE5_PRED = multiprocessing.Process(target=MGE5_prediction) #def
    PLASCAD_PRED = multiprocessing.Process(target=plascad_prediction) #plascad
    VIR_PRED = multiprocessing.Process(target=VIR_prediction) #VFDB&PATH
    
    #Write in abs conda path
    path_bin="abs_path_to_conda_bin"
    conda_path_str="".join(path.read_conda_path("CompRanking",path_bin,yaml_path)) #record abs path of conda bin
    print("The absolute path to conda bin is:{0}".format(conda_path_str)) 
    
    ################################### Step 1 Preprocessing ################################
    start_all = datetime.datetime.now() 
    start_prepro = datetime.datetime.now() #time start
    print("Preprocessing fasta files...")
    subprocess.call(["bash", PREPROCESSING, "-i", input_dir, "-p", project_prefix, "-t", threads,"-l",filterLength]) 
    #split files
    # for i in file_name_base:
    #     file_abs_pth_dir=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs")
    #     input_splitfile=os.path.join(file_abs_pth_dir,i +"_5M_contigs.fa")
    #     output_splitfile=input_dir+"/"+project_prefix+"/CompRanking_intermediate/preprocessing/5M_contigs/split/"+i+"_5M_contigs"
    #     print(output_splitfile)
        
    #     try:
    #         os.stat(output_splitfile)
    #     except:
    #         os.system("mkdir -p " + output_splitfile)
    #     os.system("python split.py "+ input_splitfile + " " + splitLength + " " + output_splitfile)
        
    end_prepro = datetime.datetime.now() #time end
    print("Step Preprocessing cost time: {}".format(end_prepro-start_prepro))
    
    ################################## Step 2 multiprocessing ############################################
    #ARG search
    #time start
    start = datetime.datetime.now() 
    AMR_PRED1.start()
    AMR_PRED2.start()
    AMR_PRED3.start()
    AMR_PRED1.join() #start rgi 
    AMR_PRED2.join() #start deeparg 
    AMR_PRED3.join() #start sarg
    #time end
    end = datetime.datetime.now() 
    print("ARG search cost time: {}".format(end-start))
    
    #### VF prediction ####
    start_VF = datetime.datetime.now() #time start
    VIR_PRED.start()
    VIR_PRED.join() #start vf
    end_VF = datetime.datetime.now() #time end
    print("VF and Pathogen prediction cost: {}".format(end_VF-start_VF))
    
    #### MGE prediction ####
    start_MGE = datetime.datetime.now() #time start
    PLASCAD_PRED.start() #start mobplas
    PLASCAD_PRED.join() 
    # MGE4_PRED.start() #start seeker
    MGE3_PRED.start() #start mog
    MGE5_PRED.start() #start def
    MGE5_PRED.join()
    MGE3_PRED.join()
  
    
    
    # MGE1_PRED.start()
    MGE2_PRED.start() #dvf
    MGE2_PRED.join()
    
    # MGE1_PRED.join()
    time.sleep(0.5)
    
    #run phage prediction
    MGE4_prediction()
    time.sleep(0.5)
    
    end_MGE = datetime.datetime.now() #time end
    print("MGE prediction cost: {}".format(end_MGE-start_MGE))
    end_all = datetime.datetime.now() #time end
    print("All prediction cost time: {}".format(end_all-start_all))
    
    ###################check output########################
    #check file completeness
    # file_list=file_abs_path_list_generation(input_dir)
    # base_list=file_base_acquire(file_list)
    # print(file_list)
    # print(base_list)
    # yt.check_file_completness()
    
    ###################rankARG########################
    ##fixed settings
    input_argrank="databases/SARG/ARG_rank.txt"
    input_sarg_length="databases/SARG/SARG.db.fasta.length"
    input_sarg_structure="databases/SARG/SARG.structure.txt" #"/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.structure.txt"
    SARG_output=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking") #fixed path
    file_abs_path=path.file_abs_path_list_generation(input_dir) #fixed path
    file_name_base = path.file_base_acquire(file_abs_path) #fixed path
    
    # #arg rank processing
    # for i in file_name_base:
    #     input_sarg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking", i+"_5M_contigs_SARG_Protein_diamond.txt")
    #     print(i)
    #     if not os.path.getsize(input_sarg) != 0:
    #         # print(f"[DEBUG] Running ARG_ranker on {input_sarg} ...")
    #         ARG_ranker.arg_rank(input_sarg, input_sarg_length,input_sarg_structure, input_argrank,i, SARG_output)
    #         # print(f"[DEBUG] Finished ARG_ranker on {input_sarg}.")
    #     else:
    #         print(f"[WARNING] Skipping {input_sarg} because it does not exist or is empty.")
    #         continue
    
    
    # ARG rank processing
for i in file_name_base:
    input_sarg = os.path.join(input_dir, project_prefix, "CompRanking_intermediate/AMR/ARGranking", i + "_5M_contigs_SARG_Protein_diamond.txt")
    print(f"[DEBUG] Processing: {input_sarg}")

    # 先检查文件是否存在
    if os.path.exists(input_sarg):
        try:
            # 尝试打开文件并读取第一行，确保不是空文件
            with open(input_sarg, "r") as f:
                first_line = f.readline().strip()

            if first_line:  # 如果第一行有内容，则执行 ARG_ranker
                print(f"[DEBUG] Running ARG_ranker on {input_sarg} ...")
                ARG_ranker.arg_rank(input_sarg, input_sarg_length, input_sarg_structure, input_argrank, i, SARG_output)
                print(f"[DEBUG] Finished ARG_ranker on {input_sarg}.")
            else:
                print(f"[WARNING] Skipping {input_sarg} because it is empty (no content).")
        
        except Exception as e:
            print(f"[ERROR] Failed to read {input_sarg}: {e}")
    else:
        print(f"[ERROR] {input_sarg} does not exist. Skipping ARG ranking.")
        
        
    ###################Combine AMR########################
    #gloab settings
    input_cpr_VF_sum="../databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" #fixed
    input_cpr_VF_sum=os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" )#fixed
    output=os.path.join(input_dir,project_prefix,"CompRanking_result")
    
    AMR_combine=AMRCombined()
    
    start_sum = datetime.datetime.now() #time start 
    AMR_combine_main()
    # for i in file_name_base: #if file large it will be very slow
    #     #set ARG&MGE input
    #     input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
    #     input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs_DeepARG.out.mapping.ARG")
    #     input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein60_Result.tsv")
    #     input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
    #     input_dvf=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DVF",i+"_5M_contigs.fa_gt"+str(filterLength)+"bp_dvfpred.txt")
    #     input_plasflow=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Plasflow",i+"_5M_contigs_plasflow_predictions.tsv")
    #     seeker_table=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Seeker","seeker_"+i+"_5M_contigs_output.txt")
    #     input_mobileOG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/MobileOG",i+"_5M_contigs_mobileOG_diamond.txt")
    #     input_mob_conj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_Conj_plasmids_id_out")
    #     input_mob_unconj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_mob_unconj_plasmids_id_out")
    #     input_def=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DEF",i+"_5M_contigs.fa_pred_one-hot_hybrid.txt")
        
    #     #set VF input
    #     input_contig=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
    #     input_ERR_VFDB_output=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/VFDB",i+"_5M_contigs_VFDB.out")
    #     input_patric=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/PATRIC",i+"_5M_contigs_PATRIC.out")
        
    #     #generate sum table without mob ref
    #     df_AMR_annotate_contig=AMR_combine.AMR_combined(input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_def,seeker_table,input_mobileOG)
    #     #generate sum table with mob ref
    #     df_AMR_annotate_MOB_contig=MOB_concat.plasMOB_concat(input_mob_conj,input_mob_unconj,df_AMR_annotate_contig)
    #     df_AMR_annotate_MOB_refFilter_contig=AMR_combine.refFilter(df_AMR_annotate_MOB_contig)
    #     #save ARG to tsv
    #     df_AMR_annotate_MOB_refFilter_contig.to_csv(output + "/CompRanking_" + i + "_AMR_MOB_prediction.tsv", sep="\t", index=None)
        
    #     #vf processing
    #     df_VFs_PATH_contig=Virulence_processing.VF_processing(input_contig, input_ERR_VFDB_output,input_cpr_VF_sum,input_patric)
    #     #save vf to tsv
    #     df_VFs_PATH_contig.to_csv(output + "/CompRanking_"+ i +"_Virulence_Pathogenic_prediction.tsv",sep="\t",index=None)

    #     #summary
    #     df_sum=summary_all.sum_all(df_AMR_annotate_MOB_refFilter_contig,df_VFs_PATH_contig)
    #     df_sum.to_csv(output + "/CompRanking_" + i + "_Summary.tsv", sep="\t", index=None)
        
    end_sum = datetime.datetime.now() #time end
    end_all = datetime.datetime.now()
    
    print("Summary ouput takes time {}: ".format(end_sum-start_sum))
    print("All process takes time {}: ".format(end_all-start_all))
    
    
    
    
    
    
    
    
      #hmm processing
        #acquire base
        #input_dir="/lomi_home/gaoyang/software/CompRanking/test"
        # file_abs_path=path.file_abs_path_list_generation(input_dir)
        # file_name_base = path.file_base_acquire(file_abs_path)
        
        ##generate hmm.csv
        # input="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_itermediate/Virulence/111.hmm.txt"
        # suffix="_5M_contigs"
        # for i in file_name_base:
        #     hmm_file=os.path.join(input_dir,project_prefix,"CompRanking_itermediate","Virulence",i+suffix+"_tmp_hmm.txt")
        #     output_file=os.path.join(input_dir,project_prefix,"CompRanking_itermediate","Virulence","hmm_result",i+suffix+"_hmm.csv")
        #     hmm_processing.change_tab_hmmscan(input_hmmscan=hmm_file,output_hmm_csv=output_file)
        
        # ##hmmcsv2predcsv
        # positive_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/models_and_domains/positive_domains.tsv") 
        # for i in file_name_base:
        #     input=os.path.join(input_dir,project_prefix,"CompRanking_itermediate","Virulence","hmm_result",i+suffix+"_hmm.csv")
        #     output=os.path.join(input_dir,project_prefix,"CompRanking_itermediate","Virulence","hmm_result",i+suffix+"_VF_Prediction.csv")
        #     hmm_processing.VF_predition(input_hmmcsv=input,positive_domains=positive_path,name_VF_output=output)
    
    
        
        
    
    






