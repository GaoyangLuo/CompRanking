#!/usr/env/bin python
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
#parser.add_option("-o", "--out", action = "store", type = "string", dest = "output_dir",
									#default='./', help = "output directory")

(options, args) = parser.parse_args()

####################################preseting#################################### 
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
output_dir=options.output_dir
threads=options.threads
config_path=options.config_file


#scrip path
PREPROCESSING="./scripts/preprocessing_run.sh"
AMR_PREDICTION="./scripts/AMR_run.sh"
MGE_PREDICTION="./scripts/MGE_run.sh"
MGE2_PREDICTION="./scripts/MGE_run_sec_round.sh"
VIRULENCE_PREDICTION="./scripts/Virulence_run.sh"
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.threads is None):
    threads = "12" #default threads

#===============================================================================
####################################Function####################################
#configeration
def read_conda_path(name,path,yaml_path):
    CONDA_PATH=[]
    #open yaml file
    with open(yaml_path,"r") as f:
        data=yaml.load(f,Loader=yaml.FullLoader)
        # print(data)
        try:
            if name in data.keys():
                if len(data[name][path]) > 0:
                    # print(data[name][path])
                    CONDA_PATH.append(data[name][path])
                    # print(CONDA_PATH)
                else:
                    raise Exception("Empty or wrong conda bin path")
        except:
            raise Exception("Wrong configeration, \
                   please firstly write your conda bin path into the configeratoion file")
        if len(CONDA_PATH) != 0:
            try:
                for i in CONDA_PATH:
                    if os.path.exists(i):
                        # print("Conda bin path exist")
                        return CONDA_PATH
                    else:
                        if i.startwith("/"):
                            print("plese check your conda path")
                        else:
                            print("please use absolute path started with / ")
            except:
                raise IndexError("Conda bin path don't existed")
                sys.exit('Something happened')

#acquire abs path of sample
def file_abs_path_list_generation(input):
    input=input_dir+"/*fa"
    file_abs_path_list=glob.glob(input)
    return file_abs_path_list

#acquire sample name base
def file_base_acquire(file_abs_path): #get names of all the sample files
    base_list=[]
    #Acquire file_base
    try:
        for i in file_abs_path:
            fileBase_1=i.split("/")[-1]
            fileBase_2=fileBase_1.strip(".fa")
            base_list.append(fileBase_2)
    except:
        raise IndexError("Faile to acquire base name of samples...")
    return base_list

#mutiprocessing
##AMR prediction
def ARG_prediction():
    subprocess.call(["bash", AMR_PREDICTION, 
                     "-i", rgi_input, "-t", threads, "-p", project_prefix, "-c", conda_path_str])    
##MGE prediction
def MGE_prediction():
    subprocess.call(["bash", MGE_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-c", conda_path_str])
##MGE2 prediction
def MGE2_prediction():
    subprocess.call(["bash", MGE2_PREDICTION, 
                     "-i", input_dir, "-t", threads, "-p", project_prefix, "-c", conda_path_str])  
##Virlence prediction
def VIR_prediction():
    subprocess.call(["bash", VIRULENCE_PREDICTION, 
                     "-i", VF_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])

#===============================================================================
####################################Get Started#################################
if __name__ == '__main__': 
    #### Step 0 Presetting ####
    yaml_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),config_path)
    rgi_input=os.path.join(input_dir,project_prefix,"preprocessing/5M_contigs")
    VF_input=os.path.join(input_dir,project_prefix,"preprocessing/5M_contigs")
    print(yaml_path)
    #multiprocessing
    AMR_PRED = multiprocessing.Process(target=ARG_prediction)
    MGE_PRED = multiprocessing.Process(target=MGE_prediction)
    MGE2_PRED = multiprocessing.Process(target=MGE2_prediction)
    VIR_PRED = multiprocessing.Process(target=VIR_prediction)
    
    #Write in abs conda path
    name="CompRanking"
    path="abs_path_to_conda_bin"
    conda_path_str="".join(read_conda_path(name,path,yaml_path)) #record abs path of conda bin
    print("The absolute path to conda bin is:{0}".format(conda_path_str)) 
    # subprocess.call(["bash", VIRULENCE_PREDICTION, 
    #                  "-i", VF_input, "-t", threads, "-p", project_prefix, "-m", conda_path_str])

    
    
    #### Step 1 Preprocessing ####
    # start = datetime.datetime.now() #time start
    # subprocess.call(["bash", PREPROCESSING, "-i", input_dir]) 
    # end = datetime.datetime.now() #time end
    # print(end-start)
    
    
    
    # AMR_PRED.start()
    VIR_PRED.start()
    # AMR_PRED.join()
    VIR_PRED.join()
    # MGE_PRED.start()
    # MGE2_PRED.start()
    # MGE_PRED.join()
    # MGE2_PRED.join()
    #### Step 2 ARG Prediction ####
    # start = datetime.datetime.now() #time start
    # # rgi_input=os.path.join(input_dir,project_prefix,"preprocessing/5M_contigs")
    # # print(rgi_input)
    # # subprocess.call(["bash", AMR_PREDICTION, "-i", rgi_input, "-t", threads])
    # ARG_prediction()
    # end = datetime.datetime.now() #time end
    # print(end-start)

    #### MGE prediction ####
    # start = datetime.datetime.now() #time start
    # subprocess.call(["bash", MGE_PREDICTION, "-i", input_dir, "-t", threads])
    # end = datetime.datetime.now() #time end
    # print(end-start)

    #### VF prediction ####
    # start = datetime.datetime.now() #time start
    # VF_input=os.path.join(input_dir,project_prefix,"preprocessing/5M_contigs")
    # subprocess.call(["bash", VIRULENCE_PREDICTION, "-i", VF_input, "-t", threads])
    # end = datetime.datetime.now() #time end
    # #print(end-start)
    #print(type(input_dir))
    #print(type(project_prefix))
    #print(type("/preprocessing/5M_contigs"))
    #print(VF_input)
    
    
    #check file completeness
    # file_list=file_abs_path_list_generation(input_dir)
    # base_list=file_base_acquire(file_list)
    # print(file_list)
    # print(base_list)
    # yt.check_file_completness()





