#!/usr/env/bin python
#encoding: utf-8
import yaml
import os
import optparse
import sys
import glob
import path

# parser = optparse.OptionParser()
# parser.add_option("-c", "--config_file", action = "store", type = "string", dest = "config_file", 
#                   help = "file contains basic configeration information")
# parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
#                   help = "director contained input fasta files")

# (options, args) = parser.parse_args()

# config_path=options.config_file
# input_dir=options.input_dir

# yaml_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),config_path)
# print("The yaml file path is:" + yaml_path)
class checkPoints():
    def __init__(self, workdir, preprocessing_suffix, RGI_suffix, DeepARG_suffix, 
                 SARG_suffix, DVF_suffix, seeker_suffix, plasflow_suffix,hmm_suffix):
        self.preprocessing_suffix=preprocessing_suffix
        self.RGI_suffix=RGI_suffix #RGI
        self.DeepARG_suffix=DeepARG_suffix #DeepARG
        self.SARG_suffix=SARG_suffix #ARG_ranking
        self.DVF_suffix=DVF_suffix #DVF
        self.seeker_suffix=seeker_suffix #SEEKER
        self.plasflow_suffix=plasflow_suffix #PLASFLOW
        self.hmm_suffix=hmm_suffix #hmm
        self.workdir=workdir
    
    def check_done():
        print("yes")
        
    def check_file_completness(self,file_basename,inputDir,subDir,prefix): #generate specific name list
        existed=[]
        sample_count=len(file_basename)
            #check AMR
        try:
            if subDir == "AMR": 
                tar=os.path.join(inputDir,prefix,"CompRanking_intermidate",subDir)
                if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".RGI.done")):
                    print("true")
                    if sample_count == (len(os.listdir(tar))/2):
                        print("The AMR prediction has competed with file number of {0}".format(sample_count)) 
                    else:
                        print("Wrong")  
                else:
                    print("AMR prediction is not completed")    
        except:
            raise Exception("Don't existed...") 
        #check MGE
        try:
            if subDir == "MGE":
                DVF_prefix="DVF"
                plasflow_prefix="Plasflow"
                seeker_prefix="Seeker"
                #check DVF
                tar=os.path.join(inputDir,prefix,"CompRanking_intermidate",subDir,DVF_prefix)
                if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".DVF.done")):
                    print("true")
                    try:
                        if sample_count == (len(os.listdir(tar))/2):
                            print("The DVF prediction has competed with file number of {0}".format(sample_count*2))
                    except:
                        raise Exception("DVF prediction is not completed...")
                #check plasflow
                tar=os.path.join(inputDir,prefix,"CompRanking_intermidate",subDir,plasflow_prefix)
                if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".PLASFLOW.done")):
                    try:
                        if sample_count == (len(os.listdir(tar,))/4):
                            print("The PLS prediction has competed with file number of {0}".format(sample_count*4))
                    except:
                        raise Exception("PLS prediction is not completed...")
                #check seeker
                tar=os.path.join(inputDir,prefix,"CompRanking_intermidate",subDir,seeker_prefix)
                if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".SEEKER.done")):    
                    try:
                        if sample_count == (len(os.listdir(tar))/2):
                            print("The Seeker prediction has competed with file number of {0}".format(sample_count*2))
                    except:
                        raise Exception("Seeker prediction is not completed...")
                else:
                    print("MGE prediction is not completed") 
        except:
            raise Exception("MGE is not complete...")
            # #check VF
            # if subDir == "Virulence":
            #     #check VFDB
            #     if os.path.exists(os.path.join(inputDir,prefix,prefix+".VFDB.done")):
            #         tar=os.path.join(inputDir,prefix,subDir)
            #     try:
            #         if sample_count == (len(os.listdir(tar))/2):
            #             print("The VF prediction has competed with file number of {0}".format(sample_count*2))
            #     except:
            #         raise Exception("VF prediction is not completed...")               
        

if __name__ == '__main__':
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    workdir="/lomi_home/gaoyang/software/CompRanking"
    name="CompRanking"
    prefix="CompRanking"
    path_bin="abs_path_to_conda_bin"
    yaml_path=os.path.join(workdir,"test_yaml.yaml")
    conda_path_str="".join(path.read_conda_path(name,path_bin,yaml_path))
    print(type(conda_path_str))
    file_list=path.file_abs_path_list_generation(input_dir)
    base_list=path.file_base_acquire(file_list)
    print(file_list)
    print(base_list)
    
    checkPoints=checkPoints()
    checkPoints.check_done()
    checkPoints.check_file_completness(base_list,input_dir,"MGE",prefix)