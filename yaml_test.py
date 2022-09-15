#!/usr/env/bin python
#encoding: utf-8
import yaml
import os
import optparse
import sys
import glob

parser = optparse.OptionParser()
parser.add_option("-c", "--config_file", action = "store", type = "string", dest = "config_file", 
                  help = "file contains basic configeration information")
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")

(options, args) = parser.parse_args()

config_path=options.config_file
input_dir=options.input_dir

yaml_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),config_path)
print("The yaml file path is:" + yaml_path)


def read_conda_path(name,path,yaml_path):
    CONDA_PATH=[]
    #open yaml file
    with open(yaml_path,"r") as f:
        data=yaml.load(f)
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
    for i in file_abs_path:
        fileBase_1=i.split("/")[-1]
        fileBase_2=fileBase_1.strip(".fa")
        base_list.append(fileBase_2)
    return base_list

#check file completeness
def check_file_completness(file_basename,inputDir,subDir,prefix): #generate specific name list
    existed=[]
    sample_count=len(file_basename)
        #check AMR
    try:
        if subDir == "AMR":
            tar=os.path.join(inputDir,prefix,subDir)
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
            tar=os.path.join(inputDir,prefix,subDir,DVF_prefix)
            if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".DVF.done")):
                print("true")
                try:
                    if sample_count == (len(os.listdir(tar))/2):
                        print("The DVF prediction has competed with file number of {0}".format(sample_count*2))
                except:
                    raise Exception("DVF prediction is not completed...")
            #check plasflow
            tar=os.path.join(inputDir,prefix,subDir,plasflow_prefix)
            if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),prefix+".PLASFLOW.done")):
                try:
                    if sample_count == (len(os.listdir(tar,))/4):
                        print("The PLS prediction has competed with file number of {0}".format(sample_count*4))
                except:
                    raise Exception("PLS prediction is not completed...")
            #check seeker
            tar=os.path.join(inputDir,prefix,subDir,seeker_prefix)
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
        #check VF
        if subDir == "Virulence":
            #check VFDB
            if os.path.exists(os.path.join(inputDir,prefix,prefix+".VFDB.done")):
                tar=os.path.join(inputDir,prefix,subDir)
            try:
                if sample_count == (len(os.listdir(tar))/2):
                    print("The VF prediction has competed with file number of {0}".format(sample_count*2))
            except:
                raise Exception("VF prediction is not completed...")               
      

if __name__ == '__main__':
    name="CompRanking"
    prefix="CompRanking"
    path="abs_path_to_conda_bin"
    conda_path_str="".join(read_conda_path(name,path,yaml_path))
    print(type(conda_path_str))
    file_list=file_abs_path_list_generation(input_dir)
    base_list=file_base_acquire(file_list)
    print(file_list)
    print(base_list)
    
    check_file_completness(base_list,input_dir,"MGE",prefix)
# print(os.path.exists('/lomi_home/gaoyang/miniconda/bin'))


            
        


# fo = open(a,"r")#fileobject
# #使用yaml方法获取数据
# res = yaml.load(fo)
# print(res)


