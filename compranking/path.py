#!/usr/env/bin python
import yaml, glob,os

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
def file_abs_path_list_generation(input_dir):
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

#get AGS file abs path and extract prefix
def getPrefix(input_AGS_dir):
    file_prefix=[]
    file_list=glob.glob(input_AGS_dir+"/*.AGS.txt")
    for i in file_list:
        prefix=((os.path.basename(i)).rstrip("AGS.txt"))
        file_prefix.append(prefix)
    
    return file_list, file_prefix
