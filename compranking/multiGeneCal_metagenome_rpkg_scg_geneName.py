#!/usr/bin/env python
# title             :multiGeneCal_metagenome_rpkg_scg_geneName.py -> Gene_cal.ipynb
# description       :calculate relative abundance of each sub genes
#                    This version is used to calculate rpkg abundance (beta version)
#                    and single copy gene (cell copy)
# author            :Gaoyang Luo
# date              :20240624
# version           :1.0
# usage             :python compranking/multiGeneCal_metagenome_rpkg_scg_geneName.py -i <input_dir>
#                                                                        -p <project_prefix>
#                                                                        -n <normalization_base> #AGS or scg
#                                                                        -t <threads>
#                                                                        -d <pth2KK2db>
#python compranking/multiGeneCal_metagenome_rpkg_scg_geneName.py -i /lomi_home/gaoyang/software/CompRanking/tmp_DSR -n AGS -p DSR
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import numpy as np
import os
import path
import subprocess
import multiprocessing
import optparse
import glob
from Bio import SeqIO

parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest = "input_dir", 
                  help = "director contained input fasta files")
parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "project_prefix",
				 help = "set your project name as global prefix")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_dir",
				 help = "director tontained output files")
parser.add_option("-t", "--threads", action = "store", type = "string", dest = "threads",
				 help = "how many cpus you want use")   
parser.add_option("-n", "--normalization_base", action = "store", type = "string", dest = "normalization_base",
				 help = "scg or average genome length")   
parser.add_option("-c", "--config_file", action = "store", type = "string", dest = "config_file", 
                  help = "file contains basic configeration information, defult: test_yaml.yaml")           
parser.add_option("-d", "--database", action = "store", type = "string", dest = "database",
				  help = "The path to Kranken2 database")
# parser.add_option("-r", "--restart", action = "store", type = "string", dest = "restart",
# 									default='1', help = "restart all the processs")


(options, args) = parser.parse_args()
#path configeration
input_dir=options.input_dir
project_prefix=options.project_prefix
output=options.output_dir
threads=options.threads
normalization_base=options.normalization_base
config_path=options.config_file
database=options.database
output=os.path.join(input_dir,project_prefix,"CompRanking_result")
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.threads is None):
    threads = "24" #default threads
if (options.config_file is None):
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"test_yaml.yaml") #"../test_yaml.yaml" default config_file path
if (options.database is None):
    database = "/lomi_home/gaoyang/db/kraken2/202203"#default config_file path
if (options.normalization_base is None):
    normalization_base = "AGS"
if (options.output_dir is None):
    output = os.path.join(input_dir,project_prefix,"CompRanking_result") #default output directory
if options.normalization_base =="AGS":
    cell_suffix="Cell" #genome equivalents = sample size / AGS
elif options.normalization_base =="16S":
    cell_suffix="16S" 

def get_DB_DeepARG_len(input_deeparg_length):
    #load_Deeparg_structure
    df_deepARG_strucure=pd.read_csv(input_deeparg_length,sep="\t",header=None)
    df_deepARG_strucure
    df_deepARG_strucure["id"]=df_deepARG_strucure[0].str.split("|",expand=True)[0]
    #write dic
    DB_deepARG_length={}
    for i, name in df_deepARG_strucure.iterrows():
        DB_deepARG_length.setdefault(str(name["id"]), name[1])
    return DB_deepARG_length
    
    
def get_DB_SARG_len(input_sarg_structure):
    #load SARG_structure
    df_sarg_structure=pd.read_csv(input_sarg_structure, sep="\t", header=None)
    df_sarg_structure=df_sarg_structure[[1,2]]
    #write dic
    DB_SARG_length={}
    for i, name in df_sarg_structure.iterrows():
        DB_SARG_length.setdefault(str(name[1]), name[2])
    return DB_SARG_length     

def get_MobilOG_len(input_mobileOG_structure):
    #load mobileOG_structure
    df_mobileOG_structure=pd.read_csv(input_mobileOG_structure,sep="\t", header=0)
    #write dic
    DB_MobileOG_length={}
    for i, name in df_mobileOG_structure.iterrows():
        DB_MobileOG_length.setdefault(str(name["mobileOG_ID"]), name["length"])
    return DB_MobileOG_length
    
def getPrefix(input_AGS_dir):
    file_prefix=[]
    file_list=glob.glob(input_AGS_dir+"/*.AGS.txt")
    for i in file_list:
        prefix=((os.path.basename(i)).rstrip("AGS.txt"))
        file_prefix.append(prefix)
    
    return file_list, file_prefix

def get_genome_len(input_AGS, prefix_list):
    genome_length_dic={}
    for index, j in enumerate(input_AGS): 
        for lines in open(j,'r'):
                            if lines.startswith('genome_equivalents'):
                                lines_set = lines.split('\n')[0].split('\t')
                                genome_length = float(lines_set[1])
                                genome_length_dic.setdefault(prefix_list[index],float(genome_length))
                                print("The Average Genome Length of file {} is {}".format(prefix_list[index],genome_length ))
    return genome_length_dic  

def RB_gene_sum(DB_deepARG_length,DB_SARG_length, DB_MobileOG_length, 
                input_AMR_sum,input_kk2,input_deeparg_sure,
                input_rgi,input_SARG,input_scg,input_rpkm,input_indexFile,genome_length,filebase):
    if normalization_base =="AGS":
        cell_suffix="Cell"
    elif normalization_base =="16S":
        cell_suffix="16S"   
    #load final output
    df_AMR_sum=pd.read_csv(input_AMR_sum,sep="\t",header=0)
    df_AMR_hit=df_AMR_sum[df_AMR_sum.ARG_prediction != "-"]
    df_AMR_hit1=df_AMR_hit[["ORF_ID","ARG_prediction","ARG_class","Database","CompRanking_MGE_prediction"]]
    df_AMR_hit1["db_final"]=df_AMR_hit1["Database"].str.split("/", expand=True)[0]
    
    #record hit database and orf_id
    Record_db_orf={}
    Record_ARG_name_orf={}
    for i, name in df_AMR_hit1.iterrows():
        Record_db_orf.setdefault(str(name["ORF_ID"]), str(name["db_final"]))
        Record_ARG_name_orf.setdefault(str(name["ORF_ID"]), [str(name["ARG_prediction"]), str(name["ARG_class"]),str(name["CompRanking_MGE_prediction"])])
    #Load normalization base 16s or AGS
    #load kk2
    kraken=input_kk2
    scg=input_scg
    copy_16S = 1
    if normalization_base == "AGS":
        gene_length = genome_length
    else:
        gene_length = 1550

    #metagenomes_kk2_16s_bases
    for lines in open(kraken,'r'):
        content = lines.split('\n')[0].split('\t')
        if 'Bacteria' in lines:
            copy_16S = float(content[1]) # pair end
            break
    
    #scg_bases
    ##scg filter
    df_scg=pd.read_csv(scg, sep="\t",header=None)
    df_scg.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']

    ##filter out contigs identity under 30
    df_scg_iden50 = df_scg[df_scg.identity > 30]
    ##filter out contigs length larger than 0
    df_scg_len = df_scg_iden50[df_scg_iden50.alignLen > 0]
    ##filter out bit larger than 50
    df_scg_bit50 = df_scg_len[df_scg_len.bit > 50]
    
    num_scg=len(df_scg_bit50)
    
    #load rpkm alinged reads number and mapped reads number
    # input_rpkm="/lomi_home/gaoyang/software/CompRanking/tmp_DSR/DSR/CompRanking_intermediate/preprocessing/5M_contigs/cov/S0PCL_clean.sorted_filtered.rpkm"
    #get reads numbers
    for lines in open(input_rpkm,'r'):
        #content = lines.split('\n')[0].split('\t')
        content = lines.split('\t')
        if '#Mapped' in lines:
            mapped_reads = float(content[1]) # pair end
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
    for i in array:
        rpkm_dic.setdefault(str(i[0]),float(i[4]))
    #change orf id to contigs id
    df_index=pd.read_csv(input_indexFile,sep="\t",header=None)
    array_indexFile=np.array(df_index)
    array_indexFile=array_indexFile.tolist()
    index_dic={}
    for i in array_indexFile:
        index_dic.setdefault(str(i[1]),str(i[0]))
        
    
    #load deeparg
    df_deeparg_sure=pd.read_csv(input_deeparg_sure, sep="\t")
    #load rgi
    df_RGI=pd.read_csv(input_rgi, sep="\t")
    df_RGI=df_RGI.fillna("-")
    #load sarg
    df_SARG=pd.read_csv(input_SARG, sep="\t", header=0,index_col=0)
    df_SARG["class"]=df_SARG["class"].str.upper()
    df_SARG.columns=["ORF","query","class","Phenotype","ARG_rank"]

    #get each ARG output dic
    deepARG_output_dic={}
    deepARG_output_name_dic={}
    SARG_output_dic={}
    SARG_output_name_dic={}
    for i, name in df_deeparg_sure.iterrows():
        deepARG_output_dic.setdefault(str(name["read_id"]), str(name["best-hit"].split("|")[0]))
        deepARG_output_name_dic.setdefault(str(name["read_id"]), str(name["#ARG"]))
    for i, name in df_SARG.iterrows():
        SARG_output_dic.setdefault(str(name["ORF"]), str(name["query"]))
        SARG_output_name_dic.setdefault(str(name["ORF"]), str(name["Phenotype"]))

    #generate DB_xxx_length_res
    DB_deepARG_length_res={}
    DB_SARG_length_res={}
    DB_CARD_length_res={}
    RGI_output_name_dic={}
    #deeparg_output_ARG_cal
    for i in deepARG_output_dic:
        DB_deepARG_length_res.setdefault(i,DB_deepARG_length[deepARG_output_dic[i]])
    #SARG_output_ARG_cal
    for i in SARG_output_dic:
        DB_SARG_length_res.setdefault(i,DB_SARG_length[SARG_output_dic[i]])
    #RGI_output_ARG_cal
    for i, name in df_RGI.iterrows():
        DB_CARD_length_res.setdefault(str(name["ORF_ID"]), len(str(name["CARD_Protein_Sequence"])))
        RGI_output_name_dic.setdefault(str(name["ORF_ID"]),str(name["Best_Hit_ARO"]))
    
    #cal ARGs relative abundance 16S
    #cal ARGs relative abundance RPKM
    #RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
    """
    Normalization for comparing gene coverage values. 
    RPKM corrects differences in both: sample sequencing depth and gene length.
    RPKM is introduced in
    http://www.ncbi.nlm.nih.gov/pubmed/18516045
    
    Here, we use kilo base, RPKK, defined as follows:
    
    RPKG - reads per kilo base per average genome size
    RPKG = numReads / ( geneLength / 1,000 * totalNumReads / genome equivalents )
    cell copy = (mapped_reads / DB_deepARG_length_res[orf]) / num_scg

        numReads                - number of reads mapped to a gene sequence
        geneLength              - length of the gene sequence
        totalNumReads           - total number of mapped reads of a sample
        genome equivalents      - sample bases / average genome size
        num_scg                 - abundance of single copy gene
    """
    abundance_arg_16S=0
    abundance_arg_RPKM=0
    RPKM_ARG={}
    TAXO_ARG={}
    RPKG_ARG_NAME={}
    num_contigs=len(df_AMR_sum) #totalNumReads
    num_mapped_reads= mapped_reads
    for orf in Record_db_orf:
        find_db=''
        contig_orf=index_dic[orf]
        print("The orf {0} corresponding contig number is {1}: ".format(
            orf, index_dic[orf]))
        if Record_db_orf[orf]:
            find_db=Record_db_orf[orf]
            ARG_name=Record_ARG_name_orf[orf][0].split("/")[0]
            ARG_class=Record_ARG_name_orf[orf][1].split("/")[0]
            MGE_type=Record_ARG_name_orf[orf][2]
            if ARG_name:
                print(ARG_name)
            else:
                print("No ARG name")
            if ARG_class:
                print(ARG_class)
            else:
                print("No ARG class")
            
            if contig_orf in rpkm_dic:
                mapped_reads=rpkm_dic[contig_orf]
            else:
                mapped_reads=1
            if find_db=="DeepARG":
                abundance_arg_16S += (mapped_reads / DB_deepARG_length_res[orf]) / num_scg #scg
                abundance_arg_RPKM += mapped_reads / (DB_deepARG_length_res[orf] / 1000 * gene_length)#rpkg
                TAXO_ARG.setdefault(str(orf), float((mapped_reads / DB_deepARG_length_res[orf]) / num_scg)) #scg
                RPKM_ARG.setdefault(str(orf), float(mapped_reads / (DB_deepARG_length_res[orf] / 1000 * gene_length))) #rpkg
                RPKG_ARG_NAME.setdefault(str(orf), [str(ARG_name), str(ARG_class), float(mapped_reads / (DB_deepARG_length_res[orf] / 1000 * gene_length)), str(find_db), str(MGE_type)]) #rpkg
            elif find_db=="RGI":
                abundance_arg_16S += (mapped_reads / DB_CARD_length_res[orf])/ num_scg #scg
                abundance_arg_RPKM  += mapped_reads / (DB_CARD_length_res[orf] / 1000 * gene_length)#rpkg
                TAXO_ARG.setdefault(str(orf), float((mapped_reads / DB_CARD_length_res[orf])/ num_scg)) #scg
                RPKM_ARG.setdefault(str(orf), float(mapped_reads / (DB_CARD_length_res[orf] / 1000 * gene_length)))#rpkg
                RPKG_ARG_NAME.setdefault(str(orf), [str(ARG_name), str(ARG_class), float(mapped_reads / (DB_CARD_length_res[orf] / 1000 * gene_length)), str(find_db), str(MGE_type)]) #rpkg
            elif find_db=="SARG":
                abundance_arg_16S += (mapped_reads / DB_SARG_length_res[orf]) / num_scg #scg
                abundance_arg_RPKM  += mapped_reads / (DB_SARG_length_res[orf] / 1000 * gene_length)#rpkg
                TAXO_ARG.setdefault(str(orf), float((mapped_reads / DB_SARG_length_res[orf]) / num_scg)) #scg
                RPKM_ARG.setdefault(str(orf), float(mapped_reads / (DB_SARG_length_res[orf] / 1000 * gene_length)))#rpkg
                RPKG_ARG_NAME.setdefault(str(orf), [str(ARG_name), str(ARG_class), float(mapped_reads / (DB_SARG_length_res[orf] / 1000 * gene_length)), str(find_db), str(MGE_type)]) #rpkg
            else:
                continue
    # print(abundance_arg_16S, abundance_arg_RPKM)   
    print(f"The relative abundance of ARG by {cell_suffix}(AGS) is: {abundance_arg_16S}")
    print("The relative abundance of ARG by RPKG is: {}".format(abundance_arg_RPKM))
    
    #Transfer RPKG_ARG_NAME to {ARG: [class, abundance]}
    # Initialize new dictionary
    RPKG_ARG_NAME_abundance = {}

    # Traverse the original dictionary
    for orf, values in RPKG_ARG_NAME.items():
        arg_name = values[0]
        arg_class = values[1]
        value = values[2]
        find_db = values[3]
        mge_type = values[4]
        
        if arg_name in RPKG_ARG_NAME_abundance:
            # If ARG already exists, accumulate the value
            RPKG_ARG_NAME_abundance[arg_name][1] += value
        else:
            # If ARG does not exist, create a new entry
            RPKG_ARG_NAME_abundance[arg_name] = [arg_class, value, find_db, mge_type]
    
    #cal ARG subtype
    """
    subtypes including:
        All, multidrug, beta-lactam, aminoglycoside,
        tetracycline, sulfonamide, MLS, bacitracin,
        chloramphenicol, quinlone, fosmidomycin, trimethoprim,
        kasugamycin, vancomycin, rifamycin, fosfomycin, belomycin, unclassified...
    """
    abundance_ARG_subtype_16S={}
    abundance_ARG_subtype_RPKM={}
    for i, name in df_AMR_hit.iterrows():
        tmp_16s=0
        tmp_rpkm=0
        #16s
        if abundance_ARG_subtype_16S.get(name["ARG_class"].split("/")[0].split(":")[0].strip(";")):
            tmp_16s=abundance_ARG_subtype_16S.get(name["ARG_class"].split("/")[0].split(":")[0].strip(";")) + TAXO_ARG[name["ORF_ID"]]
            abundance_ARG_subtype_16S[name["ARG_class"].split("/")[0].split(":")[0].strip(";")] = tmp_16s
        else:
            abundance_ARG_subtype_16S.setdefault(str(name["ARG_class"].split("/")[0].split(":")[0].strip(";")), float(TAXO_ARG[name["ORF_ID"]]))
        #rpkm
        if abundance_ARG_subtype_RPKM.get(name["ARG_class"].split("/")[0].split(":")[0].strip(";")):
            tmp_rpkm=abundance_ARG_subtype_RPKM.get(name["ARG_class"].split("/")[0].split(":")[0].strip(";")) + RPKM_ARG[name["ORF_ID"]]
            abundance_ARG_subtype_RPKM[name["ARG_class"].split("/")[0].split(":")[0].strip(";")] = tmp_rpkm    
        else:
            abundance_ARG_subtype_RPKM.setdefault(str(name["ARG_class"].split("/")[0].split(":")[0].strip(";")), float(RPKM_ARG[name["ORF_ID"]]))   
    
    
    
    
    ###################### MGE relative abundance calculation####################
    #get DB_mobile_OG_len_dic
    df_mobileOG_structure=pd.read_csv(input_mobileOG_structure,sep="\t", header=0)
    #get MGE reference length
    DB_MobileOG_length={}
    for i, name in df_mobileOG_structure.iterrows():
        DB_MobileOG_length.setdefault(str(name["mobileOG_ID"]), name["length"])
    #load MGE result
    df_MGE_hit=df_AMR_sum[df_AMR_sum.mobileOG_ID != "-"]
    #generate DB_MobileOG_length_res
    DB_MobileOG_length_res={}
    for i, name in df_MGE_hit.iterrows():
        DB_MobileOG_length_res.setdefault(str(name["ORF_ID"]),DB_MobileOG_length[name["mobileOG_ID"]])
    
    #cal MGEs relative abundance 
    abundance_MGE_16S=0
    abundance_MGE_RPKM=0
    RPKM_MGE={}
    TAXO_MGE={}
    RPKG_MGE_NAME={}
    for orf_MGE in DB_MobileOG_length_res:
        contig_orf=index_dic[orf_MGE]
        print("The orf {0} corresponding contig number is {1}: ".format(
            orf_MGE, index_dic[orf_MGE]))
        mapped_reads=rpkm_dic[contig_orf]
        if contig_orf in rpkm_dic:
            mapped_reads=rpkm_dic[contig_orf]
        else:
            mapped_reads=1
        abundance_MGE_16S += (mapped_reads / DB_MobileOG_length_res[orf_MGE]) / num_scg
        abundance_MGE_RPKM += mapped_reads / (DB_MobileOG_length_res[orf_MGE] / 1000 * gene_length)#rpkg
        TAXO_MGE.setdefault(str(orf_MGE),float((mapped_reads / DB_MobileOG_length_res[orf_MGE]) / num_scg))
        RPKM_MGE.setdefault(str(orf_MGE),float(mapped_reads / (DB_MobileOG_length_res[orf_MGE] / 1000 * gene_length)))#rpkg
    
    print(f"The relative abundance of MGE by {cell_suffix} is: {abundance_MGE_16S}")
    print("The relative abundance of MGE by RPKM is: {}".format(abundance_MGE_RPKM))
    print(TAXO_MGE)
    print(RPKM_MGE)
    #cal MGE subtype
    Insertion_Sequences_db=["ISFinder"]
    Integrative_Elements_db=["AICE","ICE","CIME","IME","immedb"]
    Plasmids_db=["COMPASS","Plasmid RefSeq"]
    Bacteriophages_db=["pVOG","GPD"]
    Multiple_db=["ACLAME", "Multiple"]
    abundance_MGE_subtype_16S={}
    abundance_MGE_subtype_RPKM={}
    for i,name in df_MGE_hit.iterrows():
        tmp_16s=0
        tmp_rpkm=0
        #record phage into dic
        if name["MGE_Database"] in Bacteriophages_db:
            #16s
            if abundance_MGE_subtype_16S.get("phage"):
                tmp_16s = abundance_MGE_subtype_16S.get("phage") + TAXO_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_16S["phage"] = tmp_16s                                 
            else:
                abundance_MGE_subtype_16S.setdefault(str("phage"), float(TAXO_MGE[name["ORF_ID"]]))
            #rpkm   
            if abundance_MGE_subtype_RPKM.get("phage"):
                tmp_rpkm = abundance_MGE_subtype_RPKM.get("phage") + RPKM_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_RPKM["phage"] = tmp_rpkm  
            else:
                abundance_MGE_subtype_RPKM.setdefault(str("phage"), float(RPKM_MGE[name["ORF_ID"]]))
        
        #record plasmid into dic
        elif name["MGE_Database"] in Plasmids_db:
            #16s
            if abundance_MGE_subtype_16S.get("plasmid"):
                tmp_16s = abundance_MGE_subtype_16S.get("plasmid") + TAXO_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_16S["plasmid"] = tmp_16s
            else:
                abundance_MGE_subtype_16S.setdefault(str("plasmid"), float(TAXO_MGE[name["ORF_ID"]]))
            #rpkm   
            if abundance_MGE_subtype_RPKM.get("plasmid"):
                tmp_rpkm = abundance_MGE_subtype_RPKM.get("plasmid") + RPKM_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_RPKM["plasmid"] = tmp_rpkm  
            else:
                abundance_MGE_subtype_RPKM.setdefault(str("plasmid"), float(RPKM_MGE[name["ORF_ID"]]))
        
        #record Insertion_Sequences into dic
        elif name["MGE_Database"] in Insertion_Sequences_db:
            #16s
            if abundance_MGE_subtype_16S.get("Insertion_Sequences"):
                tmp_16s = abundance_MGE_subtype_16S.get("Insertion_Sequences") + TAXO_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_16S["Insertion_Sequences"] = tmp_16s
            else:
                abundance_MGE_subtype_16S.setdefault(str("Insertion_Sequences"), float(TAXO_MGE[name["ORF_ID"]]))
            #rpkm   
            if abundance_MGE_subtype_RPKM.get("Insertion_Sequences"):
                tmp_rpkm = abundance_MGE_subtype_RPKM.get("Insertion_Sequences") + RPKM_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_RPKM["Insertion_Sequences"] = tmp_rpkm  
            else:
                abundance_MGE_subtype_RPKM.setdefault(str("Insertion_Sequences"), float(RPKM_MGE[name["ORF_ID"]]))
        
        #record Integrative_Elements into dic
        elif name["MGE_Database"] in Integrative_Elements_db:
            #16s
            if abundance_MGE_subtype_16S.get("Integrative_Elements"):
                tmp_16s = abundance_MGE_subtype_16S.get("Integrative_Elements") + TAXO_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_16S["Integrative_Elements"] = tmp_16s
            else:
                #16s
                abundance_MGE_subtype_16S.setdefault(str("Integrative_Elements"), float(TAXO_MGE[name["ORF_ID"]]))
            #rpkm   
            if abundance_MGE_subtype_RPKM.get("Integrative_Elements"):
                tmp_rpkm = abundance_MGE_subtype_RPKM.get("Integrative_Elements") + RPKM_MGE[name["ORF_ID"]]
                abundance_MGE_subtype_RPKM["Integrative_Elements"] = tmp_rpkm
            else:
                abundance_MGE_subtype_RPKM.setdefault(str("Integrative_Elements"), float(RPKM_MGE[name["ORF_ID"]]))
                
        #record phage and plasmid from ACLAME into dic
        elif name["MGE_Database"] in Multiple_db:
            if name["Taxonomy"] == "phage":
                #16s
                if abundance_MGE_subtype_16S.get("phage"):
                    tmp_16s = abundance_MGE_subtype_16S.get("phage") + TAXO_MGE[name["ORF_ID"]]
                    abundance_MGE_subtype_16S["phage"] = tmp_16s
                else:
                    abundance_MGE_subtype_16S.setdefault(str("phage"), float(TAXO_MGE[name["ORF_ID"]]))
                #rpkm   
                if abundance_MGE_subtype_RPKM.get("phage"):
                    tmp_rpkm = abundance_MGE_subtype_RPKM.get("phage") + RPKM_MGE[name["ORF_ID"]]
                    abundance_MGE_subtype_RPKM["phage"] = tmp_rpkm  
                else:
                    abundance_MGE_subtype_RPKM.setdefault(str("phage"), float(RPKM_MGE[name["ORF_ID"]]))
            else: #name["Taxonomy"] != "phage"
                #16s
                if abundance_MGE_subtype_16S.get("plasmid"):
                    tmp_16s = abundance_MGE_subtype_16S.get("plasmid") + TAXO_MGE[name["ORF_ID"]]
                    abundance_MGE_subtype_16S["plasmid"] = tmp_16s
                else:
                    abundance_MGE_subtype_16S.setdefault(str("plasmid"), float(TAXO_MGE[name["ORF_ID"]]))    
                #rpkm  
                if abundance_MGE_subtype_RPKM.get("plasmid"):
                    tmp_rpkm = abundance_MGE_subtype_RPKM.get("plasmid") + RPKM_MGE[name["ORF_ID"]]
                    abundance_MGE_subtype_RPKM["plasmid"] = tmp_rpkm  
                else:
                    abundance_MGE_subtype_RPKM.setdefault(str("plasmid"), float(RPKM_MGE[name["ORF_ID"]]))
    print(abundance_MGE_subtype_16S, abundance_MGE_subtype_RPKM)
    
    ###################combine it using a list##########################
    result=[abundance_arg_16S,abundance_arg_RPKM, abundance_MGE_16S, abundance_MGE_RPKM]
    df_ARG_subtype_16S = pd.DataFrame(pd.Series(abundance_ARG_subtype_16S))
    df_ARG_subtype_RPKM = pd.DataFrame(pd.Series(abundance_ARG_subtype_RPKM))
    df_MGE_subtype_16S = pd.DataFrame(pd.Series(abundance_MGE_subtype_16S))
    df_MGE_subtype_RPKM = pd.DataFrame(pd.Series(abundance_MGE_subtype_RPKM))
    
    return result, df_ARG_subtype_16S,df_ARG_subtype_RPKM, df_MGE_subtype_16S, df_MGE_subtype_RPKM, RPKG_ARG_NAME, RPKG_ARG_NAME_abundance

def kk2(file_name_base):
    #run kranken2
    i=file_name_base
    if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle.txt"):
        print("It seems that we have already done the {} KK2 taxonomy annotation...".format(i))
        pass
    else:
        print("KK2 mpaStyle output don't exist... {}".\
            format(os.path.join(input_dir, project_prefix, "CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle.txt"))
        subprocess.call(["bash", kk2_script, 
            "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str, "-d", database, "-n", i])
        
def cov_rpkm(file_name_base):
    #run cov_rpkm_calculation
    i=file_name_base
    if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs/cov")+"/"+i+".rpkm"):
        print("It seems that we have already done the {} rpkm files...".format(i))
        pass
    if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs/cov")+"/"+i+".cov"):
        print("It seems that we have already done the {} cov files...".format(i))
        pass
    else:
        print("KK2 mpaStyle output don't exist... {}".\
            format(os.path.join(input_dir, project_prefix, "CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle_16S.txt"))
        subprocess.call(["bash", cov_rpkm_script, 
            "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str])


def Calculation(file_name_base):
    #calculate relative abundance of functional genes
    i=file_name_base
    if normalization_base =="AGS":
        cell_suffix="Cell"
    elif normalization_base =="16S":
        cell_suffix="16S" 
    try:
        #load ARGs result and relative files
        input_rgi=os.path.join(input_dir,project_prefix,
                            "CompRanking_intermediate/AMR/RGI",
                                    i+"_5M_contigs.RGI.out.txt")
        input_SARG=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/AMR/ARGranking",
                                    i+"_SARGrank_Protein60_Result.tsv")
        input_deeparg_sure=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/AMR/DeepARG", 
                                    i+"_5M_contigs_DeepARG.out.mapping.ARG")
        if normalization_base == "AGS":
            input_kk2=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_report_kk2_mpaStyle.txt")
        elif normalization_base == "16S":
            input_kk2=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_report_kk2_mpaStyle_16S.txt")
        else:
            raise ImportError("No potinted normalization type")
        input_AMR_sum=os.path.join(input_dir,project_prefix,
                                "CompRanking_result",
                                    "CompRanking_"+i+"_AMR_MOB_prediction.tsv")
        input_scg=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/cov",
                                    i+"_5M_contigs_scg_Protein_dimond.txt")
        input_rpkm=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/cov",
                                    i+"_5M_contigs_gene.rpkm")
        input_indexFile=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_5M_contigs.fna2faa.index")
    except:
        raise ValueError("Missing the output...")
    
    try:
        #load reference length
        DB_deepARG_length = get_DB_DeepARG_len(input_deeparg_length)
        DB_SARG_length =get_DB_SARG_len(input_sarg_structure)
        DB_MobileOG_length=get_MobilOG_len(input_mobileOG_structure)
    except:
        raise SystemError("Can't load reference length, please check the original files...")
    
    #Find average genome length
    try:
        #load AGS
        input_AGS_dir=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/AGS")
        file_list,prefix=getPrefix(input_AGS_dir)
        genome_length_dic=get_genome_len(file_list,prefix)
    except:
        raise ValueError("Can't load AGS file, please check the original files,\
                         or check the name of your input file and see whether meet our name rules...")
    
    try:
        #load subtype dataframe
        result, \
        df_ARG_subtype_16S, \
        df_ARG_subtype_RPKM, \
        df_MGE_subtype_16S, \
        df_MGE_subtype_RPKM, \
        RPKG_ARG_NAME, \
        RPKG_ARG_NAME_abundance = RB_gene_sum(DB_deepARG_length,
                                                DB_SARG_length, 
                                                DB_MobileOG_length, 
                                                input_AMR_sum,
                                                input_kk2,
                                                input_deeparg_sure,
                                                input_rgi,
                                                input_SARG,
                                                input_scg,
                                                input_rpkm,
                                                input_indexFile,
                                                genome_length_dic[i],
                                                i)
        #output total relative abundance in a list    
        output_abundance="\t".join(map(str, result))
        #save output as tmp file
        df_ARG_subtype_16S.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            i+"_ARG_"+cell_suffix+"Abu_tmp.txt"),
                                sep="\t",header=False)
        df_ARG_subtype_RPKM.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            i+"_ARG_rpkmAbu_tmp.txt"),
                                sep="\t",header=False)
        df_MGE_subtype_16S.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            i+"_MGE_"+cell_suffix+"Abu_tmp.txt"),
                                sep="\t",header=False)
        df_MGE_subtype_RPKM.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            i+"_MGE_rpkmAbu_tmp.txt"),
                                sep="\t",header=False)
        
        with open(os.path.join(input_dir,project_prefix,"CompRanking_result/Gene_Abundance_Sum_scg-rpkg.txt"), "a") as f:
            f.write("\n" + i + "\t" + output_abundance)
        
        # Convert RPKG_ARG_NAME to {orf: [ARG, class, abundance]}
        RPKG_ARG_NAME_tsv_data = "orf\tARG_name\tARG_class\tvalue\n"  # Add header
        for orf, values in RPKG_ARG_NAME.items():
            RPKG_ARG_NAME_tsv_data += f"{orf}\t{values[0]}\t{values[1]}\t{values[2]}\n" # 0ARG 1class 2abundance
        # Optionally, write to a file
        with open(os.path.join(input_dir,project_prefix,"CompRanking_result","Abundance_orf_gene",i+"_Gene_Abundance_ORF_geneName_class_Cell(GE).txt"), "w") as f:
            f.write(RPKG_ARG_NAME_tsv_data)
        
        # Convert RPKG_ARG_NAME_abundance to {ARG: [class, abundance]}
        RPKG_ARG_NAME_abundance_tsv_data = "ARG_name\tClass\tfind_db\tMGE_type\tValue\n" # Add header
        # Traverse the new dictionary and append data to TSV string
        for arg, values in RPKG_ARG_NAME_abundance.items():
            RPKG_ARG_NAME_abundance_tsv_data += f"{arg}\t{values[0]}\t{values[2]}\t{values[3]}\t{values[1]}\n" # 0class 1abundance
        # Optionally, write the TSV data to a file
        with open(os.path.join(input_dir,project_prefix,"CompRanking_result",i+"_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt"), "w") as f:
            f.write(RPKG_ARG_NAME_abundance_tsv_data)
        
    except:
        raise ValueError("Write to summary abundacne file failed...")
    
    #check abu tmp files
    check_point_list=[]
    try:
        if os.path.exists(os.path.join(input_dir,project_prefix,"CompRanking_result",i +"_ARG_"+cell_suffix+"Abu_tmp.txt")):
            pass
        else:
            print("ARG subtype abundance scg cal file doesn't exit...")
            exit(1)
        if os.path.exists(os.path.join(input_dir,project_prefix,"CompRanking_result",i +"_ARG_rpkmAbu_tmp.txt")):
            pass
        else:
            print("ARG subtype abundance rpkg cal file doesn't exit...")
            exit(1)
        if os.path.exists(os.path.join(input_dir,project_prefix,"CompRanking_result",i +"_MGE_"+cell_suffix+"Abu_tmp.txt")):
            pass
        else:
            print("MGE subtype abundance scg cal file doesn't exit...")
            exit(1)
        if os.path.exists(os.path.join(input_dir,project_prefix,"CompRanking_result",i +"_MGE_rpkmAbu_tmp.txt")):
            pass
        else:
            print("MGE subtype abundance rpkg cal file doesn't exit...")
            exit(1)
        if os.path.exists(os.path.join(input_dir,project_prefix,"CompRanking_result",i+"_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt")):
            pass
        else:
            print("geneName file doesn't exit...")
            exit(1)
    except:
        raise FileNotFoundError("subtype tmp cal file miss...")

def multiKK2():    
    openthreads = len(file_name_base) 
    exfiles = []
    for i in range(openthreads):
        worker = multiprocessing.Process(target=kk2,args=([file_name_base[i]]))
        worker.start()
        print("Now processing:{}".format(file_name_base[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()  


def multiCalculation():    
    openthreads = len(file_name_base) 
    exfiles = []
    for i in range(openthreads):
        worker = multiprocessing.Process(target=Calculation,args=([file_name_base[i]]))
        worker.start()
        print("Now processing:{}".format(file_name_base[i]))
        exfiles.append(worker)

    for worker in exfiles:
        worker.join()  

if __name__ == "__main__":
    
    #global settings
    #gloab settings
    config_path=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"test_yaml.yaml")
    # input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    # output=os.path.join(input_dir,"CompRanking/CompRanking_result")
    # project_prefix="CompRanking"
    # database="/lomi_home/gaoyang/db/kraken2/202203"
    # threads="24"
    kk2_script=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"scripts/kk2_run_single_16S.sh")
    cov_rpkm_script=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"scripts/cov_rpkm_calculation.sh")
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    yaml_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),config_path)
    #Write in abs conda path
    path_bin="abs_path_to_conda_bin"
    conda_path_str="".join(path.read_conda_path("CompRanking",path_bin,yaml_path)) #record abs path of conda bin
    print("The absolute path to conda bin is:{0}".format(conda_path_str)) 
    
    
    #input_db
    # input_deeparg_length="/lomi_home/gaoyang/software/CompRanking/databases/deepargdata1.0.2/database/v2/features.gene.length"
    # input_sarg_structure="/lomi_home/gaoyang/software/CompRanking/databases/SARG/SARG.db.fasta.length"
    # input_mobileOG_structure="/lomi_home/gaoyang/software/CompRanking/databases/MobileOG-db/MobileOG-db_structure.tsv"
    input_deeparg_length=os.path.join(
                            os.path.dirname(
                                os.path.dirname(
                                    os.path.abspath(__file__))),
                                "databases/deepargdata1.0.2/database/v2/features.gene.length")
    input_sarg_structure=os.path.join(
                            os.path.dirname(
                                os.path.dirname(
                                    os.path.abspath(__file__))),
                                "databases/SARG/SARG.db.fasta.length")
    input_mobileOG_structure=os.path.join(
                                os.path.dirname(
                                    os.path.dirname(
                                        os.path.abspath(__file__))),
                                "databases/MobileOG-db/MobileOG-db_structure.tsv")
    
    #input_result
    # input_AMR_sum="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_result/CompRanking_ERR1191817.contigs_AMR_prediction.tsv"
    # input_kk2="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_report_kk2_mpaStyle.txt"
    # input_deeparg_sure="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/DeepARG/ERR1191817.contigs_5M_contigs_DeepARG.out.mapping.ARG"
    # input_rgi="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/RGI/ERR1191817.contigs_5M_contigs.RGI.out.txt"
    # input_SARG="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/ARGranking/ERR1191817.contigs_SARGrank_Protein60_Result.tsv"
    
    #run kranken2
    if normalization_base == "AGS":
        for i in file_name_base:
            if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle.txt"):
                print("It seems that we have already done the {} KK2 taxonomy annotation...".format(i))
                continue
            else:
                print("KK2 mpaStyle output don't exist... {}".\
                    format(os.path.join(input_dir, project_prefix, "CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle.txt"))
                subprocess.call(["bash", kk2_script, 
                    "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str, "-d", database, "-n", i])
    else:
        for i in file_name_base:
            if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle_16S.txt"):
                print("It seems that we have already done the {} KK2 taxonomy annotation...".format(i))
                continue
            else:
                print("KK2 mpaStyle output don't exist... {}".\
                    format(os.path.join(input_dir, project_prefix, "CompRanking_intermediate/preprocessing/5M_contigs")+"/"+i+"_report_kk2_mpaStyle_16S.txt"))
                subprocess.call(["bash", kk2_script, 
                    "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str, "-d", database, "-n", i])
        
    #run cov_rpkm
    for i in file_name_base:
        if os.path.exists(os.path.join(input_dir, project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs/cov")+"/"+i+"_5M_contigs_gene.rpkm"):
            print("It seems that we have already done the {} rpkm files...".format(i))
            continue
        else:
            print("reads alignment file don't exist... {}".\
                format(os.path.join(input_dir, project_prefix, "CompRanking_intermediate/preprocessing/5M_contigs/cov")+"/"+i+"_5M_contigs_gene.rpkm"))
            subprocess.call(["bash", cov_rpkm_script, 
                "-i", input_dir, "-t", threads, "-p", project_prefix, "-m", conda_path_str])
    #cov_rpkm_calculation.sh
    #multiCalculation
    # multiKK2()
    #multiTask
    multiCalculation() 
    
    # #########concat all the abu result############
    #concat ARG result
    #concat scg
    if normalization_base =="AGS":
        cell_suffix="Cell" #scg+rpkg
    elif normalization_base =="16S":
        cell_suffix="16S"
        
    name_list_16S=[]
    for i in file_name_base:
        name_list_16S.append(i+"_ARG_"+cell_suffix+"Abu_tmp.txt")
    init=0
    df_main=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_16S[0]),sep="\t", header=None)
    df_main.columns=["type",name_list_16S[0]]
    for i,name in enumerate(name_list_16S):
        if i < len(name_list_16S)-1:
            init+=1
            if name_list_16S[init]:
                df_2=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_16S[init]),sep="\t", header=None)
                df_2.columns=["type",name_list_16S[init]]
                df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
    #save 16s subtype abu
    df_main.to_csv(os.path.join(
                input_dir, project_prefix,
                    "CompRanking_result",
                        project_prefix+"_Abundance_ARGs_subtypes_"+cell_suffix+".txt"),sep="\t",index=None)
    #cal rpkm
    name_list_rpkm=[]
    for i in file_name_base:
        name_list_rpkm.append(i+"_ARG_rpkmAbu_tmp.txt")
    init=0
    df_main=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_rpkm[0]),
                        sep="\t", header=None)
    df_main.columns=["type",name_list_rpkm[0]]
    for i,name in enumerate(name_list_rpkm):
        if i < len(name_list_16S)-1:
            init+=1
            if name_list_rpkm[init]:
                df_2=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_rpkm[init]),
                                 sep="\t", header=None)
                df_2.columns=["type",name_list_rpkm[init]]
                df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
    #save rpkm subtype abu
    df_main.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            project_prefix+"_Abundance_ARGs_subtypes_rpkg.txt"),sep="\t",index=None)
    #concat MGE result
    #concat 16S
    name_list_16S=[]
    for i in file_name_base:
        name_list_16S.append(i+"_MGE_"+cell_suffix+"Abu_tmp.txt")
    init=0
    df_main=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_16S[0]),sep="\t", header=None)
    df_main.columns=["type",name_list_16S[0]]
    for i,name in enumerate(name_list_16S):
        if i < len(name_list_16S)-1:
            init+=1
            if name_list_16S[init]:
                df_2=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_16S[init]),sep="\t", header=None)
                df_2.columns=["type",name_list_16S[init]]
                df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
    #save 16s subtype abu
    df_main.to_csv(os.path.join(
                input_dir,project_prefix,
                    "CompRanking_result",
                        project_prefix+"_Abundance_MGEs_subtypes_"+cell_suffix+".txt"),sep="\t",index=None)
    #cal rpkm
    name_list_rpkm=[]
    for i in file_name_base:
        name_list_rpkm.append(i+"_MGE_rpkmAbu_tmp.txt")
    init=0
    df_main=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_rpkm[0]),
                        sep="\t", header=None)
    df_main.columns=["type",name_list_rpkm[0]]
    for i,name in enumerate(name_list_rpkm):
        if i < len(name_list_16S)-1:
            init+=1
            if name_list_rpkm[init]:
                df_2=pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",name_list_rpkm[init]),
                                 sep="\t", header=None)
                df_2.columns=["type",name_list_rpkm[init]]
                df_main=pd.merge(df_main,df_2,left_on="type",right_on="type",how="outer")
    #save rpkm subtype abu
    df_main.to_csv(os.path.join(
                    input_dir,project_prefix,
                        "CompRanking_result",
                            project_prefix+"_Abundance_MGEs_subtypes_rpkg.txt"),sep="\t",index=None)
    
    
    #concat Sub ARG Name result
    #os.path.join(input_dir,project_prefix,"CompRanking_result",i+"_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt"
    # 假设样本文件存放在当前目录下
    directory = os.path.join(input_dir,project_prefix,"CompRanking_result")
    sample_files = [f for f in os.listdir(directory) if f.endswith('_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt')]

    # 创建一个字典来保存 ARG_name 到 Class 的映射
    arg_class_mapping = {}
    find_db_mapping= {}
    mge_type_mapping= {}
    # 创建一个空的 DataFrame 来存放合并后的数据
    merged_df = pd.DataFrame()

    for sample_file in sample_files:
        # 获取样本名
        sample_name = sample_file.split('_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt')[0]
        
        # 读取当前样本的 TSV 文件
        sample_df = pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",sample_file), sep='\t')
        
        # 检查列名
        if 'ARG_name' not in sample_df.columns or 'Class' not in sample_df.columns or 'Value' not in sample_df.columns or 'find_db' not in sample_df.columns:
            raise ValueError(f"File {sample_file} does not contain the required columns: 'ARG_name', 'Class', 'find_db', 'Value'")
        
        # 更新字典 arg_class_mapping
        for _, row in sample_df.iterrows():
            if row['ARG_name'] not in arg_class_mapping:
                arg_class_mapping[row['ARG_name']] = row['Class']
        
        # 更新字典 find_db_mapping
        for _, row in sample_df.iterrows():
            if row['ARG_name'] not in find_db_mapping:
                find_db_mapping[row['ARG_name']] = row['find_db']
        
        # 更新字典 mge_type_mapping
        for _, row in sample_df.iterrows():
            if not row['MGE_type'].startswith("ambiguous") and not row['MGE_type'].startswith("unclassified"):
                if row['ARG_name'] not in mge_type_mapping :
                    mge_type_mapping[row['ARG_name']] = row['MGE_type']
                else:
                    if row['MGE_type'] in mge_type_mapping[row['ARG_name']]:
                        continue
                    else:
                        mge_type_mapping[row['ARG_name']] = mge_type_mapping[row['ARG_name']]+"/"+row['MGE_type']
        
        # 将当前样本的值加入到新的列中
        sample_value_df = sample_df[['ARG_name', 'Value']]
        sample_value_df.columns = ['ARG_name', sample_name]
        
        # 设置 ARG_name 为索引以便合并
        sample_value_df.set_index('ARG_name', inplace=True)
        
        # 合并数据
        if merged_df.empty:
            merged_df = sample_value_df
        else:
            merged_df = merged_df.join(sample_value_df, how='outer')
    
    # 对 mge_type_mapping 中的值进行排序
    for arg_name in mge_type_mapping:
        if pd.isna(mge_type_mapping[arg_name]) or mge_type_mapping[arg_name] == '':
            mge_type_mapping[arg_name] = '-'
        else:
            mge_types = mge_type_mapping[arg_name].split('/')
            mge_type_mapping[arg_name] = '/'.join(sorted(mge_types))
            
    # 合并后重新添加 Class 列
    merged_df.reset_index(inplace=True)
    merged_df['Class'] = merged_df['ARG_name'].map(arg_class_mapping)
    merged_df['Database'] = merged_df['ARG_name'].map(find_db_mapping)
    merged_df['MGE_type'] = merged_df['ARG_name'].map(mge_type_mapping)
   
    # 调整列顺序
    merged_df = merged_df[['ARG_name', 'Class','Database', 'MGE_type'] + [col for col in merged_df.columns if col not in ['ARG_name', 'Class','Database', 'MGE_type']]]
    merged_df_fillZero=merged_df.copy()
    merged_df_fillZero['MGE_type'].fillna('Unknown', inplace=True)
    merged_df_fillZero.fillna(0, inplace=True)
    # 重置索引并保存到新的 TSV 文件
    merged_df.to_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",project_prefix+'_merged_samples_with_class.tsv'), sep='\t', index=False)    
    merged_df_fillZero.to_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",project_prefix+'_merged_samples_with_class_fillZero.tsv'), sep='\t', index=False)
    print("合并后的文件已保存为 merged_samples.tsv")    
    
    os.system("rm " + os.path.join(input_dir,project_prefix,"CompRanking_result/*tmp*"))

    
#python Genecal.py -i /lomi_home/gaoyang/software/CompRanking/tmp_test -p DSR -t 20
        
        
    