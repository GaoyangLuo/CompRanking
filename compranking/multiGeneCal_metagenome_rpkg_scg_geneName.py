#!/usr/bin/env python
# title             :multiGeneCal_metagenome_rpkg_scg_geneName.py -> Gene_cal.ipynb
# description       :calculate relative abundance of each sub genes
#                    This version is used to calculate rpkg abundance (beta version)
#                    and single copy gene (cell copy)
# author            :Gaoyang Luo
# date              :20240701
# version           :1.0
# usage             :python compranking/multiGeneCal_metagenome_rpkg_scg_geneName.py -i <input_dir>
#                                                                        -p <project_prefix>
#                                                                        -n <normalization_base> #AGS or scg
#                                                                        -t <threads>
#                                                                        -d <pth2KK2db>
#                                                                        -k # Enable Kraken2
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
parser.add_option("-k", "--run_kk2", action = "store_true", dest = "run_kk2", default=False,
                  help = "Enable Kraken2 taxonomy annotation (default: False)")
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
run_kk2=options.run_kk2
output=os.path.join(input_dir,project_prefix,"CompRanking_result")
#default parameters
if (options.project_prefix is None):
    project_prefix="CompRanking" #default project name
if (options.threads is None):
    threads = "24" #default threads
if (options.config_file is None):
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"test_yaml.yaml") #"../test_yaml.yaml" default config_file path
if (options.database is None):
    database = "/g/data1b/mp96/JXT/db/kraken2db/minikraken2_v1_8GB"#default config_file path
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


###update

def RB_gene_sum(DB_deepARG_length, DB_SARG_length, DB_MobileOG_length, 
                input_AMR_sum, input_kk2, input_deeparg_sure,
                input_rgi, input_SARG, input_scg, input_rpkm, input_indexFile, input_indexFile2, genome_length, filebase, normalization_base="AGS"):
    
    if normalization_base == "AGS":
        cell_suffix = "Cell"
    elif normalization_base == "16S":
        cell_suffix = "16S"   
    else:
        cell_suffix = "Unknown"
    
    # ==========================================
    # 1. Load AMR results
    # ==========================================
    # print(f"Loading AMR file: {input_AMR_sum}")
    df_AMR_sum = pd.read_csv(input_AMR_sum, sep="\t", header=0)

    # Ensure ORF_ID is a string and strip potential whitespace
    df_AMR_sum["ORF_ID"] = df_AMR_sum["ORF_ID"].astype(str).str.strip()

    df_AMR_hit = df_AMR_sum[df_AMR_sum.ARG_prediction != "-"]
    df_AMR_hit1 = df_AMR_hit[["ORF_ID", "ARG_prediction", "ARG_class", "Database", "CompRanking_MGE_prediction"]].copy()
    df_AMR_hit1["db_final"] = df_AMR_hit1["Database"].str.split("/", expand=True)[0]

    Record_db_orf = {}
    Record_ARG_name_orf = {}
    for i, name in df_AMR_hit1.iterrows():
        # Force conversion to string
        orf_key = str(name["ORF_ID"])
        Record_db_orf.setdefault(orf_key, str(name["db_final"]))
        Record_ARG_name_orf.setdefault(
            orf_key,
            [str(name["ARG_prediction"]), str(name["ARG_class"]), str(name["CompRanking_MGE_prediction"])]
        )

    # ==========================================
    # 2. Set gene_length
    # ==========================================
    if normalization_base == "AGS":
        gene_length = genome_length
    else:
        gene_length = 1550

    # ==========================================
    # 3. Load RPKM file (key modification: ID cleanup)
    # ==========================================
    print(f"Loading RPKM file: {input_rpkm}")

    try:
        # Try to read the file; if the format is complex, the header setting may need adjustment
        df_rpkm = pd.read_csv(input_rpkm, sep="\t", header=4)
    except Exception as e:
        print(f"Error reading RPKM file: {e}")
        return None

    array_rpkm = np.array(df_rpkm)

    rpkm_dic = {}
    rpkm_len_dic = {}

    # --- Build RPKM dictionary ---
    for i in array_rpkm:
        # ID cleanup: take the first column, convert to string, and keep only the part before whitespace
        raw_id = str(i[0])
        orf_id = raw_id.split()[0].strip()

        try:
            length = float(i[1])  # Assume column 2 is length
            reads = float(i[4])   # Assume column 5 is read count
        except (ValueError, IndexError):
            continue

        rpkm_dic[orf_id] = reads
        rpkm_len_dic[orf_id] = length

    # ==========================================
    # 4. Load Index file (for ID conversion)
    # ==========================================
    print(f"Loading Index file: {input_indexFile2}")
    df_index2 = pd.read_csv(input_indexFile2, sep="\t", header=None)
    array_indexFile2 = np.array(df_index2)

    index_dic2 = {}
    for i in array_indexFile2:
        # Assume column 1 is the standard ID used in RPKM, and column 2 is the original ID used in AMR/SCG
        orf_id_standard = str(i[0]).split()[0].strip()
        raw_id_input = str(i[1]).strip()

        # Build mapping: original ID -> standard ID
        index_dic2.setdefault(raw_id_input, orf_id_standard)

    # ==========================================
    # 5. Process SCG and calculate AGS denominator (excluding zeros)
    # ==========================================
    df_scg = pd.read_csv(input_scg, sep="\t", header=None)
    df_scg.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']

    # Filter high-quality alignments
    df_scg_iden50 = df_scg[df_scg.identity > 30]
    df_scg_len = df_scg_iden50[df_scg_iden50.alignLen > 0]
    df_scg_bit50 = df_scg_len[df_scg_len.bit > 50]

    scg_depths = []
    unique_scg_ids = df_scg_bit50['id'].unique()

    print(f"Processing {len(unique_scg_ids)} unique SCG hits...")

    ignored_zeros = 0  # Counter: record how many entries were excluded because they were 0

    for scg_id in unique_scg_ids:
        # Clean ID
        scg_id_str = str(scg_id).split()[0].strip()

        # ID mapping conversion
        if scg_id_str in index_dic2:
            target_orf_id = index_dic2[scg_id_str]
        else:
            target_orf_id = scg_id_str

        # Get value
        if target_orf_id in rpkm_dic:
            val = rpkm_dic[target_orf_id]

            # Remove RPKM values equal to 0
            # Only detected SCGs are used to calculate average depth
            if val > 0:
                # [update] use RPKM directly, not divided by length
                depth = val
                scg_depths.append(depth)
            else:
                ignored_zeros += 1

    print(f"Found {len(scg_depths)} valid SCGs with non-zero RPKM.")
    print(f"Ignored {ignored_zeros} SCGs with 0 RPKM.")

    if len(scg_depths) > 0:
        # Calculate average value
        num_scg = sum(scg_depths) / len(scg_depths)
        print(f"Calculated AGS factor (avg SCG RPKM): {num_scg}")

        # Sanity check
        if num_scg < 1.0:
            print("Warning: AGS factor is still < 1.0. This suggests very low sequencing depth.")
    else:
        # If all SCGs are 0, the sample may contain no bacteria, or sequencing depth may be extremely low
        num_scg = 1
        print("Warning: No non-zero SCGs found! Setting num_scg to 1 to avoid division by zero.")

    # ==========================================
    # 6. loading AMR results (DeepARG, RGI, SARG)
    # ==========================================
    df_deeparg_sure = pd.read_csv(input_deeparg_sure, sep="\t")
    df_RGI = pd.read_csv(input_rgi, sep="\t")
    df_RGI = df_RGI.fillna("-")
    df_SARG = pd.read_csv(input_SARG, sep="\t", header=0, index_col=0)
    df_SARG["class"] = df_SARG["class"].str.upper()
    df_SARG.columns = ["ORF", "query", "class", "Phenotype", "ARG_rank"]

    deepARG_output_dic = {}
    SARG_output_dic = {}
    for i, name in df_deeparg_sure.iterrows():
        deepARG_output_dic.setdefault(str(name["read_id"]), str(name["best-hit"].split("|")[0]))
    for i, name in df_SARG.iterrows():
        SARG_output_dic.setdefault(str(name["ORF"]), str(name["query"]))

    DB_deepARG_length_res = {}
    DB_SARG_length_res = {}
    DB_CARD_length_res = {}
    
    # fill length dic
    for i in deepARG_output_dic:
        if deepARG_output_dic[i] in DB_deepARG_length:
            DB_deepARG_length_res.setdefault(i, DB_deepARG_length[deepARG_output_dic[i]])
        else:
            DB_deepARG_length_res.setdefault(i, 1000)

    for i in SARG_output_dic:
        if SARG_output_dic[i] in DB_SARG_length:
            DB_SARG_length_res.setdefault(i, DB_SARG_length[SARG_output_dic[i]])
        else:
            DB_SARG_length_res.setdefault(i, 1000)

    for i, name in df_RGI.iterrows():
        rgi_orf = str(name["ORF_ID"]).split()[0].strip()
        DB_CARD_length_res.setdefault(rgi_orf, len(str(name["CARD_Protein_Sequence"])))
    
    # ==========================================
    # 7. cal ARG abundance
    # ==========================================
    abundance_arg_16S = 0
    abundance_arg_RPKM = 0
    RPKM_ARG = {}
    TAXO_ARG = {}
    RPKG_ARG_NAME = {}
    SCG_ARG_NAME = {} # update, used to restore ORF results of Cell Copy (AGS)
    
    for orf in Record_db_orf:
        find_db = ''
        if Record_db_orf[orf]:
            find_db = Record_db_orf[orf]
            ARG_name = Record_ARG_name_orf[orf][0].split("/")[0]
            ARG_class = Record_ARG_name_orf[orf][1].split("/")[0]
            MGE_type = Record_ARG_name_orf[orf][2]
            
            if orf in index_dic2:
                target_arg_id = index_dic2[orf]
            else:
                target_arg_id = orf 
            
            if target_arg_id in rpkm_dic:
                current_mapped_reads = rpkm_dic[target_arg_id]
            else:
                current_mapped_reads = 0 
            
            arg_len = 1000
            if find_db == "DeepARG" and orf in DB_deepARG_length_res:
                arg_len = DB_deepARG_length_res[orf]
            elif find_db == "RGI" and orf in DB_CARD_length_res:
                arg_len = DB_CARD_length_res[orf]
            elif find_db == "SARG" and orf in DB_SARG_length_res:
                arg_len = DB_SARG_length_res[orf]
            
            # cal abundance
            val_scg_norm = (current_mapped_reads / arg_len) / num_scg # Cell Copy
            val_rpkg = current_mapped_reads / (arg_len / 1000 * gene_length) # RPKM/RPKG
            
            abundance_arg_16S += val_scg_norm
            abundance_arg_RPKM += val_rpkg
            
            TAXO_ARG.setdefault(str(orf), float(val_scg_norm))
            RPKM_ARG.setdefault(str(orf), float(val_rpkg))
            
            # save RPKM results
            RPKG_ARG_NAME.setdefault(str(orf), [str(ARG_name), str(ARG_class), float(val_rpkg), str(find_db), str(MGE_type)])
            
            # [update] add Cell Copy (AGS) results
            # same structure, change val_rpkg to val_scg_norm
            SCG_ARG_NAME.setdefault(str(orf), [str(ARG_name), str(ARG_class), float(val_scg_norm), str(find_db), str(MGE_type)])

    print(f"The relative abundance of ARG by {cell_suffix}(AGS) is: {abundance_arg_16S}")
    print("The relative abundance of ARG by RPKG is: {}".format(abundance_arg_RPKM))
    
    # ==========================================
    # 8. Summary ARG Subtype
    # ==========================================
    RPKG_ARG_NAME_abundance = {}
    for orf, values in RPKG_ARG_NAME.items():
        arg_name = values[0]
        arg_class = values[1]
        value = values[2]
        find_db_val = values[3]
        mge_type_val = values[4]
        
        if arg_name in RPKG_ARG_NAME_abundance:
            RPKG_ARG_NAME_abundance[arg_name][1] += value
        else:
            RPKG_ARG_NAME_abundance[arg_name] = [arg_class, value, find_db_val, mge_type_val]
    
    abundance_ARG_subtype_16S = {}
    abundance_ARG_subtype_RPKM = {}
    for i, name in df_AMR_hit.iterrows():
        orf_id_key = str(name["ORF_ID"]).strip()
        if orf_id_key not in TAXO_ARG:
            continue
        
        subtype_name = name["ARG_class"].split("/")[0].split(":")[0].strip(";")
        
        if subtype_name in abundance_ARG_subtype_16S:
            abundance_ARG_subtype_16S[subtype_name] += TAXO_ARG[orf_id_key]
        else:
            abundance_ARG_subtype_16S[subtype_name] = float(TAXO_ARG[orf_id_key])
            
        if subtype_name in abundance_ARG_subtype_RPKM:
            abundance_ARG_subtype_RPKM[subtype_name] += RPKM_ARG[orf_id_key]
        else:
            abundance_ARG_subtype_RPKM[subtype_name] = float(RPKM_ARG[orf_id_key])

    # ==========================================
    # 9. MGE cal
    # ==========================================
    try:
        pass 
    except:
        pass

    df_MGE_hit = df_AMR_sum[df_AMR_sum.mobileOG_ID != "-"]
    DB_MobileOG_length_res = {}
    for i, name in df_MGE_hit.iterrows():
        if name["mobileOG_ID"] in DB_MobileOG_length:
            DB_MobileOG_length_res.setdefault(str(name["ORF_ID"]), DB_MobileOG_length[name["mobileOG_ID"]])
        else:
            DB_MobileOG_length_res.setdefault(str(name["ORF_ID"]), 1000)
    
    abundance_MGE_16S = 0
    abundance_MGE_RPKM = 0
    RPKM_MGE = {}
    TAXO_MGE = {}
    
    for orf_MGE in DB_MobileOG_length_res:
        if orf_MGE in index_dic2:
            target_mge_id = index_dic2[orf_MGE]
        else:
            target_mge_id = orf_MGE
            
        if target_mge_id in rpkm_dic:
            current_mapped_reads = rpkm_dic[target_mge_id]
        else:
            current_mapped_reads = 0
            
        mge_len = DB_MobileOG_length_res[orf_MGE]
        
        val_scg_norm = (current_mapped_reads / mge_len) / num_scg
        val_rpkg = current_mapped_reads / (mge_len / 1000 * gene_length)
        
        abundance_MGE_16S += val_scg_norm
        abundance_MGE_RPKM += val_rpkg
        TAXO_MGE.setdefault(str(orf_MGE), float(val_scg_norm))
        RPKM_MGE.setdefault(str(orf_MGE), float(val_rpkg))
    
    print(f"The relative abundance of MGE by {cell_suffix} is: {abundance_MGE_16S}")
    print("The relative abundance of MGE by RPKM is: {}".format(abundance_MGE_RPKM))
    
    # ==========================================
    # 10. MGE Subtype
    # ==========================================
    Insertion_Sequences_db=["ISFinder"]
    Integrative_Elements_db=["AICE","ICE","CIME","IME","immedb"]
    Plasmids_db=["COMPASS","Plasmid RefSeq"]
    Bacteriophages_db=["pVOG","GPD"]
    Multiple_db=["ACLAME", "Multiple"]

    abundance_MGE_subtype_16S = {}
    abundance_MGE_subtype_RPKM = {}
    
    for i, name in df_MGE_hit.iterrows():
        if name["ORF_ID"] not in TAXO_MGE:
            continue
        
        mge_db = name["MGE_Database"]
        subtype_key = ""
        
        if mge_db in Bacteriophages_db:
            subtype_key = "phage"
        elif mge_db in Plasmids_db:
            subtype_key = "plasmid"
        elif mge_db in Insertion_Sequences_db:
            subtype_key = "Insertion_Sequences"
        elif mge_db in Integrative_Elements_db:
            subtype_key = "Integrative_Elements"
        elif mge_db in Multiple_db:
            if name["Taxonomy"] == "phage":
                subtype_key = "phage"
            else:
                subtype_key = "plasmid"
        
        if subtype_key:
            if subtype_key in abundance_MGE_subtype_16S:
                abundance_MGE_subtype_16S[subtype_key] += TAXO_MGE[name["ORF_ID"]]
            else:
                abundance_MGE_subtype_16S[subtype_key] = float(TAXO_MGE[name["ORF_ID"]])
            
            if subtype_key in abundance_MGE_subtype_RPKM:
                abundance_MGE_subtype_RPKM[subtype_key] += RPKM_MGE[name["ORF_ID"]]
            else:
                abundance_MGE_subtype_RPKM[subtype_key] = float(RPKM_MGE[name["ORF_ID"]])

    result = [abundance_arg_16S, abundance_arg_RPKM, abundance_MGE_16S, abundance_MGE_RPKM]
    df_ARG_subtype_16S = pd.DataFrame(pd.Series(abundance_ARG_subtype_16S))
    df_ARG_subtype_RPKM = pd.DataFrame(pd.Series(abundance_ARG_subtype_RPKM))
    df_MGE_subtype_16S = pd.DataFrame(pd.Series(abundance_MGE_subtype_16S))
    df_MGE_subtype_RPKM = pd.DataFrame(pd.Series(abundance_MGE_subtype_RPKM))
    
    # update SCG_ARG_NAME
    return result, df_ARG_subtype_16S, df_ARG_subtype_RPKM, df_MGE_subtype_16S, df_MGE_subtype_RPKM, RPKG_ARG_NAME, RPKG_ARG_NAME_abundance, SCG_ARG_NAME

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
                                    i+"_5M_contigs.fna2faa.RGI.out.txt")
        input_SARG=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/AMR/ARGranking",
                                    i+"_SARGrank_Protein60_Result.tsv")
        input_deeparg_sure=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/AMR/DeepARG", 
                                    i+"_5M_contigs.fna2faa_DeepARG.out.mapping.ARG")
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
                                    i+"_5M_contigs.fna2faa_scg_Protein_diamond.txt")
        input_rpkm=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs/cov",
                                    i+"_5M_contigs_gene.rpkm")
        input_indexFile=os.path.join(input_dir,project_prefix,
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_5M_contigs.fna2faa.index")
        input_indexFile2=os.path.join(input_dir,project_prefix, #K123_1234_1 -> 1_1, orf_id -> orf ID
                                "CompRanking_intermediate/preprocessing/5M_contigs", 
                                    i+"_5M_contigs.fna2faa.orf2id.index")
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
        RPKG_ARG_NAME_abundance, \
        SCG_ARG_NAME= RB_gene_sum(DB_deepARG_length,
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
                                input_indexFile2,
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
        
    
    # 1. save original RPKM/RPKG results
        RPKG_ARG_NAME_tsv_data = "orf\tARG_name\tARG_class\tvalue\n"  # Add header
        for orf, values in RPKG_ARG_NAME.items():
            RPKG_ARG_NAME_tsv_data += f"{orf}\t{values[0]}\t{values[1]}\t{values[2]}\n" 
        
        output_dir_orf = os.path.join(input_dir, project_prefix, "CompRanking_result", "Abundance_orf_gene")
        os.makedirs(output_dir_orf, exist_ok=True)
        output_file_orf = os.path.join(output_dir_orf, i + "_Gene_Abundance_ORF_geneName_class_RPKG.txt")
        
        with open(output_file_orf, "w") as f:
            f.write(RPKG_ARG_NAME_tsv_data)
            
        # [update] add Cell Copy (AGS) results
        # use same format，data sourced from SCG_ARG_NAME (value is val_scg_norm)
        SCG_ARG_NAME_tsv_data = "orf\tARG_name\tARG_class\tvalue\n"
        for orf, values in SCG_ARG_NAME.items():
            SCG_ARG_NAME_tsv_data += f"{orf}\t{values[0]}\t{values[1]}\t{values[2]}\n"
            
        # 文件名加了 _SCG 区别
        output_file_orf_scg = os.path.join(output_dir_orf, i + "_Gene_Abundance_ORF_geneName_class_SCG.txt")
        
        with open(output_file_orf_scg, "w") as f:
            f.write(SCG_ARG_NAME_tsv_data)
        
        # Save Gene Name Abundance
        RPKG_ARG_NAME_abundance_tsv_data = "ARG_name\tClass\tfind_db\tMGE_type\tValue\n" 
        for arg, values in RPKG_ARG_NAME_abundance.items():
            RPKG_ARG_NAME_abundance_tsv_data += f"{arg}\t{values[0]}\t{values[2]}\t{values[3]}\t{values[1]}\n" 
        
        output_file_gene = os.path.join(input_dir, project_prefix, "CompRanking_result", i + "_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt")
        with open(output_file_gene, "w") as f:
            f.write(RPKG_ARG_NAME_abundance_tsv_data)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise ValueError("Write to summary abundance file failed...") from e
    
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
    kk2_script=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"scripts/kk2_run_single.sh")
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
    
    #run kranken2 (Only if -k / --run_kk2 is provided)
    if run_kk2:
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
    else:
        print("Skipping Kraken2 annotation as --run_kk2 was not specified.")
        
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
    # Assume sample files are stored in the current directory
    directory = os.path.join(input_dir,project_prefix,"CompRanking_result")
    sample_files = [f for f in os.listdir(directory) if f.endswith('_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt')]

    # Create a dictionary to save the mapping from ARG_name to Class
    arg_class_mapping = {}
    find_db_mapping= {}
    mge_type_mapping= {}
    # Create an empty DataFrame to store merged data
    merged_df = pd.DataFrame()

    for sample_file in sample_files:
        # Get sample name
        sample_name = sample_file.split('_Gene_Abundance_geneName_class_Cell(GE)_tmp.txt')[0]
        
        # Read the TSV file of the current sample
        sample_df = pd.read_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",sample_file), sep='\t')
        
        # Check column names
        if 'ARG_name' not in sample_df.columns or 'Class' not in sample_df.columns or 'Value' not in sample_df.columns or 'find_db' not in sample_df.columns:
            raise ValueError(f"File {sample_file} does not contain the required columns: 'ARG_name', 'Class', 'find_db', 'Value'")
        
        # Update dictionary arg_class_mapping
        for _, row in sample_df.iterrows():
            if row['ARG_name'] not in arg_class_mapping:
                arg_class_mapping[row['ARG_name']] = row['Class']
        
        # Update dictionary find_db_mapping
        for _, row in sample_df.iterrows():
            if row['ARG_name'] not in find_db_mapping:
                find_db_mapping[row['ARG_name']] = row['find_db']
        
        # Update dictionary mge_type_mapping
        for _, row in sample_df.iterrows():
            if not row['MGE_type'].startswith("ambiguous") and not row['MGE_type'].startswith("unclassified"):
                if row['ARG_name'] not in mge_type_mapping :
                    mge_type_mapping[row['ARG_name']] = row['MGE_type']
                else:
                    if row['MGE_type'] in mge_type_mapping[row['ARG_name']]:
                        continue
                    else:
                        mge_type_mapping[row['ARG_name']] = mge_type_mapping[row['ARG_name']]+"/"+row['MGE_type']
        
        # Add the current sample's values to a new column
        sample_value_df = sample_df[['ARG_name', 'Value']]
        sample_value_df.columns = ['ARG_name', sample_name]
        
        # Set ARG_name as index for merging
        sample_value_df.set_index('ARG_name', inplace=True)
        
        # Merge data
        if merged_df.empty:
            merged_df = sample_value_df
        else:
            merged_df = merged_df.join(sample_value_df, how='outer')
    
    # Sort values in mge_type_mapping
    for arg_name in mge_type_mapping:
        if pd.isna(mge_type_mapping[arg_name]) or mge_type_mapping[arg_name] == '':
            mge_type_mapping[arg_name] = '-'
        else:
            mge_types = mge_type_mapping[arg_name].split('/')
            mge_type_mapping[arg_name] = '/'.join(sorted(mge_types))
            
    # Re-add Class column after merging
    merged_df.reset_index(inplace=True)
    merged_df['Class'] = merged_df['ARG_name'].map(arg_class_mapping)
    merged_df['Database'] = merged_df['ARG_name'].map(find_db_mapping)
    merged_df['MGE_type'] = merged_df['ARG_name'].map(mge_type_mapping)
   
    # Adjust column order
    merged_df = merged_df[['ARG_name', 'Class','Database', 'MGE_type'] + [col for col in merged_df.columns if col not in ['ARG_name', 'Class','Database', 'MGE_type']]]
    merged_df_fillZero=merged_df.copy()
    merged_df_fillZero['MGE_type'].fillna('Unknown', inplace=True)
    merged_df_fillZero.fillna(0, inplace=True)
    # Reset index and save to a new TSV file
    merged_df.to_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",project_prefix+'_merged_samples_with_class.tsv'), sep='\t', index=False)    
    merged_df_fillZero.to_csv(os.path.join(input_dir,project_prefix,"CompRanking_result",project_prefix+'_merged_samples_with_class_fillZero.tsv'), sep='\t', index=False)
    print("合并后的文件已保存为 merged_samples.tsv")    
    
    os.system("rm " + os.path.join(input_dir,project_prefix,"CompRanking_result/*tmp*"))

    
#python Genecal.py -i /lomi_home/gaoyang/software/CompRanking/tmp_test -p DSR -t 20