#!/usr/bin/env python
# title             :combined.py
# description       :process hmm output for virulence prediction
# author            :Gaoyang Luo
# date              :202209017
# version           :1.0
# usage             :import AMRcombined
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import path


class AMRCombined():
    #preset feature
    # def __init__(self, rgiInputDir, dvfInputDir,plasflowInputDir,seekerInputDir, VFInputDir, contigIndex):
    #     self.rgi_input_dir=rgiInputDir
    #     self.dvf_input_dir=dvfInputDir
    #     self.plasflow_input_dir=plasflowInputDir
    #     self.seeker_input_dir=seekerInputDir
    #     self.virulence_dir=VFInputDir
    #     self.contigIndex=contigIndex

    
    
    def AMR_combined(self, input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_plasflow,seeker_table,filebase):
        #open RGI results
        df_RGI=pd.read_csv(input_rgi, sep="\t")
        df_RGI=df_RGI.fillna("-")
        df_RGI=df_RGI[df_RGI.Nudged == "-"]
        df_RGI_tmp=df_RGI.loc[:,["ORF_ID","Best_Hit_ARO","AMR Gene Family","Resistance Mechanism","SNPs_in_Best_Hit_ARO"]]
        df_RGI_tmp["Best_Hit_ARO"]=df_RGI_tmp["Best_Hit_ARO"].str.upper()
        
        #load deeparg results
        df_deeparg_sure=pd.read_csv(input_deeparg, sep="\t")
        df_deeparg_sure1=df_deeparg_sure.loc[:,["read_id","#ARG", "predicted_ARG-class"]]
        df_RGI_deeparg=pd.merge(df_deeparg_sure1,df_RGI_tmp,left_on="read_id",right_on="ORF_ID",how="outer")
        df_RGI_deeparg=df_RGI_deeparg.fillna("-")
        df_RGI_deeparg["ARG_prediction"]="-" #add empty column
        df_RGI_deeparg["Database"]="-"
        
        #orfid_filter
        #merge RGI and DeepARG
        for index, name in df_RGI_deeparg.iterrows():
            if df_RGI_deeparg["read_id"][index]=="-":
                df_RGI_deeparg["read_id"][index]=df_RGI_deeparg["ORF_ID"][index]
            if df_RGI_deeparg["#ARG"][index]=="-":#DeepARG not found, only RGI found
                df_RGI_deeparg["#ARG"][index]=df_RGI_deeparg["Best_Hit_ARO"][index]
                df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
                df_RGI_deeparg["Database"][index]="RGI"
            if df_RGI_deeparg["#ARG"][index]!="-" and df_RGI_deeparg["Best_Hit_ARO"][index]!="-":#both found
                if df_RGI_deeparg["#ARG"][index]== df_RGI_deeparg["Best_Hit_ARO"][index]:#identical
                    df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
                    df_RGI_deeparg["Database"][index]="DeepARG/RGI"
                if df_RGI_deeparg["#ARG"][index] != df_RGI_deeparg["Best_Hit_ARO"][index]:#differ
                    df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index] + "/" + df_RGI_deeparg["Best_Hit_ARO"][index]
                    df_RGI_deeparg["Database"][index]="DeepARG/RGI"
            if df_RGI_deeparg["#ARG"][index]!="-" and df_RGI_deeparg["Best_Hit_ARO"][index]=="-":
                df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
                df_RGI_deeparg["Database"][index]="DeepARG"
        
        df_RGI_deeparg.sort_values('read_id') #sort
        #rename columns
        df_RGI_deeparg1=df_RGI_deeparg.drop(["#ARG","ORF_ID","Best_Hit_ARO","AMR Gene Family","Resistance Mechanism"],axis=1,inplace=False)
        df_RGI_deeparg1=df_RGI_deeparg1[["read_id","ARG_prediction","predicted_ARG-class","SNPs_in_Best_Hit_ARO","Database"]]
        
        ####merge SARG####
        #open SARG
        df_SARG=pd.read_csv(input_SARG, sep="\t", header=0,index_col=0)
        df_SARG["class"]=df_SARG["class"].str.upper()
        df_SARG.columns=["ORF","query","class","Phenotype","ARG_rank"]
        
        #merge SARG and RGI_DeepARG
        df_RGI_Deeparg_SARG=pd.merge(df_RGI_deeparg1,df_SARG,left_on="read_id",right_on="ORF",how="outer")
        df_RGI_Deeparg_SARG=df_RGI_Deeparg_SARG.fillna("-")
        
        #orfid_filter
        #merge SARG with RGI_DeepARG
        for index, name in df_RGI_Deeparg_SARG.iterrows():
            #add orf from SARG found
            if df_RGI_Deeparg_SARG["read_id"][index]=="-":
                df_RGI_Deeparg_SARG["read_id"][index]=df_RGI_Deeparg_SARG["ORF"][index]
            #add only SARG found
            if df_RGI_Deeparg_SARG["ARG_prediction"][index]=="-":
                df_RGI_Deeparg_SARG["ARG_prediction"][index]=df_RGI_Deeparg_SARG["class"][index] #fill arg_pred
                df_RGI_Deeparg_SARG["Database"][index]="SARG" # fill database
                df_RGI_Deeparg_SARG["predicted_ARG-class"][index]=df_RGI_Deeparg_SARG["Phenotype"][index]#fill arg class
            #add both found
            if df_RGI_Deeparg_SARG["ARG_prediction"][index]!="-" and df_RGI_Deeparg_SARG["class"][index]!="-":#both found
                #check if indentity
                if df_RGI_Deeparg_SARG["ARG_prediction"][index] == df_RGI_Deeparg_SARG["class"][index]:#identical
                    df_RGI_Deeparg_SARG["Database"][index]="DeepARG/RGI/SARG"
                if df_RGI_Deeparg_SARG["ARG_prediction"][index] != df_RGI_Deeparg_SARG["class"][index]:#differ
                    df_RGI_Deeparg_SARG["ARG_prediction"][index]=df_RGI_Deeparg_SARG["ARG_prediction"][index] + "/" + df_RGI_Deeparg_SARG["class"][index]
                    df_RGI_Deeparg_SARG["Database"][index]=df_RGI_Deeparg_SARG["Database"][index] + "/" + "SARG"
                #check class identity
                if df_RGI_Deeparg_SARG["predicted_ARG-class"][index] != df_RGI_Deeparg_SARG["Phenotype"][index]:
                    df_RGI_Deeparg_SARG["predicted_ARG-class"][index]=df_RGI_Deeparg_SARG["predicted_ARG-class"][index] + "/"+df_RGI_Deeparg_SARG["Phenotype"][index]
        
        #cut off the edges
        df_RGI_Deeparg_SARG_final=df_RGI_Deeparg_SARG.drop(["ORF","query","class","Phenotype"],axis=1, inplace=False)
        
        ####add contig ID####
        #load index
        df_contig_ID=pd.read_csv(input_contig_ID, sep="\t", header=None)
        
        #merge contig and RGI_DeepARG
        df_RGI_Deeparg_contig=pd.merge(df_contig_ID,df_RGI_Deeparg_SARG_final,left_on=1,right_on='read_id',how="outer")
        df_RGI_Deeparg_contig=df_RGI_Deeparg_contig.fillna("-")
        
        #drop duplicates
        df_RGI_Deeparg_contig=df_RGI_Deeparg_contig.drop(["read_id"],axis=1,inplace=False)
        
        ####merge MGE####
        #open df_dvf
        df_dvf=pd.read_csv(input_dvf,sep="\t",header=0)
        df_dvf["name"]=df_dvf["name"].str.split(" ", expand=True)[0]
        df_dvf["dvf_pred"]="-"
        
        #determine dvf_pred
        for index, name in df_dvf.iterrows():
            if (df_dvf["score"][index]>=0.7):
                if (df_dvf["pvalue"][index] <= 0.05):
                    df_dvf["dvf_pred"][index]="phage"
        df_dvf=df_dvf[["name","dvf_pred"]]
        
        #open df_plasflow
        df_plasflow=pd.read_csv(input_plasflow,sep="\t",header=0,index_col=0)
        df_plasflow["label"]=df_plasflow["label"].str.split(".", expand=True)[0]
        df_plasflow=df_plasflow[["contig_name","label"]]
        
        #merge plasflow and contigID
        df_plasflow_contig=pd.merge(df_contig_ID,df_plasflow,left_on=0,right_on="contig_name",how="left")

        # merge plasflow and dvf
        df_plasflow_dvf_contig=pd.merge(df_plasflow_contig,df_dvf,left_on=0,right_on="name",how="left")
        df_plasflow_dvf_contig=df_plasflow_dvf_contig.fillna("-")
        df_plasflow_dvf_contig=df_plasflow_dvf_contig.drop(["contig_name","name"],axis=1,inplace=False)
        df_plasflow_dvf_contig["tmp"]="-"
        
        #Give dvf and plasflow combined prediction
        #ambiguous (plasmid/phage)
        for i, name in df_plasflow_dvf_contig.iterrows():
            if df_plasflow_dvf_contig["dvf_pred"][i] == "phage":
                if df_plasflow_dvf_contig["label"][i] == "unclassified":
                    df_plasflow_dvf_contig["tmp"][i] = "phage"
                else:
                    df_plasflow_dvf_contig["tmp"][i] = "ambiguous " + "(" + df_plasflow_dvf_contig["label"][i] + "/" + "phage" + ")"
            else:
                df_plasflow_dvf_contig["tmp"][i] = df_plasflow_dvf_contig["label"][i]
        
        #check results
        df_plasflow_dvf_contig=df_plasflow_dvf_contig.drop(["label","dvf_pred"],axis=1,inplace=False)
        df_plasflow_dvf_contig.columns=["Contig", "ORF_ID", "MGE_prediction"]
        
        ##merge seeker results
        
        PF_res=df_plasflow_dvf_contig
        seeker_res=pd.read_csv(seeker_table, sep="\t", header=None)
        #Write seeker result into dictionary
        dic={} 
        for i, names in seeker_res.iterrows():
            keys=names[0]
            values=names[1]
            dic[keys]=values
        PF_res["seeker_res"]="0"
        #Concatenate seeker resulte to pathofact table
        seeker_res_list=[]
        for j, mcp_rows in PF_res.iterrows():
            tmp=mcp_rows["Contig"]
            tar=dic[tmp]
            seeker_res_list.append(tar)
        PF_res["seeker_res"]=seeker_res_list
        
        #Judge MGE according to peram rules
        peram_MGE_pred=[]
        for c, ult in PF_res.iterrows():
            if ult["MGE_prediction"] == "ambiguous (plasmid/phage)": 
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("plasmid")
            elif ult["MGE_prediction"] == "ambiguous (chromosome/phage)":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("chromosome")
            elif ult["MGE_prediction"] == "phage":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("unclassified")
            elif ult["MGE_prediction"] == "chromosome":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("ambiguous (chromosome/phage)")
                else:
                    peram_MGE_pred.append("chromosome")
            elif ult["MGE_prediction"] == "plasmid":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("ambiguous (plasmid/phage)")
                else:
                    peram_MGE_pred.append("plasmid")
            elif ult["MGE_prediction"] == "unclassified":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("unclassified")

        #Write peram MGE prediction into final table
        PF_res["CompRanking_MGE_prediction"]=peram_MGE_pred
        PF_res=PF_res.drop(["MGE_prediction","seeker_res"],axis=1,inplace=False)
        
        ####merge ARG and MGE####
        #merge ARG and MGE
        df_AMR_contig=pd.merge(df_RGI_Deeparg_contig,PF_res,left_on=0,right_on="Contig",how="left")
        df_AMR_contig=df_AMR_contig.drop(["Contig","ORF_ID"], axis=1, inplace=False)
        df_AMR_contig=df_AMR_contig.drop_duplicates()
        df_AMR_contig.columns=["Contig","ORF_ID","ARG_prediction","ARG_class","SNPs","Database","ARG_rank","CompRanking_MGE_prediction"]
        
        ####processing MobileOG-db
        #load MobileOG result
        df_MobileOG=pd.read_csv(input_mobileOG, sep="\t", header=None)
        df_MobileOG.columns=['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        #处理分列
        df_MobileOG_tmp=df_MobileOG["sub_id"].str.split("|",expand=True)
        df_MobileOG_tmp=df_MobileOG_tmp[[0,1,3,4,5]]
        df_MobileOG_tmp.columns=["mobileOG_ID", "Gene_Name", "Taxonomy","Major_Category","MGE_Database"]     
        #concat
        df_MobileOG_concat=pd.concat((df_MobileOG,df_MobileOG_tmp),axis=1)
        df_MobileOG_concat=df_MobileOG_concat.drop(["sub_id","mismat","gapOpens","qStart","qEnd","sStart","sEnd","eval","bit"],axis=1, inplace=False)  
        #load mobileOG structure
        input_mobileOG_structure="/lomi_home/gaoyang/software/CompRanking/databases/MobileOG-db/MobileOG-db_structure.tsv"
        df_mobileOG_structure=pd.read_csv(input_mobileOG_structure, sep="\t", header=0) 
        #creat dic index of structure
        mobilOG_structure_dic={}
        for i, name in df_mobileOG_structure.iterrows():
            mobilOG_structure_dic[df_mobileOG_structure["mobileOG_ID"][i]]=df_mobileOG_structure["length"][i]
        #cal identity and cov
        #identity60
        df_MobileOG_concat["coverage"]="-"
        df_MobileOG_concat_iden60=df_MobileOG_concat[df_MobileOG_concat.identity > 60]    
        #cal cov
        empty=[]
        for index, name in df_MobileOG_concat_iden60.iterrows():
            empty.append(df_MobileOG_concat_iden60["alignLen"][index]/mobilOG_structure_dic[df_MobileOG_concat_iden60["mobileOG_ID"][index]])
        df_MobileOG_concat_iden60["coverage"]=empty
        #coverage>0.9
        df_MobileOG_concat_iden60_cov90=df_MobileOG_concat_iden60[df_MobileOG_concat_iden60.coverage > 0.9]
        #merge mobileOG 
        df_AMR_annotate_contig=pd.merge(df_AMR_contig,df_MobileOG_concat_iden60_cov90,left_on="ORF_ID",right_on="id",how="left")
        df_AMR_annotate_contig=df_AMR_annotate_contig.drop(["id","identity","alignLen"],axis=1, inplace=False)
        df_AMR_annotate_contig=df_AMR_annotate_contig.fillna("-")
        
                                
        ####save####
        # df_AMR_annotate_contig.to_csv("/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/CompRanking_ERR1191817_AMR_result.csv",sep="\t",index=0)
        AMRfile=df_AMR_annotate_contig.to_csv(output + "/CompRanking" + filebase + "_AMR_prediction2.tsv", sep="\t", index=0)
        
        return AMRfile


if __name__ == "__main__":
    import pandas as pd
    import re
    import glob
    import os
    import path
    #gloab settings
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    output=os.path.join(input_dir,"CompRanking/CompRanking_result")
    project_prefix="CompRanking"
    
    #AMR combine
    # input_rgi="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/RGI/ERR1191817.contigs_5M_contigs.RGI.out.txt"
    # input_deeparg="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/DeepARG/ERR1191817.contigs_5M_contigs_DeepARG.out.mapping.ARG"
    
    # input_SARG="../test_SARGrank_Protein_Result_tmp.txt"
    # input_contig_ID="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.index"
    # input_dvf="../test/CompRanking/CompRanking_intermediate/MGE/DVF/ERR1191817.contigs_5M_contigs.fa_gt500bp_dvfpred.txt"
    # input_plasflow="../test/CompRanking/CompRanking_intermediate/MGE/Plasflow/ERR1191817.contigs_5M_contigs_plasflow_predictions.tsv"
    # seeker_table="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/MGE/Seeker/seeker_ERR1191817.contigs_5M_contigs_output.txt"
    # input_mobileOG="/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/MGE/MobileOG/ERR1191817.contigs_5M_contigs_mobileOG_diamond.txt"
    
    
    a=AMRCombined()
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    for i in file_name_base:
        input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
        input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs_DeepARG.out.mapping.ARG")
        input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein80_Result.tsv")
        input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
        input_dvf=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DVF",i+"_5M_contigs.fa_gt500bp_dvfpred.txt")
        input_plasflow=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Plasflow",i+"_5M_contigs_plasflow_predictions.tsv")
        seeker_table=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Seeker","seeker_"+i+"_5M_contigs_output.txt")
        input_mobileOG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/MobileOG",i+"_5M_contigs_mobileOG_diamond.txt")

        a.AMR_combined(input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_plasflow,seeker_table,i)
    
    
    
    
    
        
   
        
        
        