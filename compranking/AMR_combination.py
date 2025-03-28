#!/usr/bin/env python
# title             :AMR_combined.py
# description       :combine AMR output
# author            :Gaoyang Luo
# date              :202209017: v1; 20240630:v2; 20250423:v3
# version           :1.0
# usage             :import AMRcombined
# required packages :pandas, numpy 
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import os
# import MOB_concat
# import Virulence_processing
# import summary_all
# from compranking import path


class AMRCombined():
    #preset feature
    # def __init__(self, rgiInputDir, dvfInputDir,plasflowInputDir,seekerInputDir, VFInputDir, contigIndex):
    #     self.rgi_input_dir=rgiInputDir
    #     self.dvf_input_dir=dvfInputDir
    #     self.plasflow_input_dir=plasflowInputDir
    #     self.seeker_input_dir=seekerInputDir
    #     self.virulence_dir=VFInputDir
    #     self.contigIndex=contigIndex

    def AMR_combined(self, input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_def,seeker_table,input_mobileOG):
        #######################Processing deepARG,RGI, SARG###################################
        #open RGI results
        df_RGI=pd.read_csv(input_rgi, sep="\t")
        df_RGI=df_RGI.fillna("-")
        # df_RGI=df_RGI[df_RGI.Nudged == "-"]
        df_RGI_tmp=df_RGI.loc[:,["ORF_ID","Best_Hit_ARO","Drug Class","Resistance Mechanism","SNPs_in_Best_Hit_ARO"]]
        df_RGI_tmp["Best_Hit_ARO"]=df_RGI_tmp["Best_Hit_ARO"].str.upper()
        
        #load deeparg results
        df_deeparg_sure=pd.read_csv(input_deeparg, sep="\t")
        df_deeparg_sure1=df_deeparg_sure.loc[:,["read_id","#ARG", "predicted_ARG-class"]]
        df_RGI_deeparg=pd.merge(df_deeparg_sure1,df_RGI_tmp,left_on="read_id",right_on="ORF_ID",how="outer")
        df_RGI_deeparg=df_RGI_deeparg.fillna("-")
        df_RGI_deeparg["ARG_prediction"]="-" #add empty column
        df_RGI_deeparg["Database"]="-"
        df_RGI_deeparg["orf_final"]="-"
        #orfid_filter
        #merge RGI and DeepARG
        # for index, name in df_RGI_deeparg.iterrows(): #RGI:ORF_ID deepARG:read_id
        #     if df_RGI_deeparg["read_id"][index]=="-": #not found in deepARG
        #         df_RGI_deeparg["read_id"][index]=df_RGI_deeparg["ORF_ID"][index] #fill deeparg with RGI
        #     if df_RGI_deeparg["#ARG"][index]=="-":#DeepARG not found, only RGI found
        #         df_RGI_deeparg["Database"][index]="RGI" # #ARG:deepARG  Best_Hit_ARO:RGI
        #         df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["Best_Hit_ARO"][index]
        #     #both found
        #     if df_RGI_deeparg["#ARG"][index]!="-" and df_RGI_deeparg["Best_Hit_ARO"][index]!="-":#both found
        #         if df_RGI_deeparg["#ARG"][index]== df_RGI_deeparg["Best_Hit_ARO"][index]:#identical
        #             df_RGI_deeparg["Database"][index]="DeepARG/RGI"
        #             df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]     
        #         if df_RGI_deeparg["#ARG"][index] != df_RGI_deeparg["Best_Hit_ARO"][index]:#differ
        #             df_RGI_deeparg["Database"][index]="DeepARG/RGI"
        #             df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index] + "/" + df_RGI_deeparg["Best_Hit_ARO"][index]
        #     #only DeepARG found        
        #     if df_RGI_deeparg["#ARG"][index]!="-" and df_RGI_deeparg["Best_Hit_ARO"][index]=="-": #only DeepARG found
        #         df_RGI_deeparg["Database"][index]="DeepARG"
        #         df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
        
        for index, name in df_RGI_deeparg.iterrows(): #RGI:ORF_ID deepARG:read_id
            # if df_RGI_deeparg["read_id"][index]=="-": #not found in deepARG
            #     df_RGI_deeparg["read_id"][index]=df_RGI_deeparg["ORF_ID"][index] #fill deeparg with RGI
            #DeepARG not found, only RGI found
            if df_RGI_deeparg["read_id"][index]=="-":
                if df_RGI_deeparg["ORF_ID"][index]!="-":#DeepARG not found, only RGI found
                    df_RGI_deeparg["Database"][index]="RGI" # #ARG:deepARG  Best_Hit_ARO:RGI
                    df_RGI_deeparg["orf_final"][index]=df_RGI_deeparg["ORF_ID"][index]
                    df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["Best_Hit_ARO"][index]
                    if len(df_RGI_deeparg["Drug Class"][index].split(";")) > 1: # more than one Drug class in RGI defined as multidrug
                        df_RGI_deeparg["predicted_ARG-class"][index]= "multidrug"
                    else:
                        df_RGI_deeparg["predicted_ARG-class"][index]=df_RGI_deeparg["Drug Class"][index].split(" ")[0].strip(";")
            #when DeepARG found   
            if df_RGI_deeparg["read_id"][index]!="-": 
                #but RGI not found
                if df_RGI_deeparg["ORF_ID"][index]=="-": #only DeepARG found
                    df_RGI_deeparg["Database"][index]="DeepARG"
                    df_RGI_deeparg["orf_final"][index]=df_RGI_deeparg["read_id"][index]
                    df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
                else:
                    #if df_RGI_deeparg["read_id"][index]!="-" and df_RGI_deeparg["ORF_ID"][index]!="-"
                    #both found RGI found too
                    if df_RGI_deeparg["#ARG"][index]== df_RGI_deeparg["Best_Hit_ARO"][index]:#check identical
                        df_RGI_deeparg["Database"][index]="DeepARG/RGI"
                        df_RGI_deeparg["orf_final"][index]=df_RGI_deeparg["read_id"][index]
                        df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index]
                    
                    if df_RGI_deeparg["#ARG"][index] != df_RGI_deeparg["Best_Hit_ARO"][index]:#check differ
                        df_RGI_deeparg["Database"][index]="DeepARG/RGI"
                        df_RGI_deeparg["orf_final"][index]=df_RGI_deeparg["read_id"][index]
                        df_RGI_deeparg["ARG_prediction"][index]=df_RGI_deeparg["#ARG"][index] + "/" + df_RGI_deeparg["Best_Hit_ARO"][index]
        
        #recover the ori id
        df_RGI_deeparg=df_RGI_deeparg[["orf_final","read_id","#ARG","predicted_ARG-class","ORF_ID","Best_Hit_ARO","Drug Class","Resistance Mechanism","SNPs_in_Best_Hit_ARO","ARG_prediction","Database"]]
        df_RGI_deeparg=df_RGI_deeparg.drop(["read_id"],axis=1,inplace=False)
        df_RGI_deeparg.columns=["read_id","#ARG","predicted_ARG-class","ORF_ID","Best_Hit_ARO","Drug Class","Resistance Mechanism","SNPs_in_Best_Hit_ARO","ARG_prediction","Database"]
        
        df_RGI_deeparg.sort_values('read_id') #sort
        #rename columns
        df_RGI_deeparg1=df_RGI_deeparg.drop(["#ARG","ORF_ID","Best_Hit_ARO","Drug Class","Resistance Mechanism"],axis=1,inplace=False)
        df_RGI_deeparg1=df_RGI_deeparg1[["read_id","ARG_prediction","predicted_ARG-class","SNPs_in_Best_Hit_ARO","Database"]]
        
        ####merge SARG####
        #open SARG
        df_SARG=pd.read_csv(input_SARG, sep="\t", header=0,index_col=0)
        df_SARG["class"]=df_SARG["class"].str.upper()
        df_SARG.columns=["ORF","query","class","Phenotype","ARG_rank"]
        
        #merge SARG and RGI_DeepARG
        df_RGI_Deeparg_SARG=pd.merge(df_RGI_deeparg1,df_SARG,left_on="read_id",right_on="ORF",how="outer")
        df_RGI_Deeparg_SARG=df_RGI_Deeparg_SARG.fillna("-")
        for i ,name in df_RGI_Deeparg_SARG.iterrows():
            if name["Phenotype"]=="macrolide-lincosamide-streptogramin":
                df_RGI_Deeparg_SARG["Phenotype"][i]="MLS"
        
        #orfid_filter
        #merge SARG with RGI_DeepARG
        for index, name in df_RGI_Deeparg_SARG.iterrows():
            #add orf from SARG found
            # if df_RGI_Deeparg_SARG["read_id"][index]=="-":
            #     df_RGI_Deeparg_SARG["orf_final"][index]=df_RGI_Deeparg_SARG["ORF"][index]
            
            #add both found
            #when both found 
            if df_RGI_Deeparg_SARG["read_id"][index]!="-" and df_RGI_Deeparg_SARG["ORF"][index]!="-":#both found
                df_RGI_Deeparg_SARG["Database"][index]=df_RGI_Deeparg_SARG["Database"][index]+"/" + "SARG"
                #check if indentical
                if df_RGI_Deeparg_SARG["ARG_prediction"][index] == df_RGI_Deeparg_SARG["class"][index]:#identical
                    continue
                else: #differ
                    df_RGI_Deeparg_SARG["ARG_prediction"][index]=df_RGI_Deeparg_SARG["ARG_prediction"][index] + "/" + df_RGI_Deeparg_SARG["class"][index]+"(SARG)"
                #check class identity
                if df_RGI_Deeparg_SARG["predicted_ARG-class"][index] != df_RGI_Deeparg_SARG["Phenotype"][index]:
                    df_RGI_Deeparg_SARG["predicted_ARG-class"][index]=df_RGI_Deeparg_SARG["predicted_ARG-class"][index] + "/"+df_RGI_Deeparg_SARG["Phenotype"][index]+"(SARG)"
            #add only SARG found
            if df_RGI_Deeparg_SARG["read_id"][index]=="-" and df_RGI_Deeparg_SARG["ORF"][index]!="-" :
                df_RGI_Deeparg_SARG["read_id"][index]=df_RGI_Deeparg_SARG["ORF"][index] #fill read_id with SARG hit
                df_RGI_Deeparg_SARG["ARG_prediction"][index]=df_RGI_Deeparg_SARG["class"][index] #fill arg_pred
                df_RGI_Deeparg_SARG["Database"][index]="SARG" # fill database
                df_RGI_Deeparg_SARG["predicted_ARG-class"][index]=df_RGI_Deeparg_SARG["Phenotype"][index]#fill arg class
        
        #cut off the edges
        df_RGI_Deeparg_SARG_final=df_RGI_Deeparg_SARG.drop(["ORF","query","class","Phenotype"],axis=1, inplace=False)
        
        #add contig ID#
        #load index
        df_contig_ID=pd.read_csv(input_contig_ID, sep="\t", header=None)
        
        #merge contig and RGI_DeepARG
        df_RGI_Deeparg_contig=pd.merge(df_contig_ID,df_RGI_Deeparg_SARG_final,left_on=1,right_on='read_id',how="outer")
        df_RGI_Deeparg_contig=df_RGI_Deeparg_contig.fillna("-")
        
        #drop duplicates
        df_RGI_Deeparg_contig=df_RGI_Deeparg_contig.drop(["read_id"],axis=1,inplace=False)
        
        ###################################Processing MGE#####################################
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
        
        #open df_def
        #DeepMicrobeFinder
        df_def=pd.read_csv(input_def,header=0,sep="\t")
        Eukaryote="Eukaryote"
        EukaryoteVirus="EukaryoteVirus"
        Plasmid="Plasmid"
        Prokaryote="Prokaryote"
        ProkaryoteVirus="ProkaryoteVirus"
        #check if empty
        if len(df_def["Sequence Name"]) == 0:
            df_def["max_idx"]=[]
            res=df_def[["Sequence Name","max_idx"]]
        else:
            df_def["Sequence Name"]=df_def["Sequence Name"].str.split(" ",expand=True)[0]
            df_def1=df_def.iloc[:,0:6]
            df_def1['max_idx']=df_def.iloc[:,1:6].idxmax(axis=1)
            res=df_def1[df_def1.max_idx == "Plasmid"]
            res=res[["Sequence Name","max_idx"]]
        #merge def and dvf
        df_def_contig=pd.merge(df_contig_ID,res,left_on=0,right_on="Sequence Name",how="left")
        df_def_contig=df_def_contig.fillna("-")
        df_def_dvf_contig=pd.merge(df_def_contig,df_dvf,left_on=0,right_on="name",how="left")
        df_def_dvf_contig=df_def_dvf_contig.fillna("-")
        df_def_dvf_contig=df_def_dvf_contig.drop(["Sequence Name","name"],axis=1,inplace=False)
        df_def_dvf_contig["tmp"]="-"
        df_def_dvf_contig.columns=[0,1,"label","dvf_pred","tmp"]
        
        # #open df_plasflow
        # df_plasflow=pd.read_csv(input_plasflow,sep="\t",header=0)#index_col=0
        # df_plasflow["label"]=df_plasflow["label"].str.split(".", expand=True)[0]
        # df_plasflow=df_plasflow[["contig_name","label"]]
        
        # #merge plasflow and contigID
        # df_plasflow_contig=pd.merge(df_contig_ID,df_plasflow,left_on=0,right_on="contig_name",how="left")

        # # merge plasflow and dvf
        # df_plasflow_dvf_contig=pd.merge(df_plasflow_contig,df_dvf,left_on=0,right_on="name",how="left")
        # df_plasflow_dvf_contig=df_plasflow_dvf_contig.fillna("-")
        # df_plasflow_dvf_contig=df_plasflow_dvf_contig.drop(["contig_name","name"],axis=1,inplace=False)
        # df_plasflow_dvf_contig["tmp"]="-"
        
        #Give dvf and plasflow combined prediction
        #ambiguous (plasmid/phage)
        # for i, name in df_plasflow_dvf_contig.iterrows():
        #     if df_plasflow_dvf_contig["dvf_pred"][i] == "phage":
        #         if df_plasflow_dvf_contig["label"][i] == "unclassified":
        #             df_plasflow_dvf_contig["tmp"][i] = "phage"
        #         elif df_plasflow_dvf_contig["label"][i] == "-":
        #             df_plasflow_dvf_contig["tmp"][i] = "phage"
        #         else:                
        #             df_plasflow_dvf_contig["tmp"][i] = "ambiguous " + "(" + df_plasflow_dvf_contig["label"][i] + "/" + "phage" + ")"
        #     else:
        #         if df_plasflow_dvf_contig["label"][i] == "-":
        #             df_plasflow_dvf_contig["tmp"][i] = "unclassified"
        #         else:
        #             df_plasflow_dvf_contig["tmp"][i] = df_plasflow_dvf_contig["label"][i]
        for i, name in df_def_dvf_contig.iterrows():
            if df_def_dvf_contig["dvf_pred"][i] == "phage":
                if df_def_dvf_contig["label"][i] == "unclassified":
                    df_def_dvf_contig["tmp"][i] = "phage"
                elif df_def_dvf_contig["label"][i] == "-":
                    df_def_dvf_contig["tmp"][i] = "phage"
                else:                
                    df_def_dvf_contig["tmp"][i] = "ambiguous " + "(" + df_def_dvf_contig["label"][i] + "/" + "phage" + ")"
            else:
                if df_def_dvf_contig["label"][i] == "-":
                    df_def_dvf_contig["tmp"][i] = "unclassified"
                else:
                    df_def_dvf_contig["tmp"][i] = df_def_dvf_contig["label"][i]
        #check results
        df_def_dvf_contig=df_def_dvf_contig.drop(["label","dvf_pred"],axis=1,inplace=False)
        df_def_dvf_contig.columns=["Contig", "ORF_ID", "MGE_prediction"]
        
        ##merge seeker results
        PF_res=df_def_dvf_contig
        seeker_res=pd.read_csv(seeker_table, sep="\t", header=None)
        #Write seeker result into dictionary
        dic={} 
        for i, names in seeker_res.iterrows():
            if  seeker_res[2][i] <= 0.6:
                keys=names[0]
                values="Bacteria"
                dic[keys]=values
            else:
                keys=names[0]
                values="Phage"
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
            if ult["MGE_prediction"] == "ambiguous (Plasmid/phage)": 
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("plasmid")
            # elif ult["MGE_prediction"] == "ambiguous (chromosome/phage)":
            #     if ult["seeker_res"] == "Phage":
            #         peram_MGE_pred.append("phage")
            #     else:
            #         peram_MGE_pred.append("chromosome")
            elif ult["MGE_prediction"] == "phage":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("phage")
                else:
                    peram_MGE_pred.append("unclassified")
            # elif ult["MGE_prediction"] == "chromosome":
            #     if ult["seeker_res"] == "Phage":
            #         peram_MGE_pred.append("ambiguous (chromosome/phage)")
            #     else:
            #         peram_MGE_pred.append("chromosome")
            elif ult["MGE_prediction"] == "Plasmid":
                if ult["seeker_res"] == "Phage":
                    peram_MGE_pred.append("ambiguous (Plasmid/phage)")
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
        
        ##############################merge ARG and MGE####################################
        #merge ARG and MGE
        df_AMR_contig=pd.merge(df_RGI_Deeparg_contig,PF_res,left_on=0,right_on="Contig",how="left")
        df_AMR_contig=df_AMR_contig.drop(["Contig","ORF_ID"], axis=1, inplace=False)
        df_AMR_contig=df_AMR_contig.drop_duplicates()
        df_AMR_contig.columns=["Contig","ORF_ID","ARG_prediction","ARG_class","SNPs","Database","ARG_rank","CompRanking_MGE_prediction"]
        
        ############################processing MobileOG-db##############################
        #load MobileOG result

        df_MobileOG=pd.read_csv(input_mobileOG, sep="\t", header=None)
        df_MobileOG.columns=['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
        #处理分列
        df_MobileOG_tmp=df_MobileOG["sub_id"].str.split("|",expand=True)
        df_MobileOG_tmp
        if df_MobileOG_tmp.shape[1] < 7:
            if df_MobileOG_tmp.shape[1] == 5:
                df_MobileOG_tmp[5]="-"
            if df_MobileOG_tmp.shape[1] == 6:
                df_MobileOG_tmp[6] = "-"
        # df_MobileOG_tmp[5]="-"
        df_MobileOG_tmp=df_MobileOG_tmp[[0,1,3,4,5]]
        df_MobileOG_tmp.columns=["mobileOG_ID", "Gene_Name", "Taxonomy","Major_Category","MGE_Database"]
  
        df_MobileOG_concat=pd.concat((df_MobileOG,df_MobileOG_tmp),axis=1)
        df_MobileOG_concat=df_MobileOG_concat.drop(["sub_id","mismat","gapOpens","qStart","qEnd","sStart","sEnd","eval","bit"],axis=1, inplace=False)  
        #load mobileOG structure
        # input_mobileOG_structure="/lomi_home/gaoyang/software/CompRanking/databases/MobileOG-db/MobileOG-db_structure.tsv"
        input_mobileOG_structure=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"databases/MobileOG-db_structure.tsv")
        df_mobileOG_structure=pd.read_csv(input_mobileOG_structure, sep="\t", header=0) 
        #creat dic index of structure
        mobilOG_structure_dic={}
        mobilOG_structure_db_dic={}
        for i, name in df_mobileOG_structure.iterrows():
            mobilOG_structure_dic[df_mobileOG_structure["mobileOG_ID"][i]]=df_mobileOG_structure["length"][i]
            mobilOG_structure_db_dic[df_mobileOG_structure["mobileOG_ID"][i]]=df_mobileOG_structure["Database"][i]
        #modify database
        db_tmp=[]
        for i, name in df_MobileOG_concat.iterrows():
            db_tmp.append(mobilOG_structure_db_dic[name["mobileOG_ID"]])
        df_MobileOG_concat["MGE_Database"]=db_tmp
        #cal identity and cov
        #identity60
        df_MobileOG_concat["coverage"]="-"
        df_MobileOG_concat_iden60=df_MobileOG_concat[df_MobileOG_concat.identity > 65]    
        #cal cov
        empty=[]
        for index, name in df_MobileOG_concat_iden60.iterrows():
            empty.append(df_MobileOG_concat_iden60["alignLen"][index]/mobilOG_structure_dic[df_MobileOG_concat_iden60["mobileOG_ID"][index]])
        df_MobileOG_concat_iden60["coverage"]=empty
        #coverage>0.5
        df_MobileOG_concat_iden60_cov50=df_MobileOG_concat_iden60[df_MobileOG_concat_iden60.coverage > 0.95]
        #merge mobileOG 
        df_AMR_annotate_contig=pd.merge(df_AMR_contig,df_MobileOG_concat_iden60_cov50,left_on="ORF_ID",right_on="id",how="left")
        df_AMR_annotate_contig=df_AMR_annotate_contig.drop(["id","identity","alignLen"],axis=1, inplace=False)
        df_AMR_annotate_contig=df_AMR_annotate_contig.fillna("-")

        
                                
        ####save####
        # df_AMR_annotate_contig.to_csv("/lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/AMR/CompRanking_ERR1191817_AMR_result.csv",sep="\t",index=0)
        # AMRfile=df_AMR_annotate_contig.to_csv(output + "/CompRanking_" + filebase + "_AMR_prediction.tsv", sep="\t", index=0)
        
        # return AMRfile
        return df_AMR_annotate_contig
 
    def refFilter(self,df_AMR_annotate_MOB_contig):
        # #according to MOB to filter ambiguous (plasmid/phage), if MOB == conj, ambiguous = plasmid
        # a=df_AMR_annotate_MOB_contig[df_AMR_annotate_MOB_contig.CompRanking_MGE_prediction == "ambiguous (plasmid/phage)"]
        # b=a[a.MOB == "mob_conj"]
        # index_list=[]
        # for i, name in b.iterrows():
        #     index_list.append(i)

        # for i in index_list:
        #     df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="plasmid"
    
        #Using MobileOG result filter
        Insertion_Sequences=["ISFinder"]
        Integrative_Elements=["AICE","ICE","CIME","IME","immedb"]
        Plasmids=["COMPASS","Plasmid RefSeq"]
        Bacteriophages=["pVOG","GPD"]
        Multiple=["ACLAME", "Multiple"]
        count_ref_plasmid=[]
        count_ref_phage=[]
        for i, name in df_AMR_annotate_MOB_contig.iterrows():
            if name["MGE_Database"] in Bacteriophages:
                df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="phage"
                count_ref_phage.append(df_AMR_annotate_MOB_contig["Contig"][i])
            elif name["MGE_Database"] in Plasmids:
                df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="plasmid"
                count_ref_plasmid.append(df_AMR_annotate_MOB_contig["Contig"][i])
            elif name["MGE_Database"] in Insertion_Sequences:
                df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="IS"
            elif name["MGE_Database"] in Integrative_Elements:
                df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="IE"
            elif name["MGE_Database"] in Multiple:
                if name["Taxonomy"] != "phage":
                    df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="plasmid"
                    # count_ref_plasmid.append(df_AMR_annotate_MOB_contig["Contig"][i])
                    continue
                elif name["Taxonomy"] == "phage":
                    df_AMR_annotate_MOB_contig["CompRanking_MGE_prediction"][i]="phage"  
                    count_ref_phage.append(df_AMR_annotate_MOB_contig["Contig"][i])
                else:
                    continue
            else:
                continue
        return df_AMR_annotate_MOB_contig


if __name__ == "__main__":
    import pandas as pd
    import re
    import os
    import path
    import datetime
    import summary_all
    import MOB_concat
    import Virulence_processing
    #gloab settings
    project_prefix="CompRanking"
    input_dir="/lomi_home/gaoyang/software/CompRanking/test"
    input_cpr_VF_sum="../databases/CompRanking_VirulenceDB/CompRanking_Virulence_Summary.csv" #fixed
    output=os.path.join(input_dir,project_prefix,"CompRanking_result")
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
    
    
    AMR_combine=AMRCombined()
    file_abs_path=path.file_abs_path_list_generation(input_dir)
    file_name_base = path.file_base_acquire(file_abs_path)
    
    start = datetime.datetime.now() #time start 
    for i in file_name_base:
        #set ARG&MGE input
        # input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.RGI.out.txt")
        input_rgi=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/RGI",i+"_5M_contigs.fna2faa.RGI.out.txt")
        # input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs_DeepARG.out.mapping.ARG")
        input_deeparg=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/DeepARG",i+"_5M_contigs.fna2faa.DeepARG.out.mapping.ARG")
        input_SARG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/AMR/ARGranking",i+"_SARGrank_Protein60_Result.tsv")
        # input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
        input_contig_ID=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.fna2faa.index")
        input_dvf=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/DVF",i+"_5M_contigs.fa_gt500bp_dvfpred.txt")
        input_plasflow=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Plasflow",i+"_5M_contigs_plasflow_predictions.tsv")
        seeker_table=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/Seeker","seeker_"+i+"_5M_contigs_output.txt")
        input_mobileOG=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/MobileOG",i+"_5M_contigs_mobileOG_diamond.txt")
        input_mob_conj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_Conj_plasmids_id_out")
        input_mob_unconj=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/MGE/plascad",i+"_5M_contigs_mob_unconj_plasmids_id_out")
        
        #set VF input
        input_contig=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/preprocessing/5M_contigs",i+"_5M_contigs.index")
        input_ERR_VFDB_output=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/VFDB",i+"_5M_contigs_VFDB.out") #5M_contigs_VFDB_setA1e-5.out
        input_patric=os.path.join(input_dir,project_prefix,"CompRanking_intermediate/Virulence/PATRIC",i+"_5M_contigs_PATRIC.out")
        
        #generate sum table without mob ref
        df_AMR_annotate_contig=AMR_combine.AMR_combined(input_rgi, input_contig_ID, input_deeparg, input_SARG,input_dvf, input_plasflow,seeker_table,input_mobileOG)
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
        
    end = datetime.datetime.now() #time end
    print(end-start)
    
    # #calculate
    # start = datetime.datetime.now() #time start
    
    
    
    
        
        
    
    
    
    
    
        
   
        
        
        