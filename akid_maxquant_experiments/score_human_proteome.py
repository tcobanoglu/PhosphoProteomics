# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:51:12 2022

@author: Tugce Su
"""

import pandas as pd
import csv

KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/gbm-akid-result.txt", delimiter=r"\s+")
KS_relations["Kinase"]= KS_relations["Kinase"].str.upper()

#to choose nwkin and psp kinase
KS_relations_nwkin= KS_relations.loc[KS_relations["Source"] == "NWK"]
KS_relations_psp= KS_relations.loc[KS_relations["Source"] == "PSP"]


#akid input
c_akid['kinase_domain'] = c_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c_akid['kinase_domain']= c_akid['kinase_domain'].str.upper()
c_akid['peptide'] = c_akid['peptide'].str.replace("_", '-')

#nwk and psp check in current INKA
KS_relations_nwkin['diff_sub'] = c_akid.peptide.isin(KS_relations.Accession_AApos)      
KS_relations_psp['diff_sub'] = c_akid.peptide.isin(KS_relations.Accession_AApos) 


#merge akid kinase and nwkin kinase
KS_nwkin= pd.merge(KS_relations_nwkin, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
KS_nwkin=KS_nwkin.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)

KS_nwkin["KS"]= KS_nwkin["Accession_AApos"] +"_" +KS_nwkin["akid_kin"]
c_akid["KS"]= c_akid["peptide"] + "_" + c_akid["kinase_domain"]
KS_nwkin['diff_ks'] = c_akid.KS.isin(KS_nwkin.KS)
print(KS_nwkin["diff_ks"].value_counts())
c_akid['diff_ks'] = c_akid.KS.isin(KS_nwkin.KS)
print(c_akid["diff_ks"].value_counts())
c_nwkin_akid = c_akid[c_akid['diff_ks'] == True]#


##PSP merge 
KS_psp= pd.merge(KS_relations_psp, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
KS_psp=KS_psp.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)
KS_psp["KS"]= KS_psp["Accession_AApos"] +"_" +KS_nwkin["akid_kin"]
c_akid["KS"]= c_akid["peptide"] + "_" + c_akid["kinase_domain"]
KS_psp['diff_ks'] = c_akid.KS.isin(KS_nwkin.KS)
print(KS_psp["diff_ks"].value_counts())
c_akid['diff_ks'] = c_akid.KS.isin(KS_psp.KS)
print(c_akid["diff_ks"].value_counts())
c_psp_akid = c_akid[c_akid['diff_ks'] == True]#




###Score of nwkin and psp
nwkin_score= c_nwkin_akid["score"]
psp_score= c_psp_akid["score"]
akid_score= c_akid["score"]


#save them as csv
akid_score.to_csv("akid_gbm_score.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)
nwkin_score.to_csv("nwkin_gbm_score.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)
psp_score.to_csv("psp_gbm_score.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)
#akid_score.to_csv("akid_score.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)