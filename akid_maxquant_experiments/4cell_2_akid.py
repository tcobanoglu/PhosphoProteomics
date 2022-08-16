# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:14:21 2022

@author: Tugce Su
"""


import pandas as pd
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
KS_relations["Kinase"]= KS_relations["Kinase"].str.upper()
kin_ks=KS_relations["Kinase"].unique()


#KS_relations= KS_relations.loc[KS_relations["Source"] == "PSP"]
KS_relations= KS_relations.loc[KS_relations["Source"] == "NWK"]
#make all kinase upper case all the time
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID2/celllines2_predictions", delimiter=r"\s+")
c_akid['kinase_domain'] = c_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c_akid['kinase_domain']= c_akid['kinase_domain'].str.upper()
kinase= c_akid["kinase_domain"].unique()
l_kin= set(kinase.tolist()) #akid kinase
psp= c_akid["peptide"].unique()
akid_kin= pd.DataFrame(kinase,columns=["akid_kin"])

c_akid['peptide'] = c_akid['peptide'].str.replace("_", '-')
KS_relations['diff_sub'] = c_akid.peptide.isin(KS_relations.Accession_AApos)      
print(KS_relations["diff_sub"].value_counts())


#open new human kin file
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
#merge akid kinase and hgnc 2016
KS_relations_akid= pd.merge(KS_relations, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
print(KS_relations_akid.columns)
KS_relations_akid=KS_relations_akid.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)


#ks check 
KS_relations_akid["KS"]= KS_relations_akid["Accession_AApos"] +"_" +KS_relations_akid["akid_kin"]
c_akid["KS"]= c_akid["peptide"] + "_" + c_akid["kinase_domain"]
KS_relations_akid['diff_ks'] = c_akid.KS.isin(KS_relations_akid.KS)
print(KS_relations_akid["diff_ks"].value_counts())
c_akid['diff_ks'] = c_akid.KS.isin(KS_relations_akid.KS)
print(c_akid["diff_ks"].value_counts())
