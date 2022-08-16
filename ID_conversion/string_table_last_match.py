# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:46:08 2022

@author: Tugce Su
"""

import pandas as pd 
akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid_uni_match1.txt", sep="\t")
string= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/gene_translation_table.tsv",sep="\t",usecols= ["string_id","uniprot"])
print(akid.isnull().sum(axis = 0))


df1 = akid[akid.isna().any(axis=1)]
kinase= akid["kinase"][akid['Kinase_Uniprot'].isnull()]
kinase= kinase.unique()
kin= pd.DataFrame(kinase)
kin.columns=["id"]

phosph= akid["substrate_site"][akid['Subs_Uniprot'].isnull()]
phosph=phosph.unique()
psp= pd.DataFrame(phosph)
psp.columns=["id"]



psp_merge= pd.merge(psp,string,right_on="string_id",left_on="id",how="inner")
psp_merge=psp_merge.drop(columns=['string_id'])

kin_merge= pd.merge(kin,string,right_on="string_id",left_on="id",how="inner")
kin_merge= kin_merge.drop(columns=["string_id"])

akid_merge= pd.merge(akid,kin_merge,left_on= "kinase", right_on="id",how="left")
akid_merge['Kinase_Uniprot'].fillna(akid_merge['uniprot'], inplace=True)
akid= akid_merge.drop(columns=['uniprot',"id"])

akid_merge= pd.merge(akid,psp_merge,left_on= "substrate_site", right_on="id",how="left")
akid_merge['Subs_Uniprot'].fillna(akid_merge['uniprot'], inplace=True)
akid= akid_merge.drop(columns=['uniprot',"id"])


print(akid.isnull().sum(axis = 0))

####
#find unmatch ids after string
###

kinase1= akid["kinase"][akid['Kinase_Uniprot'].isnull()]
kinase_string= kinase1.unique()
kinase_string= pd.DataFrame(kinase_string)
kinase_string.columns=["id"]
kinase_string.to_csv("unmatch-ensembl-kin.txt",index=False,sep='\t')

phosph1= akid["substrate_site"][akid['Subs_Uniprot'].isnull()]
phosph_string=phosph1.unique()
print(phosph_string)
phosph_string= pd.DataFrame(phosph_string)
phosph_string.columns=["id"]
phosph_string.to_csv("unmatch-ensembl-psp.txt",index=False,sep='\t')

###for INKA
##
aa= akid["Subs_Uniprot"]+"-"+akid["location"]
kin=akid["kinase_name"]
d = {'Accession_AApos': aa, 'Kinase': kin}
inka = pd.DataFrame(d)

inka.to_csv("akid-inka.txt",index=False,sep='\t')
##
akid_inka= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid-inka.txt", sep="\t")
akid_inka= akid_inka.dropna()

akid_inka.to_csv("akid-inka.txt",index=False,sep='\t')










