# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:21:47 2022

@author: Tugce Su
"""

##the covarage experiment

##0.99 top10 
import pandas as pd
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/4cell-akid-result.txt", delimiter=r"\s+")
c_akid['kinase_domain'] = c_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c_akid['kinase_domain']= c_akid['kinase_domain'].str.upper()
kinase= c_akid["kinase_domain"].unique()
psp_all= c_akid["peptide"].unique()
c9_akid= c_akid.loc[c_akid["score"] >= 0.99].reset_index()
top10= c9_akid.groupby(by="peptide")["score"].nlargest(10)
top_10=pd.DataFrame(top10)
pep=c9_akid["peptide"].unique()
kin= c9_akid["kinase_domain"].unique()
ind=[]
for i in top_10.index:
    #print(i[1])
    ind.append(i[1])

c10_akid= c9_akid.iloc[ind]
    
