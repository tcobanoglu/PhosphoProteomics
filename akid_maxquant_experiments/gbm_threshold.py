# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 11:38:45 2022

@author: Tugce Su
"""

##gbm threshold cut off

import pandas as pd
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
#psp= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID2/whole_proteome_predictions", delimiter=r"\s+")
psp= psp= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/urea-akid-result.txt", delimiter=r"\s+")
psp_99= psp.loc[psp["score"] >= 0.999]
###top10 with 0.9

top10= psp_99.groupby(by="peptide")["score"].nlargest(5)
top_10=pd.DataFrame(top10)
##
ind=[]
#print(top_10.index)
for i in top_10.index:
    #print(i)
    ind.append(i[1])
    

c10_akid= psp.iloc[ind]
c10_akid['Accession_AApos'] = c10_akid['peptide'].str.replace("_", '-')


##kinase match
##HGNC 2016
c10_akid['kinase_domain'] = c10_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c10_akid['kinase_domain']= c10_akid['kinase_domain'].str.upper()

kin["HGNC 2016"]= kin['HGNC 2016'].fillna(kin['akid_kin'])

#merge akid kinase-hgnc2016
psp_akid= pd.merge(c10_akid, kin, left_on= "kinase_domain",right_on="HGNC 2016",how="inner")
print(psp_akid.columns)
psp_akid= psp_akid.drop(['kinase_domain', 'peptide', 'score', 'HGNC 2021', 'kinase(dk0)', 'Approved symbol', 'akid_kin'], axis=1)
psp_akid["Kinase.HGNC"] =psp_akid["HGNC 2016"]
psp_akid= psp_akid.drop(["HGNC 2016"],axis=1)
file = psp_akid.to_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID3/urea_0.999_top5.txt",sep="\t",index=False)

