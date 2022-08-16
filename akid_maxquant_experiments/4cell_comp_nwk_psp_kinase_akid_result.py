# -*- coding: utf-8 -*-
"""
Created on Tue May 17 11:18:02 2022

@author: Tugce Su
"""

import pandas as pd
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
KS_relations["Kinase"]= KS_relations["Kinase"].str.upper()
kin_ks=KS_relations["Kinase"].unique()
#KS_relations= KS_relations.loc[KS_relations["Source"] == "PSP"]
KS_relations= KS_relations.loc[KS_relations["Source"] == "PSP"]

#KS_relations1= KS_relations.loc[KS_relations["Kinase"] == "FGFR4"]
#make all kinase upper case all the time
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/4cell-akid-result.txt", delimiter=r"\s+")
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
#seperate ks relations dataset based on the source
KS_relations_akid["KS"]= KS_relations_akid["Accession_AApos"] +"_" +KS_relations_akid["akid_kin"]
KS_psp= KS_relations_akid.loc[KS_relations["Source"] == "PSP"]
KS_nwk= KS_relations_akid.loc[KS_relations["Source"] == "NWK"]


#compare akid result
c_akid["KS"]= c_akid["peptide"] + "_" + c_akid["kinase_domain"]
#c_akid['diff_ks'] = c_akid.KS.isin(KS_relations_akid.KS)
c_akid['diff_ks'] = c_akid.KS.isin(KS_psp.KS)
#c_akid['psp'] = c_akid.peptide.isin(KS_psp.Accession_AApos)
#print(c_akid["psp"].value_counts())
print(c_akid["diff_ks"].value_counts())
#change score

c_8= c_akid.loc[c_akid["score"] >= 0.8]
c_9= c_akid.loc[c_akid["score"] >= 0.9]
c_95= c_akid.loc[c_akid["score"]>= 0.95]
c_99= c_akid.loc[c_akid["score"]>= 0.99]
c_999= c_akid.loc[c_akid["score"]>= 0.999]
c_1= c_akid.loc[c_akid["score"]>= 1.0]

##for kinase
c_8= c_akid.loc[c_akid["score"] >=  0.8]
kin8=c_8.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'

c_9= c_akid.loc[c_akid["score"] >=  0.9]
kin9=c_9.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'

c_95= c_akid.loc[c_akid["score"] >=  0.95]
kin95=c_95.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'


c_99= c_akid.loc[c_akid["score"] >=  0.99]
kin99=c_99.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'


c_999= c_akid.loc[c_akid["score"] >=  0.999]
kin999=c_999.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'

c_1= c_akid.loc[c_akid["score"] >=  1.0]
kinn=c_1.kinase_domain.value_counts(normalize=True).mul(100).round(1).astype(str) + '%'


#KS_relations_akid['diff_ks'] = c_akid.KS.isin(KS_relations_akid.KS)
#print(KS_relations_akid["diff_ks"].value_counts())
#c_akid['diff_ks'] = c_akid.KS.isin(KS_relations_akid.KS)
print(c_akid["diff_ks"].value_counts())


##gbm true
"""
c_true = c_akid[c_akid['diff_ks'] == True]# same KS with KS relations
c_08= c_akid.loc[c_akid["score"] >= 0.8]#number of psp
c_0_psp= c_08["peptide"].unique()
c8_true = c_08[c_08['diff_ks'] == True]#match with KS 
kin_08= c_08["kinase_domain"].unique()
"""
##PSP covarage

c8_psp= c_8["peptide"].unique()
c8_true = c_8[c_8['diff_ks'] == True]   
psp_8 = len(c8_true)/len(c8_psp)*100
print(len(c8_psp)/len(psp)*100)

c9_psp= c_9["peptide"].unique()
c9_true = c_9[c_9['diff_ks'] == True]   
psp_9 = len(c9_true)/len(c9_psp)*100
print(len(c9_psp)/len(psp)*100)

c95_psp= c_95["peptide"].unique()
c95_true = c_95[c_95['diff_ks'] == True]   
psp_95 = len(c95_true)/len(c95_psp)*100
print(len(c95_psp)/len(psp)*100)


c99_psp= c_99["peptide"].unique()
#print(c_99["diff_ks"].value_counts())

c99_true = c_99[c_99['diff_ks'] == True]   
psp_99 = len(c99_true)/len(c99_psp)*100
print(len(c99_psp)/len(psp)*100)



c999_psp= c_999["peptide"].unique()
c999_true = c_999[c_999['diff_ks'] == True]   
psp_999 = len(c999_true)/len(c999_psp)*100
print(len(c999_psp)/len(psp)*100)


c1_psp= c_1["peptide"].unique()
c1_true = c_1[c_1['diff_ks'] == True]   
psp_1 = len(c1_true)/len(c1_psp)*100
print(len(c1_psp)/len(psp)*100)