# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 16:47:29 2022

@author: Tugce Su
"""
import pandas as pd

##gbm result threshold

gbm_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/gbm-akid-result.txt", delimiter=r"\s+")
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
gbm_akid['kinase_domain'] = gbm_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
gbm_akid['kinase_domain']= gbm_akid['kinase_domain'].str.upper()
KS_relations_new= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/new_inka_20220413/KB_new/nwk_mapped_20220411.tsv", sep= "\t",low_memory=False)
psp= gbm_akid["peptide"].unique()



#arrange gbm psp
gbm_akid["peptide_1"]= gbm_akid["peptide"]
gbm_akid['peptide'] = gbm_akid['peptide'].str.replace("_", '-')
KS_relations['diff_sub'] = gbm_akid.peptide.isin(KS_relations.Accession_AApos)      
print(KS_relations["diff_sub"].value_counts())

#open new human kin file
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
#merge akid kinase and hgnc 2016
KS_relations_akid= pd.merge(KS_relations, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
print(KS_relations_akid.columns)
KS_relations_akid=KS_relations_akid.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)
##merge new KS relation
KS_relations_new_akid=  pd.merge(KS_relations_new, kin, left_on= "Kinase.HGNC",right_on="HGNC 2021",how="inner")
print(KS_relations_new_akid.columns)
KS_relations_new_akid=KS_relations_new_akid.drop(['Kinase', 'Kinase.STRING.ID', 'NetworKIN.score',
       'STRING.score', 'Kinase.HGNC', 'HGNC 2016', 'HGNC 2021', 'kinase(dk0)',
       'Approved symbol'],axis=1)

#ks check 
KS_relations_akid["KS"]= KS_relations_akid["Accession_AApos"] +"_" +KS_relations_akid["akid_kin"]
gbm_akid["KS"]= gbm_akid["peptide"] + "_" + gbm_akid["kinase_domain"]
KS_relations_akid['diff_ks'] = gbm_akid.KS.isin(KS_relations_akid.KS)
print(KS_relations_akid["diff_ks"].value_counts())
gbm_akid['diff_ks'] = gbm_akid.KS.isin(KS_relations_akid.KS)
print(gbm_akid["diff_ks"].value_counts())

##new ks
KS_relations_new_akid["KS"]= KS_relations_new_akid["Accession.AApos"] +"_" +KS_relations_new_akid["akid_kin"]
gbm_akid["KS"]= gbm_akid["peptide"] + "_" + gbm_akid["kinase_domain"]
KS_relations_new_akid['diff_ks'] = gbm_akid.KS.isin(KS_relations_new_akid.KS)
print(KS_relations_new_akid["diff_ks"].value_counts())
gbm_akid['diff_ks_new'] = gbm_akid.KS.isin(KS_relations_new_akid.KS)
print(gbm_akid["diff_ks_new"].value_counts())



##gbm true
gbm_true = gbm_akid[gbm_akid['diff_ks'] == True]# same KS with KS relations
psp=gbm_true["peptide"].unique()


##############check threshold
#low 0.25 medium 0.6 high 0.8

##1
psp_all= gbm_akid["peptide"].unique()
gbm_1= gbm_akid.loc[gbm_akid["score"] >= 1]
psp_1= gbm_1["peptide"].unique()
gbm1_true = gbm_1[gbm_1['diff_ks'] == True]
kin_1= gbm_1["kinase_domain"].unique()

##0.9
gbm_9= gbm_akid.loc[gbm_akid["score"] >= 0.9]
psp_9= gbm_9["peptide"].unique()
gbm9_true = gbm_9[gbm_9['diff_ks'] == True]
kin_9= gbm_9["kinase_domain"].unique()

##0.8 
gbm_08= gbm_akid.loc[gbm_akid["score"] >= 0.8]
psp_08= gbm_08["peptide"].unique()
gbm8_true = gbm_08[gbm_08['diff_ks'] == True]
kin_08= gbm_08["kinase_domain"].unique()

##0.75
gbm_75= gbm_akid.loc[gbm_akid["score"] >= 0.75]
psp_75= gbm_75["peptide"].unique()
gbm75_true = gbm_75[gbm_75['diff_ks'] == True]

##0.7
gbm_07= gbm_akid.loc[gbm_akid["score"] >= 0.7]
psp_07= gbm_07["peptide"].unique()
gbm7_true = gbm_07[gbm_07['diff_ks'] == True]

##0.6
gbm_06= gbm_akid.loc[gbm_akid["score"] >= 0.6]
psp_06= gbm_06["peptide"].unique()
gbm6_true = gbm_06[gbm_06['diff_ks'] == True]


## group kinase as Y kinase S/T kinase

#gbm_08["psp"]= gbm_08["peptide_1"].str.split("_")
#Y Kinase
Y_kinase= gbm_08[gbm_08["peptide_1"].str.contains('_Y')]
ST_kinase= gbm_08[gbm_08["peptide_1"].str.contains('_S',"_T")]










