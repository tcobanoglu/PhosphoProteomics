# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 10:11:56 2022

@author: Tugce Su
"""

##gbm 

import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np 




KinBase= pd.read_excel("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/Kinase/Table S1.xls")
KinBase=KinBase.rename(columns = {'Name':'Kinase'})
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
P_sites = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/gbm/Phospho (STY)Sites.txt", sep='\t',low_memory=False)
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins','Sequence window']]
Ph_sites= Ph_sites.dropna()
Ph_sites= Ph_sites.values.tolist()
list1=[] #id
list2=[] #residue
list3=[] #position
list4= []#id-residue-position- sequence

for row in Ph_sites:
    list1.append((row[0].split(";")))
    list2.append(row[1])
    list3.append(row[2].split(";"))

for y, row in enumerate(list1):
    for x, unit in enumerate(row):
        list4.append([unit, list2[y],list3[y][x]])
        
        
###merge columns and count###  
#covarage of subs
      
ids = pd.DataFrame(list4)
ids.columns = ["id", "res", "pos"]
#print(ids["res"].value_counts())
ids["Accession_AApos"] = (ids["id"]+ "-"+ids["res"]+ ids["pos"]).unique()
print(ids["res"].value_counts())
ids['diff_sub'] = ids.Accession_AApos.isin(KS_relations.Accession_AApos)      
print(ids["diff_sub"].value_counts())
#pie= ids["diff_sub"].value_counts().plot(kind="pie",labels=['Phosphosite with unknown kinase','Phosphosite with known kinase'],y="Phosphosite",autopct="%1.1f%%",figsize=(16, 8),colors=["c","y"], shadow=True,legend = True)
#pie.set_title("The Covarage of GBM Dataset")
aa_unique= ids["Accession_AApos"].unique()



ids_true = ids[ids['diff_sub'] == True]#matching rows
ids_false = ids[ids['diff_sub'] == False]#unmatches
count_unknown_phos= ids_false["Accession_AApos"].unique()##unique count of unknown psp
count_known_phos= ids_true["Accession_AApos"].unique()
count_all= ids["Accession_AApos"].unique()
print(ids_true["res"].value_counts())
print(ids_false["res"].value_counts())





kin_unmatch= pd.merge(ids_false,KS_relations, on="Accession_AApos",how="right")
kin_unmatch.dropna()

kin_match= pd.merge(ids_true,KS_relations, on="Accession_AApos",how="right")
kin_match=kin_match.dropna()
#kinase=(kin_match["Kinase"].value_counts())#count of match kinases
kinase = kin_match["Kinase"].value_counts().rename_axis('Kinase').reset_index(name='counts')
source_match= kin_match["Source"].value_counts().rename_axis('Source').reset_index(name='counts')#source of match kinases


match= kin_match["Kinase"].unique()
unmatch=kin_unmatch["Kinase"].unique()
print(unmatch)

##D0 preparation
list1=[] #id
list2=[] #residue
list3=[] #position
list4= []#id-residue-position- sequence
list5=[]# sequence

for row in Ph_sites:
    list1.append((row[0].split(";")))
    list2.append(row[1])
    list3.append(row[2].split(";"))
    list5.append(row[3])

for y, row in enumerate(list1):
    for x, unit in enumerate(row):
        list4.append([unit, list2[y],list3[y][x],list5[y]])
        
##convert df for d0  
ids = pd.DataFrame(list4)
ids.columns = ["id", "res", "pos","seq"]
ids["seq"]= ids["seq"].str[8:23]

ids["res-pos"]= ids["res"] + ids["pos"]
new_ids= ids[["id","res-pos","seq"]]
new_ids.to_csv("gbm-psp", sep="\t", quoting=csv.QUOTE_NONE,index=False)# all of them


##convert for akid format
new_ids["akid"]= ">" + new_ids["id"] + "_"+ new_ids["res-pos"]
new_ids["akid-seq"]= new_ids["akid"] +"\n" +new_ids["seq"]
new_ids= new_ids.drop(["akid","id","res-pos","seq"],axis=1)

new_ids.to_csv("akid-gbm",index=False,header= 0,sep='|', encoding="utf-8", quoting=csv.QUOTE_NONE,escapechar='|')


















