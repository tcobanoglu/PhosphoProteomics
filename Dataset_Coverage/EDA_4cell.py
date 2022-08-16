# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 16:20:26 2022

@author: Tugce Su
"""

import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np 


KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
P_sites = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/new_inka_20220413/Phospho (STY)Sites.txt", sep='\t',low_memory=False)
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins','Sequence window']]
Ph_sites= Ph_sites.dropna()
Ph_sites= Ph_sites.values.tolist()
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
        
##convert df     
ids = pd.DataFrame(list4)
ids.columns = ["id", "res", "pos","seq"]
ids["seq"]= ids["seq"].str[8:23]
ids["res-pos"]= ids["res"] + ids["pos"]
print(ids["res"].value_counts())
##
##covarage
ids = pd.DataFrame(list4)
ids.columns = ["id", "res", "pos","seq"]
#print(ids["res"].value_counts())
ids["Accession_AApos"] = (ids["id"]+ "-"+ids["res"]+ ids["pos"])

ids['diff_sub'] = ids.Accession_AApos.isin(KS_relations.Accession_AApos)      
print(ids["diff_sub"].value_counts())

##graph
pie= ids["diff_sub"].value_counts().plot(kind="pie",labels=['Phosphosite with no kinase','Phosphosite with kinase'],y="Phosphosite",autopct="%1.1f%%",figsize=(16, 8),colors=["c","y"], shadow=True,textprops={'fontsize': 28})
pie.set_title("The Coverage of 4 Cell Lines Dataset",fontsize=36)
###





new_ids= ids[["id","res-pos","seq"]]
new_ids.to_csv("phosph", sep="\t", quoting=csv.QUOTE_NONE,index=False)# all of them
###conver akid convert
##convert for akid format
new_ids["akid"]= ">" + new_ids["id"] + "_"+ new_ids["res-pos"]
new_ids["akid-seq"]= new_ids["akid"] +"\n" +new_ids["seq"]
new_ids= new_ids.drop(["akid","id","res-pos","seq"],axis=1)

#new_ids.to_csv("akid-4cell_lines",index=False,header= 0,sep='|', encoding="utf-8", quoting=csv.QUOTE_NONE,escapechar='|')


























