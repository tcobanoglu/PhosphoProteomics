# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 10:49:04 2021

@author: Tugce Su
"""

import csv
import pandas as pd
#P_sites = "C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/INKA_urea/Phospho (STY)Sites.txt"
import matplotlib.pyplot as plt
import numpy as np 
#plt.style.use('classic')
####open the files####

P_sites = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/INKA_urea/Phospho (STY)Sites.txt", sep='\t',low_memory=False)
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
#print(type(KS_relations["Accession_AApos"]))
KinBase= pd.read_excel("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/Kinase/Table S1.xls")
KinBase=KinBase.rename(columns = {'Name':'Kinase'})

"""
######################


evidence = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/INKA_urea/evidence.txt",sep="\t",low_memory=False)
e_name = evidence.columns
P_name = P_sites.columns
M_peptides = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/INKA_urea/modificationSpecificPeptides.txt", sep = "\t",low_memory=False)
M_name = M_peptides.columns
"""
####Columns##########
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins']]
Ph_sites= Ph_sites.dropna()
Ph_sites= Ph_sites.values.tolist()
list1=[] #id
list2=[] #residue
list3=[] #position
list4= []#id-residue-position

for row in Ph_sites:
    list1.append((row[0].split(";")))
    list2.append(row[1])
    list3.append(row[2].split(";"))

for y, row in enumerate(list1):
    for x, unit in enumerate(row):
        list4.append([unit, list2[y],list3[y][x]])
############################


        
###merge columns and count###  
#covarage of subs
      
ids = pd.DataFrame(list4)
ids.columns = ["id", "res", "pos"]
#print(ids["res"].value_counts())
ids["Accession_AApos"] = (ids["id"]+ "-"+ids["res"]+ ids["pos"]).unique()
ids['diff_sub'] = ids.Accession_AApos.isin(KS_relations.Accession_AApos)      
print(ids["diff_sub"].value_counts())
##graph
"""

pie= ids["diff_sub"].value_counts().plot(kind="pie",labels=['Phosphosite with unknown kinase','Phosphosite with known kinase'],y="Phosphosite",autopct="%1.1f%%",figsize=(16, 8),colors=["c","y"], shadow=True,legend = True)
pie.set_title("The Covarage of Urea Dataset")
"""
####
ids_true = ids[ids['diff_sub'] == True]#matching rows
ids_false = ids[ids['diff_sub'] == False]#unmatches
count_unknown_phos= ids_false["Accession_AApos"].unique()##unique count of unknown psp
count_known_phos= ids_true["Accession_AApos"].unique()
count_all= ids["Accession_AApos"].unique()
print(ids_true["res"].value_counts())
print(ids_false["res"].value_counts())
"""
##CophosK
ids["subs"] = (ids["id"]+ " "+ids["res"]+ ids["pos"]).unique()
CoPhosK = ids.to_csv("CoPhosK_nomatch", sep="\t", columns=["subs"],index=False)
AKID= ids.to_csv("AKID",sep="\t", columns=["id"],index=False)
ids_false_ = ids_false.to_csv("unknown_phos",sep="\t", columns=["id"],index=False)

ids_true_ = ids_true.to_csv("known_phos",sep="\t", columns=["id"],index=False)
"""
##graph
"""
res_pie=ids_true["res"].value_counts().plot(kind="pie",autopct="%1.1f%%",figsize=(16, 8),colors=["c", "g", "y"], shadow=True,legend=True)
res_pie.set_title("Residues with known KS association")
"""
"""
res_pie1=ids_false["res"].value_counts().plot(kind="pie",autopct="%.2f",figsize=(16, 8),colors=["c", "g", "y"],shadow=True,legend=True)
res_pie1.set_title("Residues with unknown KS association")
"""
######
#####################################



######Find unmatch and match kinases 

kin_unmatch= pd.merge(ids_false,KS_relations, on="Accession_AApos",how="right")
kin_unmatch.dropna()

kin_match= pd.merge(ids_true,KS_relations, on="Accession_AApos",how="right")
kin_match=kin_match.dropna()
#kinase=(kin_match["Kinase"].value_counts())#count of match kinases
kinase = kin_match["Kinase"].value_counts().rename_axis('Kinase').reset_index(name='counts')
source_match= kin_match["Source"].value_counts().rename_axis('Source').reset_index(name='counts')#source of match kinases
#kinase= pd.DataFrame(kinase)
#kinase.columns= ["Kinase"]

#graph
"""
graph= kin_match["Source"].value_counts()
#fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
s_graph=graph.plot.bar(
               #ylim=(0,313),
              
               title='The Source of Reported Phosphosites',
               #legend=False,
               color="cyan",
               #alpha=0.5,
               rot=0,figsize=(12,6),edgecolor="black",xlabel="Database",ylabel="Counts")
plt.xticks(size=12)
plt.yticks(size=12)
plt.xlabel(xlabel="Database",size=16)
plt.ylabel(ylabel="Counts",size=16)
plt.title(label='The Source of Reported Phosphosites',size=20)
s_graph.set_facecolor('white')
s_graph.set_axis_bgcolor('white')
s_graph.patch.set_facecolor('xkcd:white')
"""


KS_relations['kin_check'] = KS_relations.Kinase.isin(kinase.Kinase)
KS_unmatch = KS_relations[KS_relations['kin_check'] == False]#unmatches
#no_sub=KS_unmatch["Kinase"].value_counts()#unmatch kinases
no_sub = KS_unmatch["Kinase"].value_counts().rename_axis('Kinase').reset_index(name='counts')
no_sub= pd.DataFrame(no_sub)
source_unmatch= KS_unmatch["Source"].value_counts().rename_axis('Source').reset_index(name='counts')

##########################################





#####KİnBase Comparision####Group/Superfamily


##kinase with no subs

no_sub['Kinase']=no_sub['Kinase'].astype(str)
KS_missing= pd.merge(no_sub,KinBase, on="Kinase",how="inner")
#KS_missing=no_sub.merge(KinBase, left_on='Kinase', right_on='Kinase')
KS_no_group=(KS_missing.groupby(["Group"])["counts"].sum())
KS_no_group=KS_no_group.rename_axis('Group').reset_index(name='counts')

























##kinase with subs

KS_match_group= pd.merge(kinase,KinBase, on="Kinase",how="inner")
KS_group= KS_match_group.groupby(["Group"])["counts"].sum()
KS_group=KS_group.rename_axis('Group').reset_index(name='counts')
#.plot.barh(y='Count',color="c",stacked=True,figsize=(8, 6))




#graph

main_group= pd.merge(KS_group,KS_no_group, on="Group",how="inner").rename(columns={"Group":"Group","counts_x":"Reported phosphosites","counts_y":"Unreported phosphosites"})
main_group.reset_index()
main_group['Reported phosphosites']=main_group['Reported phosphosites'].astype(int)
main_group['Unreported phosphosites']=main_group['Unreported phosphosites'].astype(int)
index=["AGC","Atypical","CAMK","CK1","CMGC","Other","STE","TK","TKL"]
#main_group.plot(kind="bar",labels=index)


a=main_group["Reported phosphosites"]

print(a)
b=main_group["Unreported phosphosites"].astype(int)
df = pd.DataFrame({"Reported Phosphosites":list(KS_group["counts"]),
                   'Unreported phosphosites':list(KS_no_group["counts"]) }, index=index)

ax = df.plot.bar(rot=0,color={"Reported Phosphosites":"cyan", "Unreported phosphosites": "purple"},figsize=(12,10))

plt.title(label='Reported and Unreported Kinase Groups with NetworKIN and PhosphoSitePlus',size=20)
plt.xlabel('Kinase Groups',size=14)
plt.ylabel('Counts',size=14)
plt.xticks(size=12)
plt.yticks(size=12)
ax.set_axisbelow(True)
ax.yaxis.grid(True, color='#EEEEEE')
ax.xaxis.grid(False)
ax.set_facecolor('white')
#fig.tight_layout()

#######################






































            
        
       