# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:16:28 2022

@author: Tugce Su
"""
import pandas as pd
import csv
akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/AKID-Ensembl.csv", sep=";")
kinase= akid["kinase"].unique()
substrate= akid["substrate_site"]
akid[['substrate_site', 'location']] = akid['substrate_site'].str.split('_', 1, expand=True)
subs= akid["substrate_site"].unique()
name= akid["kinase_name"].unique()
kin_uni= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/kinase_akid_uniprot.txt",sep= "\t")
subs_uni= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/subs_akid_uniprot.txt",sep= "\t")

print(akid.columns)

##merge ids
akid_kin= pd.merge(akid,kin_uni,right_on='From',left_on="kinase",how="left")
akid_kin = akid_kin.drop(columns=['From'])
akid_kin=akid_kin.rename(columns={"To": "Kinase_Uniprot"})
akid_sub= pd.merge(akid_kin,subs_uni,right_on='From',left_on="substrate_site",how="left")
akid_sub = akid_sub.drop(columns=['From'])
akid_uni= akid_sub.rename(columns={"To": "Subs_Uniprot"})
akid_uni.to_csv("akid_id.txt",sep='\t',index=False)
uni= akid_uni["Kinase_Uniprot"].unique()
print(akid_uni.isnull().sum(axis = 0))


no_ids=[]
#no_ids= akid_uni["substrate_site"].iloc[:, akid_uni["Subs_Uniprot"].isna().any()]
#null_data = akid_uni[akid_uni.isnull().any(axis=1)]

subs=akid_uni[akid_uni['Subs_Uniprot'].isnull()]
subs_ensembl=subs["substrate_site"].unique()
subs_ensembl=subs_ensembl.tolist()
kin=akid_uni[akid_uni['Kinase_Uniprot'].isnull()]
kin_ensembl=kin["kinase"].unique()
kin_HGNC= kin["kinase_name"].unique()
kin_ensembl=kin_ensembl.tolist()
no_ids= subs_ensembl + kin_ensembl
with open('ensembl_no_match', 'w') as f:
    write = csv.writer(f)
      
    write.writerow(no_ids)
    
#match ensembl no macth with HGNC uniprot

ids= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/no_symbol.txt", sep='\t')   
print(ids.columns)
ids= ids.drop(columns=["Approved name","Chromosome","NCBI Gene ID"])   

kin_merge= pd.merge(kin,ids,left_on="kinase_name",right_on='Approved symbol',how="left")    
subs_merge=pd.merge(subs,ids,left_on="substrate_name",right_on='Approved symbol',how="left")     
 
##save merges as csv
print(kin_merge.columns)
kin_merge=kin_merge.drop(columns=['substrate_site', 'substrate_name',
       'AKID_score', 'Sequence', 'location', 'Kinase_Uniprot', 'Subs_Uniprot',
       'Approved symbol', 'Ensembl gene ID'])

print(subs_merge.columns)
subs_merge=subs_merge.drop(columns=['kinase', 'kinase_name', 
       'AKID_score', 'Sequence', 'location', 'Kinase_Uniprot', 'Subs_Uniprot',
       'Approved symbol', 'Ensembl gene ID'])

###

subs_merge.to_csv("akid-hgnc-subs.txt",index=False,sep='\t')
kin_merge.to_csv("akid-hgnc-kin.txt",index=False,sep='\t')

#####read the files 
akid_uni1 = open('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid_id.txt', 'r')
akid_lines = akid_uni1.readlines()
akid_list = [line.split("\t") for line in akid_lines]

subs = open('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid-hgnc-subs.txt', 'r')
subs_lines = subs.readlines()
subs_list = [line.split("\t") for line in subs_lines]

kin = open('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid-hgnc-kin.txt', 'r')
kin_lines = kin.readlines()
kin_list = [line.split("\t") for line in kin_lines]
#####

for i in range(len(akid_list)):
    #print(i)
    for y in range(len(kin_list)):
        if len(akid_list[i][7]) < 4 :
            #print ("a")
            if kin_list[y][1] == akid_list[i][2]:
                akid_list[i][7] = kin_list[y][2]
                #print(akid_list[i])
                
                
for i in range(len(akid_list)):
    for x in range(len(subs_list)):
        if len(akid_list[i][8]) < 4:
            #print ("a")
            if subs_list[x][1] == akid_list[i][3]:
                akid_list[i][8] = subs_list[x][2]
               
    #print(akid_list[i])            

df = pd.DataFrame(akid_list)
df1 = df.transpose()

df.to_csv("akid_uni_hgnc.txt",sep='\t',index=False)
#print(df.isnull().sum(axis = 0))    
df= df.replace(" ", "nan")
count_kin=0
count_sub=0
for i in range(len(akid_list)):

   
    if len(akid_list[i][7]) < 3 :
        #print(akid_list[i][7])
        count_kin= count_kin +1
    if len(akid_list[i][8]) < 3 :
        count_sub= count_sub +1
        
"""
##match HGNC 2016
kinase_ids= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/kinase_ids.txt",sep= "\t")
print(kinase_ids.columns)
HGNC_kin = kinase_ids.drop(columns=['Protein.Name.HGNC21Apr2016',
       'Chrom.Location.HGNC21Apr2016', 'DeepKinZero', 'HGNC 2021',
       'Approved name', 'Chromosome','Ensembl gene ID'])
HGNC_akid= pd.merge(HGNC_kin, akid_uni,left_on='uniprot',right_on= "Kinase_Uniprot",how="left")
###
"""

