# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 11:27:41 2022

@author: Tugce Su
"""
import pandas as pd
P_sites = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/new_inka_20220413/Phospho (STY)Sites.txt", sep='\t',low_memory=False)

#SKMEL
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins',"Score SKMel28_2","Score SKMel28_1"]]
Ph_sk1 = Ph_sites[Ph_sites["Score SKMel28_1"].notna()]
Ph_sk2 = Ph_sites[Ph_sites["Score SKMel28_2"].notna()]

##H22
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins',"Score H2228_1","Score H2228_2"]]
Ph_h1 = Ph_sites[Ph_sites["Score H2228_2"].notna()]
Ph_h2 = Ph_sites[Ph_sites["Score H2228_2"].notna()]

##Score K562_1
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins',"Score K562_1","Score K562_2"]]
Ph_k1 = Ph_sites[Ph_sites["Score K562_1"].notna()]
Ph_k2 = Ph_sites[Ph_sites["Score K562_2"].notna()]

##Score HCC827ER3_2
Ph_sites= P_sites[['Proteins','Amino acid','Positions within proteins',"Score HCC827ER3_1","Score HCC827ER3_2"]]
Ph_hc1 = Ph_sites[Ph_sites["Score HCC827ER3_1"].notna()]
Ph_hc2 = Ph_sites[Ph_sites["Score HCC827ER3_2"].notna()]


def psp_combine(Ph1,Ph2,cell_type):
    
    Ph_sites= pd.concat([Ph1,Ph2])
    Ph_sites= Ph_sites.dropna()
    Ph_sites= Ph_sites.values.tolist()
    list1=[] #id
    list2=[] #residue
    list3=[] #position
    list4= []#id-residue-position
    
    
    for row in Ph_sites:
        str(row[0])
        str(row[1])
        str(row[2])
        list1.append((row[0].split(";")))
        list2.append(row[1])
        list3.append(row[2].split(";"))

    for y, row in enumerate(list1):
        for x, unit in enumerate(row):
            list4.append([unit, list2[y],list3[y][x]])
            
    ids = pd.DataFrame(list4)
    ids.columns = ["id", "res", "pos"]
    #print(ids["res"].value_counts())
    ids["Accession_AApos" ] = (ids["id"]+ "-"+ids["res"]+ ids["pos"])
    ids.columns += cell_type
    return(ids)
#####
##
hc= psp_combine(Ph_hc1,Ph_hc2,"_HCC827ER3_2")
sk= psp_combine(Ph_sk1,Ph_sk2,"_SKMEL")
h2= psp_combine(Ph_h1,Ph_h2,"_H22")
k5= psp_combine(Ph_k1,Ph_k2,"_K562_1")           


##AXL Kinase sites: 
    
import pandas as pd
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/4cell-akid-result.txt", delimiter=r"\s+")
c_akid['kinase_domain'] = c_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c_akid['kinase_domain']= c_akid['kinase_domain'].str.upper()
#kinase= c_akid["kinase_domain"].unique()

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

##choose AXL 

axl = c10_akid[c10_akid["kinase_domain"] == "AXL"]
axl['peptide'] = axl['peptide'].str.replace("_", '-')
h2["axl"]= h2.Accession_AApos_H22.isin(axl.peptide)

h2_axl= h2[h2["axl"] == True]
h2_pep= h2_axl["id_H22"].unique()
h2_pep_count= h2_axl["id_H22"].value_counts()

























