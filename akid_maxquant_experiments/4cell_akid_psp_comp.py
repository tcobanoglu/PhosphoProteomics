# -*- coding: utf-8 -*-
"""
Created on Mon May 23 09:55:30 2022

@author: Tugce Su
"""
###get akid psp giher than 0.99
##Venn diagram of shared psp
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
c10_akid= c9_akid.iloc[ind]## has aa_pos to compare
#######

####get psp of each cell line
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


###HCC
hc['diff_sub'] = hc.Accession_AApos_HCC827ER3_2.isin(h2.Accession_AApos_H22)
print(hc['diff_sub'].value_counts())

#after 0.99 threshold


#c10_akid["peptide"] = c10_akid["peptide"].str.replace('_','-')
#hc_akid= pd.merge(c10_akid, hc, how='left', right_on='Accession_AApos_HCC827ER3_2', left_on='peptide')
#psp_count=(hc_akid['diff_sub'].value_counts())
#print(psp_count)
#hc_h2 = hc_akid[hc_akid['diff_sub'] == True]
#pep_comm= hc_h2["peptide"].unique()
#hc_h2 = hc_akid[hc_akid['diff_sub'] == False]

#pep_diff= hc_h2["peptide"].unique()

##H22

h2['diff_sub'] = h2.Accession_AApos_H22.isin(hc.Accession_AApos_HCC827ER3_2)
print(h2['diff_sub'].value_counts())

##make list for venn diagram

l_hc= hc['Accession_AApos_HCC827ER3_2'].tolist()
l_h2= h2["Accession_AApos_H22"].tolist()
l_sk= sk["Accession_AApos_SKMEL"].tolist()
l_k5= k5["Accession_AApos_K562_1"].tolist()

#venn

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles

from matplotlib import pyplot as plt
import venn

##
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of HCC827ER3 and H2228", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(l_hc), set(l_h2)],set_labels = ('HCC827ER3', 'H2228'),set_colors=("cyan", "pink"))
#venn2([(l_hc), (l_h2)],set_labels = ('HCC827ER3', 'H2228'),set_colors=("cyan", "pink"))


###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of SKMEL and K562", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_k5), set(l_sk)],set_labels = ('K562', 'SKMEL'),set_colors=("cyan", "pink"))

###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of 'HCC827ER3' and K562_1", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_hc), set(l_k5)],set_labels = ('HCC827ER3', 'K562'),set_colors=("cyan", "pink"))

####

plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of H2228 and SKMEL", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_h2), set(l_sk)],set_labels = ('H2228', 'SKMEL'),set_colors=("cyan", "pink"))


###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of HCC827ER3 and SKMEL", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_hc), set(l_sk)],set_labels = ('HCC827ER3', 'SKMEL'),set_colors=("cyan", "pink"))



###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of K562 and H2228", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_k5), set(l_h2)],set_labels = ('K562', 'H2228'),set_colors=("cyan", "pink"))

####after threshold

#hc
c10_akid["peptide"] = c10_akid["peptide"].str.replace('_','-')
hc['diff_sub'] = hc.Accession_AApos_HCC827ER3_2.isin(c10_akid.peptide)
print(hc['diff_sub'].value_counts())
hc= hc[hc['diff_sub'] == True]

#sk
c10_akid["peptide"] = c10_akid["peptide"].str.replace('_','-')
sk['diff_sub'] = sk.Accession_AApos_SKMEL.isin(c10_akid.peptide)
print(sk['diff_sub'].value_counts())
sk= sk[sk['diff_sub'] == True]
 

#h2
h2['diff_sub'] = h2.Accession_AApos_H22.isin(c10_akid.peptide)
print(h2['diff_sub'].value_counts())
h2= h2[h2['diff_sub'] == True]

#_K562_1
k5['diff_sub'] = k5.Accession_AApos_K562_1.isin(c10_akid.peptide)
print(k5['diff_sub'].value_counts())
k5= k5[k5['diff_sub'] == True]

##find similarity between datasets(kinase included)

hc_merge= pd.merge(hc, c10_akid, left_on="Accession_AApos_HCC827ER3_2",right_on="peptide", how= "left" )
sk_merge= pd.merge(sk, c10_akid, left_on="Accession_AApos_SKMEL",right_on="peptide", how= "left" )
h2_merge= pd.merge(h2, c10_akid, left_on="Accession_AApos_H22",right_on="peptide", how= "left" )
k5_merge= pd.merge(k5, c10_akid, left_on="Accession_AApos_K562_1",right_on="peptide", how= "left" )
#merge kinase
hc_merge["KS"]= hc_merge["Accession_AApos_HCC827ER3_2"] +"-" + hc_merge["kinase_domain"]
sk_merge["KS"]= sk_merge["Accession_AApos_SKMEL"] +"-" + sk_merge["kinase_domain"]
h2_merge["KS"]= h2_merge["Accession_AApos_H22"] +"-" + h2_merge["kinase_domain"]
k5_merge["KS"]= k5_merge["Accession_AApos_K562_1"] +"-" + k5_merge["kinase_domain"]

#check KS similarity
sk_merge["diff_ks"]= sk_merge.KS.isin(hc_merge.KS)
print(sk_merge["diff_ks"].value_counts(normalize=True))

h2_merge["diff_ks"]= h2_merge.KS.isin(hc_merge.KS)
print(h2_merge["diff_ks"].value_counts(normalize=True))

k5_merge["diff_ks"]= k5_merge.KS.isin(hc_merge.KS)
print(k5_merge["diff_ks"].value_counts(normalize=True))

#after threshold how many of them are each other
#hc["diff_sk"]= hc.Accession_AApos_HCC827ER3_2.isin(sk.Accession_AApos_SKMEL)
#print(hc["diff_sk"].value_counts(normalize=True))
##HCC AKID accession id result similarity
sk["diff_hc"]= sk.Accession_AApos_SKMEL.isin(hc.Accession_AApos_HCC827ER3_2)
print(sk["diff_hc"].value_counts(normalize=True))

k5["diff_hc"]= k5.Accession_AApos_K562_1.isin(hc.Accession_AApos_HCC827ER3_2)
print(k5["diff_hc"].value_counts(normalize=True))

h2["diff_hc"]= h2.Accession_AApos_H22.isin(hc.Accession_AApos_HCC827ER3_2)
print(h2["diff_hc"].value_counts(normalize=True))

"""
venn3(([set(l_hc), set(l_h2),set(l_k5)]), set_labels = ('HCC827ER3_2', 'H22','K562_1'))


venn4(([set(l_hc), set(l_h2),set(l_k5),set(l_sk)]), set_labels = ('HCC827ER3_2', 'H22','K562_1','SKMEL'))


a= venn.venn4(([set(l_hc), set(l_h2),set(l_k5),set(l_sk)]), names = ['HCC827ER3_2', 'H22','K562_1','SKMEL'])
"""


l_hc= hc['Accession_AApos_HCC827ER3_2'].tolist()
l_h2= h2["Accession_AApos_H22"].tolist()
l_sk= sk["Accession_AApos_SKMEL"].tolist()
l_k5= k5["Accession_AApos_K562_1"].tolist()



##
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of HCC827ER3_2 and H22", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(l_hc), set(l_h2)],set_labels = ('HCC827ER3_2', 'H22'),set_colors=("cyan", "pink"))

###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of SKMEL and K562_1", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_k5), set(l_sk)],set_labels = ('K562_1', 'SKMEL'),set_colors=("cyan", "pink"))

###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of 'HCC827ER3_2' and K562_1", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_hc), set(l_k5)],set_labels = ('HCC827ER3_2', 'K562_1'),set_colors=("cyan", "pink"))

####

plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of 'H22' and K562_1", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_h2), set(l_sk)],set_labels = ('H22', 'K562_1'),set_colors=("cyan", "pink"))


###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of 'HCC827ER3_2' and SKMEL", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_hc), set(l_sk)],set_labels = ('HCC827ER3_2', 'SKMEL'),set_colors=("cyan", "pink"))



###
plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("Phosphosites of K562_1 and H22", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour

venn2([set(l_k5), set(l_h2)],set_labels = ('K562_1', 'H22'),set_colors=("cyan", "pink"))

def jaccard_sim(a, b):
      a= set(a)
      b= set(b)
      c = a.intersection(b)
      return float(len(c)) / (len(a) + len(b) - len(c))

hc_h2=  jaccard_sim(l_hc,l_h2)
sk_h2=  jaccard_sim(l_sk,l_h2)
k5_h2=  jaccard_sim(l_k5,l_h2)
sk_hc=  jaccard_sim(l_sk,l_hc)
sk_k5=  jaccard_sim(l_sk,l_k5)
hc_k5=  jaccard_sim(l_hc,l_k5)

print(hc_h2)
print(sk_h2)
print(k5_h2)
print(sk_hc)
print(sk_k5)
print(hc_k5)

##calculate the percentile of score
from scipy import stats



























