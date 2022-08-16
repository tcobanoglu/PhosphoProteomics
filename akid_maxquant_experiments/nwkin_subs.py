# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:03:12 2022

@author: Tugce Su
"""

##discover NWKIN substrate

import pandas as pd
from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt
from scipy import stats


nwk_subs= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID2/NWK_substrates_4cell.txt", sep= "\t")# 0.999 top5
akid_subs= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA_code/inka-4cell/InKA_20220330/NWK_substrates_1.1_nwkin_threshold_akidonly.txt", sep= "\t")#0.999 top10
#print(nwk_subs.columns)


def file_each_psp(col,inp): 
    usecols= [col,'Accession_AApos']
    file= pd.read_csv(inp, sep= "\t",usecols=usecols)

    file= file[file[col] != 0]
    count= file[col].sum()
    print(count)
    return(file)
akid= "C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA_code/inka-4cell/InKA_20220330/NWK_substrates_1.1_nwkin_threshold_akidonly.txt"#0.999 top10
nwkin= "C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID2/NWK_substrates_4cell.txt"# 0.999 top5

h22_akid= file_each_psp('Spectral.Count.H2228_1',akid)
hcc_akid= file_each_psp('Spectral.Count.HCC827ER3_1',akid)
k5_akid= file_each_psp('Spectral.Count.K562_1',akid)
sk_akid= file_each_psp('Spectral.Count.SKMel28_1',akid)

h22_nwkin= file_each_psp('Spectral.Count.H2228_1',nwkin)
hcc_nwkin= file_each_psp('Spectral.Count.HCC827ER3_1',nwkin)
k5_nwkin= file_each_psp('Spectral.Count.K562_1',nwkin)
sk_nwkin= file_each_psp('Spectral.Count.SKMel28_1',nwkin)
##AKID
h22_la= h22_akid['Accession_AApos'].tolist()
hcc_la= hcc_akid['Accession_AApos'].tolist()
k5_la= k5_akid['Accession_AApos'].tolist()
sk_la= sk_akid['Accession_AApos'].tolist()





plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of HCC827ER3 and H2228", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_la), set(h22_la)],set_labels = ('HCC827ER3', 'H2228'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_hcc_h22.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of HCC827ER3 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_la), set(sk_la)],set_labels = ('HCC827ER3', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_hcc_sk.pdf', format='pdf')   




plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of HCC827ER3 and K562", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_la), set(k5_la)],set_labels = ('HCC827ER3', 'K562'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_hcc_k5.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of H2228 and K562", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(h22_la), set(k5_la)],set_labels = ('H2228', 'K562'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_k5_h22.pdf', format='pdf')   


plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of H2228 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(h22_la), set(sk_la)],set_labels = ('H2228', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_sk_h22.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("AKID Psp of K562 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(k5_la), set(sk_la)],set_labels = ('K562', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/AKID_k5_sk.pdf', format='pdf')   

##NWKIN
h22_ln= h22_nwkin['Accession_AApos'].tolist()
hcc_ln= hcc_nwkin['Accession_AApos'].tolist()
k5_ln= k5_nwkin['Accession_AApos'].tolist()
sk_ln= sk_nwkin['Accession_AApos'].tolist()


plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of HCC827ER3 and H2228", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_ln), set(h22_ln)],set_labels = ('HCC827ER3', 'H2228'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_hc_h2.pdf', format='pdf')   




plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of HCC827ER3 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_ln), set(sk_ln)],set_labels = ('HCC827ER3', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_hc_sk.pdf', format='pdf')   




plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of HCC827ER3 and K562", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(hcc_ln), set(k5_ln)],set_labels = ('HCC827ER3', 'K562'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_hc_k5.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of H2228 and K562", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(h22_ln), set(k5_ln)],set_labels = ('H2228', 'K562'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_k5_h2.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of H2228 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(h22_ln), set(sk_ln)],set_labels = ('H2228', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_sk_h2.pdf', format='pdf')   



plt.figure(figsize = (6, 5))
font1 = {'family':'serif','color':'black','size':20} 
plt.title("NWKIN Psp of K562 and SKMEL28", fontdict=font1)
font2 = {'family': 'Comic Sans MS', 'size': 18} # use for labels
plt.rc('font', **font2) # sets the default font 
plt.rcParams['text.color'] = 'black' # changes default text colour
venn2([set(k5_ln), set(sk_ln)],set_labels = ('K562', 'SKMEL28'),set_colors=("cyan", "pink"))
#plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/nwk_k5_sk.pdf', format='pdf')   



def jaccard_sim(a, b):
      a= set(a)
      b= set(b)
      c = a.intersection(b)
     # print(float(len(c)) / (len(a) + len(b) - len(c)))
      return float(len(c)) / (len(a) + len(b) - len(c))

hc_h2=  jaccard_sim(hcc_la,h22_la)
sk_h2=  jaccard_sim(sk_la,h22_la)
k5_h2=  jaccard_sim(k5_la,h22_la)
sk_hc=  jaccard_sim(sk_la,hcc_la)
sk_k5=  jaccard_sim(sk_la,k5_la)
hc_k5=  jaccard_sim(hcc_la,k5_la)

print(hc_h2)
print(sk_h2)
print(k5_h2)
print(sk_hc)
print(sk_k5)
print(hc_k5)


import seaborn as sns

#sns.clustermap(akid_subs)
##one hot encoding

kin_dummy= pd.get_dummies(akid_subs["Kinase.HGNC"])


akid_10 = pd.merge(
    left=akid_subs,
    right=kin_dummy,
    left_index=True,
    right_index=True,
)

#categorical encoding/ label
#from sunbird.categorical_encoding import frequency_encoding
akid_10["Accession_AApos"] = akid_10["Accession_AApos"].astype("category")
akid_10["Accession_AApos"] = akid_10["Accession_AApos"].cat.codes
#akid_10= akid_10.drop(columns=["Kinase.HGNC"])

akid_10["Accession_AApos"] = akid_10["Accession_AApos"].astype("category")
akid_10["Accession_AApos"] = akid_10["Accession_AApos"].cat.codes
sns.clustermap(akid_10)




kin_dummy= pd.get_dummies(akid_subs["Accession_AApos"])


akid_10 = pd.merge(
    left=akid_subs,
    right=kin_dummy,
    left_index=True,
    right_index=True,
)

#categorical encoding/ label
#from sunbird.categorical_encoding import frequency_encoding
import fastcluster
"""
from sklearn.cluster import KMeans
akid_10= akid_subs
akid_10["Kinase.HGNC"] = akid_10["Kinase.HGNC"].astype("category")
akid_10["Kinase.HGNC"] = akid_10["Kinase.HGNC"].cat.codes
akid_10["Accession_AApos"] = akid_10["Accession_AApos"].astype("category")
akid_10["Accession_AApos"] = akid_10["Accession_AApos"].cat.codes
#akid_10= akid_10.drop(columns=["Accession_AApos"])
c1= akid_10.loc[1:3000]
sns.clustermap(c1)


k1= akid_subs[1:300]
kmeans= KMeans(n_clusters=3)
kmeans.fit(akid_10)
centroids = kmeans.cluster_centers_
print(centroids)
plt.scatter(akid_10['Accession_AApos'], akid_10['Kinase.HGNC'], c= kmeans.labels_.astype(float), s=50, alpha=0.5)
plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
plt.show()

""""""

























