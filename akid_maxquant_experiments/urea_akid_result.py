# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 14:00:30 2022

@author: Tugce Su
"""
#akid urea cell lines
import pandas as pd
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
KS_relations["Kinase"]= KS_relations["Kinase"].str.upper()
kin_ks=KS_relations["Kinase"].unique()

#make all kinase upper case all the time
urea_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/urea-akid-result.txt", delimiter=r"\s+")
urea_akid['kinase_domain'] = urea_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
urea_akid['kinase_domain']= urea_akid['kinase_domain'].str.upper()
kinase= urea_akid["kinase_domain"].unique()
l_kin= set(kinase.tolist()) #akid kinase
psp= urea_akid["peptide"].unique()
akid_kin= pd.DataFrame(kinase,columns=["akid_kin"])

urea_akid['peptide'] = urea_akid['peptide'].str.replace("_", '-')
KS_relations['diff_sub'] = urea_akid.peptide.isin(KS_relations.Accession_AApos)      
print(KS_relations["diff_sub"].value_counts())


#open new human kin file
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
#merge akid kinase and hgnc 2016
KS_relations_akid= pd.merge(KS_relations, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
print(KS_relations_akid.columns)
KS_relations_akid=KS_relations_akid.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)

#ks check 
KS_relations_akid["KS"]= KS_relations_akid["Accession_AApos"] +"_" +KS_relations_akid["akid_kin"]
urea_akid["KS"]= urea_akid["peptide"] + "_" + urea_akid["kinase_domain"]
KS_relations_akid['diff_ks'] = urea_akid.KS.isin(KS_relations_akid.KS)
print(KS_relations_akid["diff_ks"].value_counts())
urea_akid['diff_ks'] = urea_akid.KS.isin(KS_relations_akid.KS)
print(urea_akid["diff_ks"].value_counts())


u= urea_akid["score"].describe()
##gbm true
urea_true = urea_akid[urea_akid['diff_ks'] == True]# same KS with KS relations
#psp=urea_true["peptide"].unique()


urea_group= urea_akid.groupby(["peptide"])["kinase_domain"].count()

###urea threshold
urea_08= urea_akid.loc[urea_akid["score"] >= 0.8]#number of psp
urea_0_psp= urea_08["peptide"].unique()
cov= len(urea_0_psp)/len(psp) *100
urea8_true = urea_08[urea_08['diff_ks'] == True]#match with KS 
kin_08= urea_08["kinase_domain"].unique()


##covarage of 1/3 gbm
import matplotlib.pyplot as plt
import numpy as np

y=np.array([8.48,91.52])
mylabels = ["Identical KS", "Different KS"]
mycolors = ["c","y"]
plt.pie(y, labels = mylabels, colors = mycolors,autopct='%1.1f%%',shadow=True, startangle=90)
plt.title("The Covarage of Annotated KS in GBM Dataset", pad=32)
plt.show() 



from matplotlib.ticker import StrMethodFormatter
ax = urea_true.hist(column='score', bins=25, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

ax = ax[0]
for x in ax:

    # Despine
    x.spines['right'].set_visible(False)
    x.spines['top'].set_visible(False)
    x.spines['left'].set_visible(False)

    # Switch off ticks
    x.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on")

    # Draw horizontal axis lines
    vals = x.get_yticks()
    for tick in vals:
        x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)

    # Remove title
    x.set_title("The Histogram of Common KS Associations",size=20)

    # Set x-axis label
    x.set_xlabel("AKID Score", labelpad=20, weight='bold', size=16)

    # Set y-axis label
    x.set_ylabel("Number of Hits", labelpad=20, weight='bold', size=16)

    # Format y-axis label
    x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))


##dist of whole gbm prediction

from matplotlib.ticker import StrMethodFormatter
ax = urea_akid.hist(column='score', bins=100, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

ax = ax[0]
for x in ax:

    # Despine
    x.spines['right'].set_visible(False)
    x.spines['top'].set_visible(False)
    x.spines['left'].set_visible(False)

    # Switch off ticks
    x.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on")

    # Draw horizontal axis lines
    vals = x.get_yticks()
    for tick in vals:
        x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)

    # Remove title
    x.set_title("The Histogram of KS Association Predictions",size=20)

    # Set x-axis label
    x.set_xlabel("AKID Score", labelpad=20, weight='bold', size=16)

    # Set y-axis label
    x.set_ylabel("Number of Hits", labelpad=20, weight='bold', size=16)

    # Format y-axis label
    x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))



#### kinase distribution 
from matplotlib import pyplot as plt
kin_count= urea_true["kinase_domain"].value_counts().head(20)
print(kin_count)
kin_count= kin_count.sort_values()
kin_count.plot(kind='barh',color='#86bf91', figsize=(12,8))
plt.title("Top 20 Kinase",size=20)
plt.xlabel("Number of Hits", labelpad=20, weight='bold', size=18)
plt.ylabel("Kinase", labelpad=20, weight='bold', size=16)












