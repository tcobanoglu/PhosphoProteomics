# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:23:29 2022

@author: Tugce Su
"""
#http://kinase.com/kinbase/FastaFiles/
import pandas as pd
#gbm 
#make all kinase upper case all the time
gbm_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/gbm-akid-result.txt", delimiter=r"\s+")
gbm_akid['kinase_domain'] = gbm_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
gbm_akid['kinase_domain']= gbm_akid['kinase_domain'].str.upper()
kinase= gbm_akid["kinase_domain"].unique()
l_kin= set(kinase.tolist()) #akid kinase
psp= gbm_akid["peptide"].unique()
akid_kin= pd.DataFrame(kinase,columns=["akid_kin"])

#KS relations
KS_relations = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/psp_and_nwk_unique_KSrelations_inka1.0.tsv", sep= "\t",low_memory=False)
KS_relations["Kinase"]= KS_relations["Kinase"].str.upper()
kin_ks=KS_relations["Kinase"].unique()
l_kin_ks= set(kin_ks.tolist())

##ks difference
ks_diff= l_kin - l_kin_ks
ks_same= l_kin & l_kin_ks

## are they human kinase or group? 

hum_kin = pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/dk0-2016-2021.txt", sep= "\t")
hum_kin["HGNC 2021"]= hum_kin["HGNC 2021"].str.upper()
hum_kin_2021= set(hum_kin["HGNC 2021"].tolist())
same= hum_kin_2021 & l_kin
diff= l_kin - hum_kin_2021
hum_kin["kinase(dk0)"]=hum_kin["kinase(dk0)"].str.upper()
##D0 AKID kinase
d0_kin= hum_kin["kinase(dk0)"].str.upper()
d0_kin= d0_kin.unique()
d0_kin= set(d0_kin.tolist())
d0_same= d0_kin & l_kin
d0_diff= l_kin- d0_kin
#match known D0 kinase and AKID
akid_merge= pd.merge(hum_kin,akid_kin, left_on="kinase(dk0)",right_on="akid_kin",how="right" )
file = akid_merge.to_csv("akid-d0-hgnc-kinase",sep="\t",index=False)




#####  MATCH HGNC 



##akid kinase names==allias symbol need to match them approved name
alias= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/allias-approved-hgnc.txt",header=None,sep="\t",skiprows=1)
alias_pd=alias.iloc[:,0].str.split(" ",n=1, expand=True)



ali=open("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/allias-approved-hgnc.txt","r")
line=ali.readlines()
#line=line.split(" ")
#lines=[]
lines=[l.split(" ") for l in line]
lines = [' '.join(i).split() for i in lines]
ids={}#keys are approved name and values are alias
for l in lines:
    key,value=l[0],l[1:]
    ids[key]= value
l_kin=list(l_kin)
akid_kin=[]
for x,y in ids.items():
   for i in y:
       for l in l_kin:
           if i == l:
               akid_kin.append([i,x])
              
          
#####    
#Finally compare KS relation and akid
####

#arrange gbm psp
gbm_akid['peptide'] = gbm_akid['peptide'].str.replace("_", '-')
KS_relations['diff_sub'] = gbm_akid.peptide.isin(KS_relations.Accession_AApos)      
print(KS_relations["diff_sub"].value_counts())

#open new human kin file
kin= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/akid-d0-hgnc-kinase.txt",sep="\t")
#merge akid kinase and hgnc 2016
KS_relations_akid= pd.merge(KS_relations, kin, left_on= "Kinase",right_on="HGNC 2016",how="inner")
print(KS_relations_akid.columns)
KS_relations_akid=KS_relations_akid.drop(['Kinase','HGNC 2016','HGNC 2021', 'kinase(dk0)', 'Approved symbol'],axis=1)

#ks check 
KS_relations_akid["KS"]= KS_relations_akid["Accession_AApos"] +"_" +KS_relations_akid["akid_kin"]
gbm_akid["KS"]= gbm_akid["peptide"] + "_" + gbm_akid["kinase_domain"]
KS_relations_akid['diff_ks'] = gbm_akid.KS.isin(KS_relations_akid.KS)
print(KS_relations_akid["diff_ks"].value_counts())
gbm_akid['diff_ks'] = gbm_akid.KS.isin(KS_relations_akid.KS)
print(gbm_akid["diff_ks"].value_counts())



##gbm true
gbm_true = gbm_akid[gbm_akid['diff_ks'] == True]# same KS with KS relations
psp=gbm_true["peptide"].unique()

##covarage of 1/3 gbm
import matplotlib.pyplot as plt
import numpy as np

y=np.array([59.80,40.19])
mylabels = ["Identical KS", "Different KS"]
mycolors = ["c","y"]
plt.pie(y, labels = mylabels, colors = mycolors,autopct='%1.1f%%',shadow=True, startangle=90)
plt.title("The Covarage of Annotated KS in GBM Dataset", pad=32)
plt.show() 


##distribution of known 1/3
from matplotlib.ticker import StrMethodFormatter
ax = gbm_true.hist(column='score', bins=25, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

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
ax = gbm_akid.hist(column='score', bins=100, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9, range=[0,1.5])

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

plt.savefig('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/Figures/gbm_score_range.pdf', format='pdf')
#

### boxplot to see the range of data
import seaborn as sns
bplot = sns.boxplot(x='score', y="kinase_domain",
                 data=gbm_true, 
                 width=0.3,
                 palette="colorblind")


#### kinase distribution 

kin_count= gbm_true["kinase_domain"].value_counts().head(20)
print(kin_count)
kin_count= kin_count.sort_values()
kin_count.plot(kind='barh',color='#86bf91', figsize=(12,8))
plt.title("Top 20 Kinase",size=20)
plt.xlabel("Number of Hits", labelpad=20, weight='bold', size=18)
plt.ylabel("Kinase", labelpad=20, weight='bold', size=16)

##kin dist for whole dataset
kin_count= gbm_akid["kinase_domain"].value_counts()
kin_count= kin_count.sort_values()
kin_count.plot(kind='barh',color='#86bf91')
plt.title("Top 20 Kinase",size=20)
plt.xlabel("Number of Hits", labelpad=20, weight='bold', size=18)
plt.ylabel("Kinase", labelpad=20, weight='bold', size=16)

##look ind kinase


top5=gbm_true[~gbm_true["kinase_domain"].isin(["EGFR","MET","FYN","ABL1","FGR"])]
top5.boxplot(column="score",by="kinase_domain")



import seaborn as sns
bplot = sns.boxplot(x='score',
                 data=gbm_akid, 
                 width=0.3,
                 palette="colorblind")


egfr= gbm_true.loc[gbm_true["kinase_domain"]=="ABL1"]

ax = egfr.hist(column='score', bins=50, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

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
    x.set_title("ABL1",size=20)

    # Set x-axis label
    x.set_xlabel("AKID Score", labelpad=20, weight='bold', size=16)

    # Set y-axis label
    x.set_ylabel("Number of Hits", labelpad=20, weight='bold', size=16)

    # Format y-axis label
    x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))






#let's say 0.8
#low 0.25 medium 0.6 high 0.8
gbm_08= gbm_akid.loc[gbm_akid["score"] >= 1]
psp_hg= gbm_08["peptide"].unique()


















