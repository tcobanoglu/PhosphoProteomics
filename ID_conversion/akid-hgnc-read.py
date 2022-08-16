# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 11:05:04 2022

@author: Tugce Su
"""

import pandas as pd
import csv

akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid-inka.txt",sep='\t')
print(akid.loc[176,:])

print(akid.isnull().sum(axis = 0))


akid1 = open('C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/akid-inka.txt', 'r')
akid1 = akid1.readlines()
akid_list = [line.split("\t") for line in akid1]

count_kin=0
count_sub=0
kinase= 0
for i in range(len(akid_list)):

   
    if len(akid_list[i][7]) < 3 :
        #print(akid_list[i][7])
        count_kin= count_kin +1
       #print(akid_list[i])
    if len(akid_list[i][8]) < 3 :
        count_sub= count_sub +1
        print(akid_list[i])

i = 2
for i in range(len(akid_list)):
    if len(akid_list[i][1])< 5:
        if len(akid_list[i][1]) > 2: 
            akid_list[i -1][8] = akid_list[i][1]
            del akid_list[i]
            print(akid_list[i-1])
y=2  
for y in range(len(akid_list)):
    if len(akid_list[y][1]) < 3:
        del akid_list[y]
        

d=" " 
list1=[] 
for a in akid_list:
    l=d.join(a)
    list1.append(l)

print(list1[1])   
#d= {"col1":akid_list}
df=pd.DataFrame(akid_list)
df['col1'] = df['col1'].astype("string")
#print(df.columns)
df=df["col1"].str.split(' ', expand=True,n=8)
df = df.iloc[:,:-4]

print(df.iloc[29534,:])


df.to_csv("akid_uni_match1.txt",index=False,sep='\t')

df.columns = df.iloc[0] 
df = df[1:]
print(df.isnull().sum(axis = 0))












