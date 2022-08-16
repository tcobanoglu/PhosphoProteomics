# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:24:39 2021

@author: Tugce Su
"""

#kinase id match
import pandas as pd
from csv import reader
import csv
"""
deepkinzero= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/DeepKinZero-master/Data/AllKinases.txt",sep='\t',header=None)
HGN_kinase=pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/kinases_HGNC_2016.tsv",sep='\t')
KinBase= pd.read_excel("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/Kinase/Table S1.xls")
KinBaseS2= pd.read_excel("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/Kinase/Table S2.xls", usecols=["Name","Ensembl"]).dropna()
"""
#HGN= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/kinase-id.txt", sep='\t',usecols=[""])

"""
deepkinzero["check"]=deepkinzero[0].isin(KinBase["Name"])
HGN_kinase["check_deepkinzero"]=HGN_kinase["Gene.Symbol.HGNC21Apr2016"].isin(deepkinzero[0])
HGN_kinase["check_kinbase"]=HGN_kinase["Gene.Symbol.HGNC21Apr2016"].isin(KinBase["Name"])
"""
path= "C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/chr-2021.txt"
with open(path, 'r') as inp:
    lines = inp.readlines()
    
with open('no_symbol.txt', 'w') as out:
    for line in lines:
        if not 'symbol withdrawn' in line:
            out.write(line)
            
