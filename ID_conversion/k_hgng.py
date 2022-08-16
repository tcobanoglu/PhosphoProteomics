# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:18:18 2021

@author: Tugce Su
"""

import pandas as pd
import csv

ids= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/no_symbol.txt", sep='\t')
deepkinzero= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA/week4/DeepKinZero-master/Data/AllKinases.txt",sep='\t',header=None)
deepkinzero= deepkinzero.rename(columns={0: "kinase", 1: "number", 2: "uniprot"})
HGN_kinase=pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/kinases_HGNC_2016.tsv",sep='\t')
new_match= pd.merge(deepkinzero,ids, left_on="uniprot",right_on="UniProt ID(supplied by UniProt)",how="left")
"""
print(new_match.columns)
new_match= new_match.drop(columns=['NCBI Gene ID', 'UniProt ID(supplied by UniProt)','number','Chromosome'])
new_match=new_match.rename(columns={'kinase': "DeepKinZero", 'Approved symbol': "HGNC 2021"})
new_match.to_csv("deepkinzero_HGNC_2021.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)
#print(HGN_kinase.columns)
#print(new_match.columns)
"""
merge_ids= pd.merge(HGN_kinase,ids, left_on='Protein.Name.HGNC21Apr2016',right_on='Approved name',how="left")
print(merge_ids.columns)

approved_id= pd.merge(HGN_kinase,ids, left_on='Gene.Symbol.HGNC21Apr2016',right_on='Approved symbol',how="left")
approved_id.to_csv("2016-2021.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)
"""
merge_ids=merge_ids.rename(columns={ 'Approved symbol': "HGNC 2021",'Gene.Symbol.HGNC21Apr2016':"HGNC 2016"})
kinase_ids= merge_ids.drop(columns=['NCBI Gene ID', 'UniProt ID(supplied by UniProt)','number','in_PSP', 'in_NWK',
       'Kinase_with_Substrate', 'Kinase_without_Substrate'])

#kinase_ids.to_csv("kinase_ids.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)

#match without deepkinzero
"""

##match 2016-2021-deepkinzeroids

kinase_ids=pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-2/kinase_ids202116.txt",sep='\t')
dk0_2021= new_match[["kinase",'Approved symbol']]
dk0_2021_2016= pd.merge(kinase_ids,dk0_2021, left_on='HGNC 2021',right_on='Approved symbol',how="outer")
dk0_2021_2016.to_csv("dk0-2016-2021.txt", sep="\t", quoting=csv.QUOTE_NONE,index=False)