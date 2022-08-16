# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:53:31 2022

@author: Tugce Su
"""

import pandas as pd
from mlxtend.evaluate import permutation_test
c_akid= pd.read_csv("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/AKID/4cell-akid-result.txt", delimiter=r"\s+")
c_akid['kinase_domain'] = c_akid['kinase_domain'].str.replace("_Hsap_domain1", '')
c_akid['kinase_domain']= c_akid['kinase_domain'].str.upper()
kinase= c_akid["kinase_domain"].unique()
psp_all= c_akid["peptide"].unique()
c9_akid= c_akid.loc[c_akid["score"] >= 0.99].reset_index()
c0_akid= c_akid.loc[c_akid["score"] < 0.99].reset_index()


#http://rasbt.github.io/mlxtend/user_guide/evaluate/permutation_test/
treatment= c0_akid["score"].tolist()
control= c9_akid["score"].tolist()

p_value = permutation_test(treatment, control,
                           method='approximate',
                           num_rounds=10000,
                           seed=0)
print(p_value) #9.999000099990002e-05

