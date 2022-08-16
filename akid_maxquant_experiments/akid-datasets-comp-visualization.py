# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:41:19 2022

@author: Tugce Su
"""

# library
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Use the venn2 function
plt.figure(figsize = (6, 5))

v= venn2(subsets = (1402971, 189043,27113), set_labels = ('KS Relations Dataset', 'AKID Dataset'),alpha=0.9,set_colors=("pink", "skyblue", "green"))
plt.title("Unique Phosphosites Comparison")