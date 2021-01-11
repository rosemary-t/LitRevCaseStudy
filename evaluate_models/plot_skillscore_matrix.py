# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 15:47:21 2021

@author: rosemaryt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True) #LaTeX fonts on plots for compatibility with thesis etc
plt.rc('font', family='serif', size=16)
fileloc = r'C:\Users\rosemaryt\Documents\GitHub\LitRevCaseStudy\evaluate_models'

matrixdata = pd.read_csv(f'{fileloc}\\skillscore_matrix_h2PB.csv')
matrixdata.set_index('Unnamed: 0', inplace=True)
axlabs = matrixdata.columns.values
plotdata = np.array(matrixdata)

fig, ax = plt.subplots()
im = ax.imshow(plotdata, cmap='RdYlGn')
ax.set_xticks(np.arange(len(axlabs)))
ax.set_yticks(np.arange(len(axlabs)))
ax.set_xticklabels(axlabs)
ax.set_yticklabels(axlabs)


for i in range(len(axlabs)):
    for j in range(len(axlabs)):
        text = ax.text(j, i, round(plotdata[i, j],2),
                       ha="center", va="center", color='k')
plt.show()
