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

zone = 4
horizon = 2
reldata = pd.read_csv(f'{fileloc}\\zone{zone}_h{horizon}_reliability.csv')

models = reldata.Model.unique()
colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
#colours = [332288, 88CCEE, 44AA99, 117733, 999933, DDCC77, CC6677, 882255, AA4499]
#colours = ['C{}'.format(i) for i in range(len(models))]

plt.figure()
plt.axhline(0, c='k')
plt.xlim(0,1)
plt.ylim(-0.15, 0.15)
plt.xlabel('Nominal')
plt.ylabel('Relative Empirical')

for i in range(len(models)):
    m = models[i]
    plotcol = colours[i]
    plotdata = reldata[reldata['Model']==m]
    plt.fill_between(x=plotdata['Nominal'], y1=plotdata['flat_lower'], y2=plotdata['flat_upper'], alpha=0.7, color=plotcol, label=m, ls='solid', lw=0)
    plt.scatter(plotdata['Nominal'], plotdata['flat_empirical'], color=plotcol, s=5)
    

plt.legend(loc=4)
plt.show()