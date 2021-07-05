# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:32:24 2021

@author: rosemaryt
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True) #LaTeX fonts on plots
plt.rc('font', family='serif', size=16)

models = ["VAR", "MC", "EMD"]
colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
ss_type = "RMSE"
if ss_type=="pb":
    long_ss = "Pinball"
else: 
    long_ss = ss_type

path = r"C:\Users\rosemaryt\Documents\GitHub\LitRevCaseStudy\evaluate_models"

plt.figure()

plt.axhline(0, c='k')

for i in range(len(models)):
    m = models[i]
    plotcol = colours[i+1] # so the colours of each model match that seen in reliability diagram.
    plotdata = pd.read_csv(f'{path}/{m}_skillscores_horizons.csv')
    horizons = plotdata[plotdata['type']=='mean']['Horizon']
    plt.fill_between(x=horizons, y1=plotdata[plotdata['type']=='lower'][f'{ss_type}_ss'], y2=plotdata[plotdata['type']=='upper'][f'{ss_type}_ss'], alpha=0.7, color=plotcol, label=m, ls='solid', lw=0)
    plt.scatter(horizons, plotdata[plotdata['type']=='mean'][f'{ss_type}_ss'], color=plotcol, s=5)
    
plt.xlabel('Horizon (hrs)')
plt.ylabel(f'{long_ss} skill score')
plt.legend(loc=4)

