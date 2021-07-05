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

zone = 1
horizons = range(1,7)
colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']


gbm_rel = pd.read_csv(f'C:/Users/rosemaryt/Documents/GitHub/Forecast-combination/day ahead forecasts/zone{zone}_reliability.csv')


ax_1_coords = np.array([0,0,0,1,1,1])
ax_2_coords = np.array([0,1,2,0,1,2]) 
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10,5))

for ih in range(0, len(horizons)):
    h = horizons[ih]
    
    l0 = ax[ax_1_coords[ih],ax_2_coords[ih]].fill_between(gbm_rel['Nominal'], y1=gbm_rel['flat_lower'], y2=gbm_rel['flat_upper'], color=colours[4],alpha=0.7, ls='solid', lw=0)
    
    
    reldata = pd.read_csv(f'{fileloc}\\zone{zone}_h{h}_reliability.csv')
    plotdata = reldata[(reldata['Model']=='VAR')]
    l1 = ax[ax_1_coords[ih],ax_2_coords[ih]].fill_between(plotdata['Nominal'], y1=plotdata['flat_lower'], y2=plotdata['flat_upper'], color=colours[1],alpha=0.7, ls='solid', lw=0)
   
    plotdata = reldata[(reldata['Model']=='MC')]
    l2 = ax[ax_1_coords[ih],ax_2_coords[ih]].fill_between(plotdata['Nominal'], y1=plotdata['flat_lower'], y2=plotdata['flat_upper'], color=colours[2],alpha=0.7, ls='solid', lw=0)
    
    
    plotdata = reldata[(reldata['Model']=='EMD')]
    l1 = ax[ax_1_coords[ih],ax_2_coords[ih]].fill_between(plotdata['Nominal'], y1=plotdata['flat_lower'], y2=plotdata['flat_upper'], color=colours[3],alpha=0.7, ls='solid', lw=0)
    

fig.legend([l0,l1,l2,l3], labels=['GBM', 'VAR','MC','EMD'], loc="center right")
plt.subplots_adjust(top=0.98,
                    bottom=0.12,
                    left=0.1,
                    right=0.830,
                    hspace=0.220,
                    wspace=0.355)

for ih in range(0, len(horizons)):
    h = horizons[ih]
    
    ax[ax_1_coords[ih],ax_2_coords[ih]].scatter(gbm_rel['Nominal'], gbm_rel['flat_empirical'], color=colours[4], s=3)
    reldata = pd.read_csv(f'{fileloc}\\zone{zone}_h{h}_reliability.csv')
    plotdata = reldata[(reldata['Model']=='VAR')]
    ax[ax_1_coords[ih],ax_2_coords[ih]].scatter(plotdata['Nominal'], plotdata['flat_empirical'],s=3, color=colours[1])
    
    plotdata = reldata[(reldata['Model']=='MC')]
    ax[ax_1_coords[ih],ax_2_coords[ih]].scatter(plotdata['Nominal'], plotdata['flat_empirical'],s=3, color=colours[2])
    
    plotdata = reldata[(reldata['Model']=='EMD')]
    ax[ax_1_coords[ih],ax_2_coords[ih]].scatter(plotdata['Nominal'], plotdata['flat_empirical'],s=3, color=colours[3])
    
    
    
    ax[ax_1_coords[ih],ax_2_coords[ih]].axhline(0, c='k')
    ax[ax_1_coords[ih],ax_2_coords[ih]].set_xlim(0,1)
    ax[ax_1_coords[ih],ax_2_coords[ih]].set_ylim(-0.15, 0.2)
    ax[ax_1_coords[ih],ax_2_coords[ih]].set_yticks([ -0.1, 0, 0.1])
    ax[ax_1_coords[ih],ax_2_coords[ih]].text(x=0.3, y=0.17, s=f'Horizon={h}',size=14)
    

    
    if ih == 3:
        ax[ax_1_coords[ih],ax_2_coords[ih]].set_xlabel('Nominal')
        ax[ax_1_coords[ih],ax_2_coords[ih]].set_ylabel('Relative Empirical')





