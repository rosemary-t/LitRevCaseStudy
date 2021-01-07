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
ss_type = "RMSE"
if ss_type=="pb":
    long_ss = "Pinball"
else: 
    long_ss = ss_type

path = r"C:\Users\rosemaryt\Documents\GitHub\LitRevCaseStudy\evaluate_models"

plt.figure()
plt.axhline(0, c='k')

for m in models:
    plotdata = pd.read_csv(f'{path}/{m}_skillscores_horizons.csv')
    error_list = [np.array(plotdata[plotdata['type']=='mean'][f'{ss_type}_ss']) - np.array(plotdata[plotdata['type']=='min'][f'{ss_type}_ss']), np.array(plotdata[plotdata['type']=='max'][f'{ss_type}_ss']) - np.array(plotdata[plotdata['type']=='mean'][f'{ss_type}_ss'])]
    # plt.plot(plotdata[plotdata['type']=='mean']['Horizon'], plotdata[plotdata['type']=='mean'][f'{ss_type}_ss'])
    plt.errorbar(plotdata[plotdata['type']=='mean']['Horizon'], plotdata[plotdata['type']=='mean'][f'{ss_type}_ss'], yerr=error_list, label=f'{m}', capsize=5)


plt.xlabel('Horizon (hrs)')
plt.ylabel(f'{long_ss} skill score')
plt.legend(loc=4)

