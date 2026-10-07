#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:07:15 2025

@author: caio
"""


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable





df_servers = pd.read_csv('../data/support/compound_servers.csv')
servers = df_servers['Name']
threshold_values = df_servers['Threshold_value']

names_metrics = ['Metric' + m for m in df_servers['Chosen_Metric'].fillna('')]

fs=40
L=5


for i in range(len(df_servers[:])):

    s = servers[i]

    print(s)
    df = pd.read_csv('../data/clean/compound/'+s+'.csv')
    
    mask = np.array(['Metric' in c for c in df.columns])
    metric_cols = np.array(df.columns[mask])
    

    
    num_cols = len(metric_cols)
    num_rows = 1
    fig = plt.figure(figsize = (4*L*num_cols,3*L*num_rows))


    index = 1
    for col in metric_cols[:]:

        
        ax = fig.add_subplot(num_rows, num_cols, index)
        index+=1
        ax.tick_params(size = fs-5)
        ax.set_ylabel('Percentage (%)', fontsize = fs)
        ax.set_xlabel(col, fontsize = fs)
        ax.set_ylim(0,100)
        ax.tick_params(labelsize = fs)
        
        metric = df[col]

        if col == names_metrics[i]:
            df.dropna(subset = ['UniProt',col], inplace = True)
            metric = df[col]
            
            p = threshold_values[i]
        
            number_of_passed = (metric > p).sum()
            percentage_of_passed =  100 * (number_of_passed / len(metric))
            text_about_discarted = '%.2f%% were above threshold'%(percentage_of_passed)
            
            weights = 100*np.ones_like(metric)/float(len(metric))
            hist = ax.hist(metric, weights=weights, color = 'r')
            count, bins = hist[0], hist[1]      
            
            
            ax.plot([p,p],[0,50],'k--',lw=5)

            if len(bins[bins>=p]) != 0:

                ax.stairs(count[bins[:-1]>=p], bins[bins>=p],
                          fill = True, linewidth = 0,
                          color = 'g')

            else:
                continue

            if len(bins[bins<p]) != 0:

                width_bins = bins[-1] - bins[-2]
                
                initial_bin = bins[bins<p][-1]
                width_cover_initial_bin = initial_bin + width_bins - p
                height_cover_initial_bin = count[bins[:-1]<p][-1]
                ax.bar( (p + width_cover_initial_bin/2),
                        height_cover_initial_bin,
                        color='g', width = width_cover_initial_bin)
        else:
            weights = 100*np.ones_like(metric)/float(len(metric))
            hist = ax.hist(metric, weights=weights, color = 'b')
            
    
        # print()
                
        
            
    
        
    plt.suptitle(s + '\n' + text_about_discarted,
                 fontsize = fs+5, y = 1.0)
    plt.tight_layout()
    plt.savefig('../data/support/histograms/compound/'+s + '.svg')
    plt.clf()

