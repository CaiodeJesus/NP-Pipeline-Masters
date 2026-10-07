#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 17:08:44 2025

@author: caio
"""


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



df_servers = pd.read_csv('../data/support/disease_servers.csv')
servers = df_servers['Name']
threshold_values = df_servers['Threshold_value']

names_metrics = ['Metric' + m for m in df_servers['Chosen_Metric'].fillna('')]

fs=40
L=5

color_diseases = ['#f8e16c',
                   '#00c49a',
                   '#156064']
other_disease = 'chagas'

for i in range(len(df_servers[:])):

    s = servers[i]

    print(s)
    d=0
    df = pd.read_csv('../data/clean/disease/leish/'+s+'.csv')
    
    mask = np.array(['Metric' in c for c in df.columns])
    metric_cols = np.array(df.columns[mask])
    

    if len(metric_cols)%2 == 0:        
        num_cols = len(metric_cols)//2
        num_rows = 2

    else:
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
            # p = 0.9
            print(p, metric.max())
        
            
            number_of_passed = (metric > p).sum()
            percentage_of_passed =  100 * (number_of_passed / len(metric))
            print(number_of_passed, len(metric))
            text_about_discarted = '%.2f%% were above threshold'%(percentage_of_passed)
            
            weights = 100*np.ones_like(metric)/float(len(metric))
            hist = ax.hist(metric, weights=weights, color = 'r', label = 'Leish out')
            count, bins = hist[0], hist[1]      
            
            
            ax.plot([p,p],[0,50],'k--',lw=5)

            if len(bins[bins>=p]) != 0:

                ax.stairs(count[bins[:-1]>=p], bins[bins>=p],
                          fill = True, linewidth = 0,
                          color = 'g', label = 'Leish in')

            else:
                continue

            if len(bins[bins<p]) != 0:

                width_bins = bins[-1] - bins[-2]
                
                initial_bin = bins[bins<p][-1]
                width_cover_initial_bin = initial_bin + width_bins - p
                height_cover_initial_bin = count[bins[:-1]<p][-1]
                ax.bar( (p + width_cover_initial_bin/2),
                        height_cover_initial_bin,
                        color='g', width = width_cover_initial_bin, label = 'Leish in')
        else:
            weights = 100*np.ones_like(metric)/float(len(metric))
            hist = ax.hist(metric, weights=weights, color = color_diseases[0], label = 'Leish')
        
        df = pd.read_csv('../data/clean/disease/'+other_disease+'/'+s+'.csv')
        metric = df[col]
        
        weights = 100*np.ones_like(metric)/float(len(metric))
        hist = ax.hist(metric, weights=weights,
                        edgecolor = color_diseases[d+1],
                        label = other_disease.capitalize(),
                        fill = False, linewidth = 5)


            
    
        print()
                
        
            
    
        plt.legend(loc='best', fontsize = fs)
    plt.suptitle(s + '\n' + text_about_discarted,
                 fontsize = fs+5, y = 1.0)
    plt.tight_layout()
    plt.savefig('../data/support/histograms/disease/'+s + '.svg')
    plt.clf()