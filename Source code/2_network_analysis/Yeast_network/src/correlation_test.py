#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 13:47:38 2025

@author: caio
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from metrics_lib import Get_metrics_data, adapt_yeast



def Tau_Heatmap(metrics, cmap,
                ax = plt.gca(), title = '',
                fs = 15, do_print = False):    
    
    columns = list(metrics.columns)
    tau_matrix = []
    for i in range(len(columns)):
        tau_line = []
        for j in range(len(columns)):
            coli, colj = columns[i], columns[j]
            # print(coli, colj)
    
            tau = stats.kendalltau(metrics[coli], metrics[colj])
            tau = tau.statistic
            tau_line.append(tau)
        tau_matrix.append(tau_line)
    
    
    tau_matrix = np.array(tau_matrix)

    
    imshow = ax.imshow(tau_matrix,
               cmap = cmap, clim = (-1,1))
    
    ax.set_title(title, size = fs+10)

    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(mappable = imshow,
                        cax = cax)
    cbar.ax.tick_params(labelsize=fs-5)
    # cbar.ax.set_ylabel('Kendall Tau coefficient',
    #                    rotation = 270, fontsize = fs)
    ax.set_xticks(np.arange(len(columns)),
                  columns[:], size = fs-5,
                  rotation = 45, rotation_mode = 'anchor',
                  horizontalalignment='right', verticalalignment='top')
    ax.set_yticks(np.arange(len(columns)),
                  columns[:], size = fs-5)
    
    results = []
    for i in range(len(columns)):
        if do_print:
            print(columns[i], np.abs(tau_matrix[i]).sum())
            print(np.abs(tau_matrix[i]))
        results.append([columns[i], np.abs(tau_matrix[i]).sum()])
        
    results = pd.DataFrame(results,
                           columns = ['Metric', 'Correlation (Sum of |Tau|)'])
    
    results.set_index('Metric', inplace = True)

                           
    return results






if __name__ == "__main__":

    # graph = 'trpA'
    # metrics = Get_metrics_data(graph)
    
    # Tau_Heatmap(metrics, cmap = 'BrBG', title = '')





    # data_folder = '../data/subnets_table'
    # subnets = os.listdir(data_folder)
    # graph = subnets[-1]


    # metrics, gabarito = adapt_yeast(graph, data_folder)
    # metrics.dropna(axis='columns', inplace = True)
    
    # results = Tau_Heatmap(metrics, cmap = 'BrBG', title = '', do_print = True)
    
    
    
    
    
    
    
    
    
    
        
    data_folder = '../metanalysis/chimeras'
    subnets = os.listdir(data_folder)
    graph = 'CL_137_18.12_combined_score.csv'

    metrics = Get_metrics_data(graph, data_folder)

    results = Tau_Heatmap(metrics, cmap = 'BrBG', title = '', do_print = True)
 