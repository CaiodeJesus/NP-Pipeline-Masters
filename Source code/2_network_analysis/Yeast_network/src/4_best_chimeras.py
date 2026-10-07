#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 14:32:40 2025

@author: caio
"""


import sys
sys.path.append("src")

#homebrew
from correlation_test import Tau_Heatmap
from contribution_test import PCA_plot
from assertivity_test import ROC_curves
from destruction_test import LCC_curves
from metrics_lib import Get_metrics_data, adapt_yeast, prep_plot

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import time








if __name__ == "__main__":


    # data_folder = '../data/subnets_table'
    # subnets = os.listdir(data_folder)
    
    
    start = time.time()

    # for graph in subnets[:]:
    data_folder = '../string_channels/results/'

    graph = 'combined_score'
    path_to_graph = f'../string_channels/{graph}.graphml'

    metrics, gabarito = adapt_yeast('chimeras_combined_score_step100'+'.csv', data_folder)
    metrics.dropna(axis='columns', inplace = True)


    G = nx.read_graphml(path_to_graph).to_undirected()
    
    
    for col in metrics.columns:
        # print(col, metrics[col].sum())
        if metrics[col].sum() == 0.0:
            metrics.pop(col)
        
    
    fs = 50
    fig, ax_list = prep_plot(2, 2, scale_fig = 20, pw = 30, ph = 7)
    
    metrics.pop('gabarito')
    
    correlation = Tau_Heatmap(metrics, ax = ax_list[0], 
                              cmap = 'BrBG',
                              title = 'Correlation\nKendal Tau Coefficient', fs = fs)
    
    print('Correlation test done!')
    
    destruction = LCC_curves(metrics, G, ax = ax_list[2],
               title = 'Diffusion\nSize of LCC', fs = fs, lw = 10,
               step = 100, do_print = True)
    
    print('Destruction test done!')
    
    
    assertivity = ROC_curves(metrics, gabarito, ax = ax_list[3],
                             title = 'Assertivity\nROC Curve',
                             fs = fs, lw = 10)
    
    print('Assertivity test done!')
    

    contribution = PCA_plot(metrics, ax = ax_list[1],
             title = 'Contribution\nPCA features', fs = fs)
    
    print('Contribution test done!')
    
    final_results = assertivity.copy()
    for df in [contribution, correlation, destruction]:
        for col in df.columns:
            final_results[col] = df[col]


    
    final_results.to_csv(f"../string_channels/results/chimera_{graph}.csv")
    # plt.savefig(f"../data/results/prancha/{graph.replace('.graphml','')}.svg")
    # plt.clf()
    
    plt.tight_layout()
    plt.savefig(f"../string_channels/results/chimeras_{graph}.svg")

    end = time.time()
    print('',"Took %f s" % ((end - start) ) )

#To do list:
##top L (num de essenciais) de acordo com cada metrica no print (numero de acertos)
##função consenso