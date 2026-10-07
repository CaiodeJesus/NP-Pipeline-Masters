#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 17:31:07 2025

@author: caio
"""

import networkx as nx
import pandas as pd
from metrics_lib import Get_metrics_data, adapt_yeast
from assertivity_test import Ranking
import matplotlib.pyplot as plt
from sklearn.metrics import auc as AUC
import os
import matplotlib as mpl
import numpy as np


def LCC_curve(G, metric, #dataframe with nodes|metric_name
              step = 1):
    metric['node'] = list(G.nodes())
    metric.sort_values(by = ['rank'], ascending = True, inplace = True)   
    ranked_by_metric =  list(metric['node'])

    n_list = []
    size_lcc_list = []
    for n in range(0,len(ranked_by_metric)-1,step):
        
        topN = ranked_by_metric[:n+1]
        
        S = G.copy()
        S.remove_nodes_from(topN)
        lcc = max(nx.connected_components(S), key=len)
        size_of_lcc = len(lcc)
        
        n_list.append(n)
        size_lcc_list.append(size_of_lcc)
    

    return n_list, size_lcc_list

def LCC_curves(metrics, G, 
               ax = plt.gca(), title = '',
               fs = 15, do_print = False,
               lw = 2, step = 1, cmap = 'tab20'):
 
    columns = list(metrics.columns)

    result = []
    if do_print:
        print()
        print('Ideal is 1.0')
        
        
    n_lines = len(columns)
    colormap = mpl.colormaps[cmap]
    colors = colormap(np.linspace(0, 1, n_lines))

    i = 0
    
    for col in columns:
        metric = pd.DataFrame(metrics[col])
        metric['rank'] = Ranking(metric[col])

        n_list, size_lcc_list = LCC_curve(G, metric,
                                          step = step)
        n_list = np.array(n_list)
        ax.plot(100*n_list /len(metrics) , size_lcc_list, '--', lw = lw,
                label = col, color = colors[i])

        i+=1
        auc = AUC(n_list, size_lcc_list)
    
        result.append([col, auc])
        if do_print:
            print(col)
            print()


        
    
    ax.set_title(title, size = fs+10)


    
    ax.tick_params(labelsize = fs-5)
    ax.set_xlabel('Fraction of nodes removed (%)', fontsize = fs)
    ax.set_ylabel('Size of the LCC', fontsize = fs)
    
    ax.legend(bbox_to_anchor=(-0.3, 1.05), fontsize = fs-5)
    
    result = pd.DataFrame(result,
                          columns = ['Metric', 'Destruction(AUC)'])
    result.set_index('Metric', inplace = True)
    
    return result
    
if __name__ == "__main__":

    # graph = 'trpA'
    # metrics = Get_metrics_data(graph)


    # path_to_graph = '../../../'
    # path_to_graph += 'graphs_python/aplicando networkx/toy_networks/graphml/'
    # # print(os.listdir(path_to_graph))
    # path_to_graph += f'{graph}.graphml'
    # G = nx.read_graphml(path_to_graph).to_undirected()
    
    
    
    # LCC_curves(metrics, G)
    
    
    
    # plt.figure(figsize = (10,4))
    
    # data_folder = '../data/subnets_table'
    # subnets = os.listdir(data_folder)
    # graph = subnets[0]


    # metrics, gabarito = adapt_yeast(graph, data_folder)
    # metrics.pop('gabarito')
    
    # graph = graph.replace('.csv','')
    # graph = '_'.join(graph.split('_')[:4])
    # path_to_graph = '../data/subnets/'
    # # print(graph)
    # for to_graph in os.listdir(path_to_graph):
    #     if graph in to_graph:
    #         # print(to_graph)
    #         graph = to_graph
        
    # path_to_graph += f'{graph}'
    # G = nx.read_graphml(path_to_graph).to_undirected()
    
  
    # result = LCC_curves(metrics, G)





    plt.figure(figsize = (10,4))


    data_folder = '../metanalysis/chimeras'
    subnets = os.listdir(data_folder)
    graph = 'CL_137_18.12_combined_score.csv'

    metrics = Get_metrics_data(graph, data_folder)

    graph = graph.replace('.csv','')
    graph = '_'.join(graph.split('_')[:4])
    path_to_graph = '../data/subnets/'
    # print(graph)
    for to_graph in os.listdir(path_to_graph):
        if graph in to_graph:
            # print(to_graph)
            graph = to_graph
        
    path_to_graph += f'{graph}'
    G = nx.read_graphml(path_to_graph).to_undirected()
    
  
    result = LCC_curves(metrics, G, step = 100, cmap = 'managua')




#Normalize by network size