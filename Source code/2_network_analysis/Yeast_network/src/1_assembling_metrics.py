#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 16:35:51 2025

@author: caio
"""

import networkx as nx
import pandas as pd
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bct
from cdlib import algorithms
import datetime
import time
from metrics_lib import *



if __name__ == "__main__":

    
    #Since metrics were already calculated (for a month)
    
    metrics_files = os.listdir('../string_channels/metrics/')
    for graph in os.listdir('../string_channels/graphs/'):
        graph = graph.replace('.graphml','')
        print(graph)
        files_of_this_graph = []
        for file in metrics_files:
            if graph in file:
                print('',file)
                files_of_this_graph.append(file)
        
        most_recent = sorted(files_of_this_graph)[-1]
        print()
        print('',most_recent)
        
        
        
        
              
    
    

    # # list_graphs = os.listdir('../data/subnets/')
    # list_graphs = ['CL_137_18.12_combined_score.graphml']
    # # graph = 'database'
    # for graph in list_graphs:
    #     start = time.time()
    #     graph = graph.replace('.graphml','')
    #     print(graph)
    #     path_to_graph = '../test/'
    #     # print(os.listdir(path_to_graph))
    #     path_to_graph += f'{graph}.graphml'
    #     G = nx.read_graphml(path_to_graph).to_undirected()
        
        
    #     node_names = list(G.nodes())
    #     # node_names = [name[1] for name in node_names]
    #     print('','size of network',len(node_names))
        
        
    #     metrics = Centrality_Metrics(G)
    #     metrics['Node_names'] = node_names
    #     metrics.set_index(['Node_names'], inplace = True)
        
        
    #     end = time.time()
    #     print('',"Took %f s" % ((end - start) ) )
    #     day = datetime.datetime.now()
    #     day = day.strftime("%y_%m_%d")
        
    #     num_metrics = len(metrics.columns)
        
    #     print(graph)
    #     print(f'Using {num_metrics} metrics in {day}.')
        
    #     out_name = f'../test/{graph}_{day}_{num_metrics}m.csv'
        
    #     print(f'Saving dataframe to {out_name}')
        
    #     metrics.to_csv(out_name)
    
    
    
    
    
    