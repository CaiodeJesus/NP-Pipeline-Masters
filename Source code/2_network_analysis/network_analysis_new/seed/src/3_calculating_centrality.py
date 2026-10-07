#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 14:56:04 2025

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
from metrics import Best_conf, Centrality_Metrics







with open("../input_config") as f:    
    input_lines = f.readlines()

metrics_to_use = []
is_metric = False
for line in input_lines:
    if is_metric:
        metric_name = line.strip('-').strip('\n')
        metrics_to_use.append(metric_name.lower())
    if 'Metrics To be Used' in line:
        is_metric = True

print(metrics_to_use)


today = datetime.datetime.now()
if __name__ == "__main__":
    
    # graph = 'database'
    for infile in os.listdir('../data/')[:]:
        
        if '.' in infile:
            continue
        
        infile_folder = f'../data/{infile}/'
        print(infile)
        
        for file in os.listdir(infile_folder + 'channels/')[:]:
            
            if '.graphml' not in file:
                continue
            
            # print(file)

            if 'subnet_tables' not in os.listdir(infile_folder):
                os.mkdir(infile_folder+'subnet_tables')
                
                
            start = time.time()
            graph = file.replace('.graphml','')
            print(graph)
            
            
            path_to_graph = infile_folder + 'channels/'
            # # print(os.listdir(path_to_graph))
            path_to_graph += f'{graph}.graphml'
            G = nx.read_graphml(path_to_graph).to_undirected()
            
            
            node_names = list(G.nodes())
            # # node_names = [name[1] for name in node_names]
            print('','size of network',len(node_names))
            
            if 'best' in metrics_to_use:
                metrics = Best_conf(G)
            
            else:
                metrics = Centrality_Metrics(G, metrics_to_use)
                
            metrics['Node_names'] = node_names
            metrics.set_index(['Node_names'], inplace = True)
            
            
            end = time.time()
            print('',"Took %f s" % ((end - start) ) )
            day = today.strftime("%y_%m_%d")
            
            num_metrics = len(metrics.columns)
            
            print(graph)
            print(f'Using {num_metrics} metrics in {day}.')
            
            channel = '_'.join(graph.split('_')[1:])
            out_name = f'{infile_folder}/subnet_tables/{channel}.csv'
           
            print(f'Saving dataframe to {out_name}')
            
            metrics.to_csv(out_name)
            

    


    
    with open('../data/centrality.log', 'w') as f:
        f.write(today.strftime("%c")+'\n')
        f.write('Metrics Used:\n-')
        f.write( '\n-'.join(list(metrics.columns)) )
        
        