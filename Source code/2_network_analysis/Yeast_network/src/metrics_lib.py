#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 14:49:59 2025

@author: caio
"""

import networkx as nx
import pandas as pd
import os
import numpy as np
import bct
from cdlib import algorithms
import matplotlib.pyplot as plt


def prep_plot(num_rows, num_cols,
              scale_fig=10, ratio = [1,1], pw = 0, ph = 0):
    
    """
    ratio = [rx, ry]
        width  ->  rx*scale_fig*num_cols
        height ->  ry*scale_fig*num_rows
    """
    
    

    fig = plt.figure(figsize = (ratio[0]*scale_fig*num_cols + pw,  #width
                                ratio[1]*scale_fig*num_rows + ph)) #height
    
    ax_list = []
    for i in range(num_cols*num_rows):
        ax_list.append(fig.add_subplot(num_rows, num_cols, i+1))

    # plt.clf()
    return fig, ax_list 


def ParticipationCoefficient(G, cd_alg = 'walktrap', communities = [], k=4):
    """
    G: Graph object of newtorkx
    cd_alg: algorithm of community detection
    communities: array of communities
        same lenght as quantity of nodes, says which group the node belongs
        communities[i] = <name of group that node i is in>
    k : girvan newman level (useful if cd_alg is girvan newman)
    
    
    
    
    EXAMPLE OF USAGE
    G = nx.read_graphml('toy_networks_nx/fig1B.graphml')
    
    #with given communities
    coms = [['0','1','2','9','10','11'],
            ['3','4','5'],
            ['6','7','8']]
    PC = list(ParticipationCoefficient(G, communities = coms)['PC'])

    
    #or with communities to be detected
    PC = list(ParticipationCoefficient(G, cd_alg = 'girvan_newman', k=3)['PC'])
    """

    node_list = list(G.nodes())

    df = pd.DataFrame()
    df['node'] = node_list

    if communities == []:
        G_num = nx.convert_node_labels_to_integers(G)
        node_list = list(G_num.nodes())

        if cd_alg == 'walktrap':
            communities = algorithms.walktrap(G_num).communities
        elif cd_alg == 'girvan_newman':
            communities = algorithms.girvan_newman(G_num,k-1).communities

        else:
            print('Unrecognized algorithm')
            
    elif len(G.nodes) != sum([len(i) for i in communities]):
        print('Check size of communities vector.')
    
    else:
        print('Tudo OK com as comunidades.')
    
    A = nx.adjacency_matrix(G).toarray()
    
    ci = []
    # print(node_list)
    # print(communities)
    for node in node_list:
        com_index = 0 
        for com in communities:
            # print(type(node), com)
            if node in com:
                 ci.append(com_index)
                 # print(node)
                 continue
            else:
                com_index +=1
                 
    # print(ci)
    participation_coeff_list = bct.centrality.participation_coef(A, ci)
    

    df['PC'] = participation_coeff_list
    
    return df



def Centrality_Metrics(G):
    centrality_metrics = []
    metric_names = []

    null_vals = np.zeros(len(G.nodes()))

    #DEGREE
    try:
        degree = nx.degree_centrality(G)
    except:
        degree = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no degree')
        
    centrality_metrics.append(degree)
    metric_names.append('degree')
    print('-',metric_names[-1])
    
    #EIGENVECTOR
    try:
        eigenvector = nx.eigenvector_centrality(G, max_iter = 800)
    except:
        degree = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no degree')        
    
    centrality_metrics.append(eigenvector)
    metric_names.append('eigenvector')
    print('-',metric_names[-1])

    #CLOSENESS
    try:
        closeness = nx.closeness_centrality(G)
    except:
        closeness = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no closeness')      
    
    centrality_metrics.append(closeness)
    metric_names.append('closeness')
    print('-',metric_names[-1])

    #CURRENT FLOW CLOSENESS
    try:
        current_flow_closeness = nx.current_flow_closeness_centrality(G)
    except:
        current_flow_closeness = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no current_flow_closeness')     
    
    centrality_metrics.append(current_flow_closeness)
    metric_names.append('current_flow_closeness')
    print('-',metric_names[-1])

    #BETWEENNESS
    try:
        betweeness = nx.betweenness_centrality(G) #shortest path
    except:
        betweeness = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no betweeness')      
    
    centrality_metrics.append(betweeness)
    metric_names.append('betweeness')
    print('-',metric_names[-1])

    #CURRENT FLOW BETWEENNESS
    try:
        cf_betweenness = nx.current_flow_betweenness_centrality(G)
    except:
        cf_betweenness = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no cf_betweenness')     
    
    centrality_metrics.append(cf_betweenness)
    metric_names.append('cf_betweenness')
    print('-',metric_names[-1])

    #COMMUNICABILITY BETWEENNESS
    try:
        comm_betweenness = nx.communicability_betweenness_centrality(G)
    except:
        comm_betweenness = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no comm_betweenness')   

    centrality_metrics.append(comm_betweenness)
    metric_names.append('comm_betweenness')
    print('-',metric_names[-1])

    #LOAD CENTRALITY
    try:
        load_centrality = nx.load_centrality(G)
    except:
        load_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no load_centrality')   

    centrality_metrics.append(load_centrality)
    metric_names.append('load_centrality')
    print('-',metric_names[-1])

    #HARMONIC CENTRALITY
    try:
        harmonic_centrality = nx.harmonic_centrality(G)
    except:
        harmonic_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no harmonic_centrality')      
    
    centrality_metrics.append(harmonic_centrality)
    metric_names.append('harmonic_centrality')
    print('-',metric_names[-1])

    #LAPLACIAN CENTRALITY
    try:
        laplacian_centrality = nx.laplacian_centrality(G)
    except:
        laplacian_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no laplacian_centrality')     
    
    centrality_metrics.append(laplacian_centrality)
    metric_names.append('laplacian_centrality')
    print('-',metric_names[-1])

    #PERCOLATION CENTRALITY
    try:
        percolation_centrality = nx.percolation_centrality(G)
    except:
        percolation_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no percolation_centrality')     
    
    centrality_metrics.append(percolation_centrality)
    metric_names.append('percolation_centrality')
    print('-',metric_names[-1])

    #SECOND ORDER CENTRALITY
    try:
        second_order_centrality = nx.second_order_centrality(G)
    except:
        second_order_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no second_order_centrality')      
    
    centrality_metrics.append(second_order_centrality)
    metric_names.append('second_order_centrality')
    print('-',metric_names[-1])
    
    #LOCAL REACHING CENTRALITY
    try:
        vals = [nx.local_reaching_centrality(G,v) for v in G.nodes()]
        local_reaching_centrality = dict(zip(G.nodes, vals ))
    except:
        local_reaching_centrality = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no degree')      
    
    centrality_metrics.append(local_reaching_centrality)
    metric_names.append('local_reaching_centrality')
    print('-',metric_names[-1])

    #PARTICIPATION COEFFICIENT
    try:
        vals = list(ParticipationCoefficient(G, cd_alg = 'girvan_newman', k=3)['PC'])
        PC = dict(zip(G.nodes, vals ))
    except:
        PC = dict(zip(G.nodes, null_vals))
        print('Algo deu errado no PC')    
    
    centrality_metrics.append(PC)
    metric_names.append('participation_coeff')
    print('-',metric_names[-1])

    
    df = pd.DataFrame()
    
    for i in range(len(metric_names)):
        # print(metric_names[i])
        df[metric_names[i]] = centrality_metrics[i]
    
    print()
    return df


def Get_metrics_data(graph, data_folder = '../data'):
    data_files = np.sort(os.listdir(data_folder))
    graph_files = data_files[ [graph in file for file in data_files] ]
    graph_files = graph_files[ [file[:len(graph)] == graph for file in graph_files]]
    chosen_file = graph_files[-1]
    
    print(f'Reading {chosen_file}')
    
    metrics_src = pd.read_csv(f'{data_folder}/{chosen_file}')
    metrics_src.set_index('Node_names', inplace = True)
    
    return metrics_src


def adapt_yeast(graph, data_folder):
    
    
    metrics = Get_metrics_data(graph, data_folder)
    
    
    gabarito = pd.read_csv('../Gabarito.csv')
    
    node_names = np.array(metrics.index)
    protein_names = np.empty(len(node_names)).astype(str)
    
    for i in range(len(node_names)):


        name = node_names[i]
        
        if '.' in name:
            name = name.split('.')[1]
        # print(name)

        
        
        protein_names[i] = name
    
    metrics['Protein_name'] = protein_names
    metrics.set_index('Protein_name', inplace = True)

    gabarito_as_centrality = gabarito.set_index('Gene')
    metrics['gabarito'] = gabarito_as_centrality['Is_essential'].astype(float)

    return metrics, gabarito