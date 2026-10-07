#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 13:36:41 2025

@author: caio
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from assertivity_test import Ranking
import matplotlib as mpl
from metrics_lib import prep_plot

if __name__ == "__main__":

    # print(os.listdir('../../gabarito/subnetwork_division_v0.csv'))

#1  #pegando a lista de clusters

    #lendo os resultados dos testes vindos de best_centralies.py
    data_folder = '../data/results/tabela/'
    # graph = subnets[0]
    graph = 'chimera_combined_score_results.csv'
    print(graph)
    print()

    final_results = pd.read_csv(f"../string_channels/results/{graph}")
    final_results = final_results.iloc[:]
    
    final_results['Destruction(AUC)'] = final_results['Destruction(AUC)'] / (6000*100/2)
    #Dividing by (network_size)*100 because that is the size of the Destruction plot square,
    #so that we get a percentage and its easier to compare
    graph = graph.replace('.csv','')
    
    
# ###############################################################################
#2  #Montando um dataframe rankeando as métricas
    #de acordo com o desempenho de cada uma nos testes
    rankings = pd.DataFrame()
    rankings['Metric'] = final_results['Metric']
    
    
    dict_order = {'Assertivity(AUC)':1, #quanto maior, melhor
                  'Contribution(sum of PCA component)':1, #quanto maior, melhor
                  'Correlation (Sum of |Tau|)':-1, #quanto menor, melhor
                  'Destruction(AUC)':-1} #quanto menor, melhor
    
    rankings['Consensus'] = np.zeros(len(rankings))

    #criando um consenso (uma visão geral dos rankings ao longo dos testes)
    for col in final_results.columns[2:]:
        # print(f"'{col}',")
        rankings[col+'_rank'] = Ranking(dict_order[col]*final_results[col])
        rankings['Consensus']+= rankings[col+'_rank']
    
    
    rankings['Consensus_rank'] = Ranking(-1*rankings['Consensus'])

# ###############################################################################
#3  #Criando um dataframe protein_rank com os rankings das proteinas
    #de acordo com cada metrica
    metric = 'c1'
    final_results.set_index('Metric', inplace = True)
    
    rank = final_results.loc[metric]['Rank']
    rank = rank.split('|')
    

    protein_ranks = pd.DataFrame()
    protein_ranks['Node'] = rank
    protein_ranks['Ranking_'+metric] = np.arange(len(rank)) + 1
    
    protein_ranks.set_index('Node', inplace = True)
    
    #um sistema de pontos pro ranking é legal, pq visualmente fica melhor
    #quanto maior a barra, mais alta a posição do nó
    protein_ranks['Consensus_pontos'] = np.zeros(len(protein_ranks))
    for metric in final_results.index[:]:
        rank = final_results.loc[metric]['Rank']
        rank = rank.split('|')
        

        rank_metric = pd.DataFrame()
        rank_metric['Node'] = rank
        rank_metric['Ranking'] = np.arange(len(rank)) + 1
        
        rank_metric.set_index('Node', inplace = True)
        
        protein_ranks[metric] = rank_metric['Ranking']
        
        
        protein_ranks[metric+'_pontos'] = len(rank_metric) - rank_metric['Ranking']
        
        protein_ranks['Consensus_pontos'] += protein_ranks[metric+'_pontos']

    
    
    protein_ranks['Consensus'] = Ranking(1*protein_ranks['Consensus_pontos'])
    
###############################################################################
#4  #Montando as chimeras
    #(pontuação do consenso até o nº lugar no ranking dos testes)


    chimeras = pd.DataFrame()
    chimeras['Node_names'] = protein_ranks.index
    
    chimeras.set_index('Node_names', inplace = True)

    max_rank = 2
    
    dict_chimeras = {}
    for max_rank in range(int(rankings['Consensus_rank'].max()) + 1):

        rankings[f'is_best_{max_rank}'] = rankings['Consensus_rank'] <= max_rank
        
        best_max_rank = rankings.query(f'is_best_{max_rank}')
        best_max_rank = np.array(best_max_rank['Metric'])
        best_max_rank = best_max_rank[best_max_rank!='gabarito'] #tirando o gabarito
        print(max_rank, ','.join(best_max_rank))
        
        if len(best_max_rank) == 0:
            continue
        
        chimeras[f'c{max_rank}'] = np.zeros(len(chimeras))
    
        for m in best_max_rank:
            chimeras[f'c{max_rank}'] += protein_ranks[m+'_pontos']

        max_points = chimeras[f'c{max_rank}'].max()
        min_points = chimeras[f'c{max_rank}'].min()
        
        chimeras[f'c{max_rank}'] = ( chimeras[f'c{max_rank}'] - min_points) / (max_points - min_points)    
        
        
        
        
        
        dict_chimeras[f'c{max_rank}'] = ','.join(best_max_rank)
        
    
    chimeras.to_csv(f'../string_channels/results/chimeras_{graph}.csv')

    
    
    
# ###############################################################################
#5  #montando a prancha que vai mostrar os resultados visualmente     
        
    
        
    fig, ax_list = prep_plot(1, 1,
                             scale_fig = 1,
                             ratio=[15,5])

    
    fs = 20
   
    names_metrics = list(final_results.index)
    n_lines = len(names_metrics)
    colors = mpl.colormaps['tab20']   #cores correspondendo a cada métrica     
    colors = colors(np.linspace(0, 1, n_lines))
    ##########################################################################
   
    #5.1 montando o plot que mostra o ranking das métricas ao longo dos testes
    
    dict_colors = {}
    for i in range(len(final_results)):
        name = names_metrics[i]
        print(name)
        dict_colors[name] = colors[i]
    
    
    rankings = rankings.sort_values(by = ['Consensus_rank'])
    for i in range(len(rankings)):
        
        line = rankings.iloc[i]
        at = line['Assertivity(AUC)_rank']
        cb = line['Contribution(sum of PCA component)_rank']
        cr = line['Correlation (Sum of |Tau|)_rank']
        dt = line['Destruction(AUC)_rank']

        media = (at+cb+cr+dt)/4
        con = line['Consensus_rank']

        # plot_line = [at,cb,cr,dt,media, con]
        plot_line = [at,cb,cr,dt,con]
        
        name = line['Metric']
        ax_list[0].plot([1,2,3,4,5], plot_line, 'o--', lw = 1, ms = 4, color = dict_colors[name])#, label = line['Metric'])
        # ax.plot([5],[media],'*',ms=10, color = colors[i])#, label = line['Metric'])
        ax_list[0].plot([5],[con],'s',ms=12, color = dict_colors[name], label = f'{con}° - {name}')
        
        
    ax_list[0].set_ylim(len(rankings),0)
    
    xticks = ['Assertivity',#'\n(AUC)',
              'Contribution',#\n(sum of PCA component)',
              'Correlation',#\n(Sum of |Tau|)',
              'Destruction',#\n(AUC)',
              'Mean Rank']
   
    ax_list[0].set_xticks([1,2,3,4,5],xticks,
                 rotation = 45, rotation_mode = 'anchor',
                 horizontalalignment='right', verticalalignment='top')
    
    ax_list[0].set_ylabel('Ranking position', fontsize = fs)
    
    ax_list[0].legend(bbox_to_anchor=(1.0, 1.05), fontsize = fs-5)
    
    # ax_list[0].set_title(graph, fontsize = fs)
    ax_list[0].tick_params(labelsize = fs-5)    
    
    
    plt.savefig('../Chimera_ranking.svg', bbox_inches="tight")
