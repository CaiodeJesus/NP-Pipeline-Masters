#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:41:01 2025

@author: caio
"""

from metrics_lib import Get_metrics_data, adapt_yeast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc as AUC
import matplotlib as mpl
import os





def Ranking(metric):
    values = np.array(list(metric))
    ranks = np.zeros(len(metric))


    value_rank = np.unique(list(metric))[::-1] #returns sorted elements
    for i in range(len(value_rank)):
        value = value_rank[i]
        rank = i+1
        is_value = values == value
        # print(value, rank, keys[is_value] )
        
        ranks[is_value] = rank

    return ranks



                 
def ROC_curve(positive, negative, 
              metric, #dataframe with nodes|metric_name
              step = 1,
              only_plot = True):
    
    positive = set(positive)
    negative = set(negative)
    
    metric.sort_values(by = ['rank'], ascending = True, inplace = True)   
    ranked_by_metric =  list(metric.index)
    
    true_pos_rate = []
    false_pos_rate = []
    for i in range(0,len(ranked_by_metric)-1,step):
         
         pred_positive = set(ranked_by_metric[:i+1])
         pred_negative = set(ranked_by_metric[i+1:])
         

         TP = positive & pred_positive
         FP = negative & pred_positive
         FN = positive & pred_negative
         TN = negative & pred_negative
         
         

         tpr = len(TP) / (len(TP) + len(FN))
         fpr = len(FP) / (len(FP) + len(TN))
         
         true_pos_rate.append(tpr)
         false_pos_rate.append(fpr)
    
    if only_plot:
        return true_pos_rate, false_pos_rate
    else:
        return true_pos_rate, false_pos_rate, '|'.join(ranked_by_metric)
         
    
    
def ROC_curves(metrics, gabarito,  
               ax = plt.gca(), title = '',
               fs = 15, do_print = False,
               lw = 2, cmap = 'tab20'):
    
    columns = list(metrics.columns)

    
    gabarito['Is_essential'] = gabarito['Is_essential'].astype(bool)
    essentials = gabarito.query('Is_essential')
    essentials = np.array(essentials['Gene'])
    
    gabarito['Not_essential'] = (1 - gabarito['Is_essential']).astype(bool)
    non_essentials = gabarito.query('Not_essential')
    non_essentials = np.array(non_essentials['Gene'])
    
    
    positive = essentials
    negative = non_essentials
    
    if do_print:
        print()
        print('Random AUC is 0.5, ideal is 1.0')
        
    result = []
    

    n_lines = len(columns)
    colormap = mpl.colormaps[cmap]
    colors = colormap(np.linspace(0, 1, n_lines))

    i = 0
    for col in columns[:]:
        
        metric = pd.DataFrame(metrics[col])
        metric['rank'] = Ranking(metric[col])
        
        TPR, FPR, rank = ROC_curve(positive, negative, metric,
                                   only_plot = False)
        
        ax.plot(FPR, TPR, label = col,
                lw = lw, color = colors[i])
        
        auc = AUC(FPR, TPR)
        if do_print:
            if auc > 0.5:
                print(f'AUC for {col} is {auc:.2f} X')
            else:
                print(f'AUC for {col} is {auc:.2f}')

        # if auc > 0.5:
        result.append([col,rank, auc])
        
        i+=1
        
    ax.plot( [0,1] , [0,1] , 'k--', label = 'random')
    
    ax.set_title(title, size = fs+10)
    
    
    
    ax.tick_params(labelsize = fs-5)
    ax.set_xlabel('False positive rate', fontsize = fs)
    ax.set_ylabel('True positive rate', fontsize = fs)
    
    ax.legend(bbox_to_anchor=(1.0, 1.05), fontsize = fs-5)

    result = pd.DataFrame(result,
                          columns = ['Metric', 'Rank', 'Assertivity(AUC)'])
    result.set_index('Metric', inplace = True)
    
    return result
    
if __name__ == "__main__":

    # graph = 'trpA'
    # metrics = Get_metrics_data(graph)
    
    
    # gabarito = pd.read_csv('../data/Toy_Gabarito_trpA.csv')
    
    
    # ROC_curves(metrics, gabarito, title = 'ROC curves')
    
    
    
  
    
    # data_folder = '../data/subnets_table'
    # subnets = os.listdir(data_folder)
    # graph = subnets[0]

    # metrics, gabarito = adapt_yeast(graph, data_folder)

    # best_metrics = ROC_curves(metrics, gabarito, title = 'ROC curves',
    #                           lw = 3, do_print = True)
    
    
    
    
        
    data_folder = '../metanalysis/chimeras'
    subnets = os.listdir(data_folder)
    graph = 'CL_137_18.12_combined_score.csv'

    metrics = Get_metrics_data(graph, data_folder)
    gabarito = pd.read_csv('../data/Gabarito.csv')
    
    gabarito_as_centrality = gabarito.set_index('Gene')
    metrics['gabarito'] = gabarito_as_centrality['Is_essential'].astype(float)

    best_metrics = ROC_curves(metrics, gabarito, title = 'ROC curves',
                              lw = 3, do_print = True, cmap = 'managua')


    