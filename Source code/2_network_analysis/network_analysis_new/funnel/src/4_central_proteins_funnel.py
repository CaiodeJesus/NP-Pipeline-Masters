#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 15:01:29 2025

@author: caio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import datetime


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



string_channels = ['neighborhood_on_chromosome',
            'gene_fusion',
            'phylogenetic_cooccurrence',
            'homology',
            'coexpression',
            'experimentally_determined_interaction',
            'database_annotated',
            'automated_textmining',
            'combined_score']


rank_by = 'Average'
n_top = 20 #We will get the top20
# rank_by = ''
# n_top = 10

with open("../input_config") as f:    
    input_lines = f.readlines()

for line in input_lines:
    if 'rank_by' in line:
        # print(line)
        rank_by = line.split(':')[1].strip(' ').strip('\n')
        
    if 'n_top' in line:
        # print(line)
        n_top = line.split(':')[1].strip(' ').strip('\n')
        n_top = int(n_top)


if 'out' not in os.listdir('../'):
    os.mkdir('../out/')
    
if 'files' not in os.listdir('../out/'):
    os.mkdir('../out/files')

today = datetime.datetime.now()

for infile in os.listdir('../data/')[:]:
    if '.' in infile:
        continue
    
    infile_folder = f'../data/{infile}/'
    print(infile)
            
            
    channel = string_channels[-1]
    metrics_table = f'../data/{infile}/subnet_tables/{channel}_graph.csv'
    protein_ranks = pd.read_csv(metrics_table)
    
    protein_ranks.set_index('Node_names', inplace = True)
    metric_names = protein_ranks.columns
    
    protein_ranks['Chimera'] = np.zeros(len(protein_ranks))
    for m in metric_names:
        protein_ranks['Rank_'+m] = Ranking(protein_ranks[m])
        protein_ranks['Chimera'] -= protein_ranks['Rank_'+m]
        
    protein_ranks['Chimera'] = Ranking(protein_ranks['Chimera'])
      
    
    all_channels = pd.DataFrame()
    all_channels['Node_names'] = protein_ranks.index
    all_channels.set_index('Node_names', inplace = True)
    all_channels[channel] = protein_ranks['Chimera']
        
    
    for channel in string_channels[:-1][:]:
        
        try:
            metrics_table = f'../data/{infile}/subnet_tables/{channel}_graph.csv'
            protein_ranks = pd.read_csv(metrics_table)
            
        except:
            continue
        
        protein_ranks.set_index('Node_names', inplace = True)
        metric_names = protein_ranks.columns
        
        protein_ranks['Chimera'] = np.zeros(len(protein_ranks))
        for m in metric_names:
            protein_ranks['Rank_'+m] = Ranking(protein_ranks[m])
            protein_ranks['Chimera'] -= protein_ranks['Rank_'+m]
            
        protein_ranks['Chimera'] = Ranking(protein_ranks['Chimera'])
          
        
        all_channels[channel] = protein_ranks['Chimera'].fillna(0)
        
    
    all_channels.fillna(0, inplace = True)
    
    
    
    average = np.zeros(len(all_channels))
    for channel in all_channels.columns:
        average -= all_channels[channel]
    
    average = Ranking(average)
    all_channels['average'] = average
    
    ##########################################################################
       
    #5.1 Guardando os rankings de cada canal
    
    
    # dict_uniprot = {}
    if 'ranks' not in os.listdir(infile_folder):
        os.mkdir(infile_folder+'ranks')
        
           
    for channel in string_channels[:]:
        channel_rank = pd.DataFrame(index = all_channels.index)
        # channel_rank['StringID'] = all_channels['Node_names']
        
        try:
            channel_rank[channel] = all_channels[channel]
        except:
            continue
        
        channel_rank.sort_values(by = channel, inplace = True)
        channel_rank.reset_index(inplace = True)

        
        dict_uniprot = {}
        with open(f'{infile_folder}ranks/{channel}.txt','w') as file:
            for i in range(len(channel_rank)):
                line = channel_rank.iloc[i]
                
                name = line['Node_names']
                # string_id = line['StringID']
                rank = line[channel]
                
                uniprot = name.split('.')[1]
                
                dict_uniprot[name] = uniprot
                
                # line['Node_names'] = uniprot
                file.write(f"{rank}º\t{name}\t")
                file.write(f"uniprot.org/uniprotkb/{uniprot}/entry\n")
                
    
    
    
    
    
    # all_channels.pop('UniProt')
    all_channels.to_csv(f'{infile_folder}all_rankings.csv')









    num_prots = len(all_channels)
    
    
    all_channels[f'is_top{n_top}'] = all_channels[rank_by] < n_top+1

    top_n = all_channels.query(f'is_top{n_top}').copy()
    top_n.sort_values(by = rank_by, inplace = True)

    top_n.pop(f'is_top{n_top}')
    top_n['UniProt'] = dict_uniprot




    fig, ax_list = prep_plot(1, 1,
                             scale_fig = 1,
                             ratio=[15,5])
    
    
    fs = 20
       
    cmap = mpl.colormaps['hsv']
    
    colors = cmap(np.linspace(0,1,len(top_n)))
    
    #cores correspondendo a cada um dos top genes    
    
    # ##########################################################################
    
    top_n.reset_index(inplace = True)
    #5.1 montando o plot que mostra o ranking das métricas ao longo dos testes
    
    text_top_n = ''
    for i in range(len(top_n)): #plotar ao contrário pros de maior pontuação ficarem na frente?
        
        line = top_n.iloc[i]
        
        up = line['UniProt']
        n = line['Node_names']
        av = line['average']
        
        rank = line[rank_by]
        
        line.pop('UniProt')
        line.pop('Node_names')
        line.pop('average')
        
        line = np.array(line)
        
        mean = line.mean()

        plot_line = list(line)+[mean, av]
        
        text_top_n += f"{rank}º\t{n}\tUniProt:{up}\n"
        
        
        ax_list[0].plot(np.arange(len(line)+2), plot_line, 'o--', lw = 2, ms = 4, color = colors[i])
        ax_list[0].plot([len(line)+1],[av],'s',ms=10, color = colors[i], label = up)
        

    
    ax_list[0].set_ylim(n_top+1,-1)
    
    xticks = list(top_n.columns)[1:-2] + ['','Average']
       
    ax_list[0].set_xticks(np.arange(len(line)+2),xticks,
                 rotation = 45, rotation_mode = 'anchor',
                 horizontalalignment='right', verticalalignment='top')
    
    ax_list[0].set_ylabel('Ranking position', fontsize = fs)
    
    ax_list[0].legend(bbox_to_anchor=(1.0, 1.05), fontsize = fs-5)
    
    title = infile
    title += f"\nTop{n_top} of {num_prots} by {rank_by.replace('_',' ').title()} ranking"
    ax_list[0].set_title(title, fontsize = fs)
    ax_list[0].tick_params(labelsize = fs-5)    
    
    plt.savefig(f'../data/{infile}/rankings.svg',bbox_inches="tight")
    plt.clf()


    print(f'Top{n_top} of {num_prots} by {rank_by.replace('_',' ').title()} ranking')
    print(text_top_n)
    print()
    print()
    
    
    with open(f'../out/infiles/{infile}.txt', 'w') as f:
        f.write('#'+today.strftime("%c")+'\n')
        f.write(f'#Top{n_top} of {num_prots} by {rank_by.replace('_',' ').title()} ranking')
        f.write('\n')
        f.write(text_top_n)
        
    








