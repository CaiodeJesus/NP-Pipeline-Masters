#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:39:03 2025

@author: caio
"""

import networkx as nx
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string_api



    

string_channels = {'neighborhood_on_chromosome':'nscore',
                   'gene_fusion':'fscore',
                   'phylogenetic_cooccurrence':'pscore',
                   'homology':'hscore',
                   'coexpression':'ascore',
                   'experimentally_determined_interaction':'escore',
                   'database_annotated':'dscore',
                   'automated_textmining':'tscore',
                   'combined_score':'score',
                   'node_1':'preferredName_A',
                   'node_2':'preferredName_B',
                   'id1':'stringId_A',
                   'id2':'stringId_B',
                   }


for infile in os.listdir('../data/')[:]:
    
    if '.' in infile:
        continue

    src_folder = f'../data/{infile}/'
    print(infile)

    raw = pd.read_csv(src_folder+'interactions.csv')


    if 'channels' not in os.listdir(src_folder):
        os.mkdir(src_folder+'channels')
    

    for channel in list(string_channels.keys())[:-4][:]:
        # channel = 'automated_textmining'
        print(infile,channel)
    
      
        columns_to_copy = ['node_1','node_2',
                           'id1', 'id2',
                           channel]
        
        df_channel = pd.DataFrame()
        for col in columns_to_copy:
            df_channel[col] = raw[string_channels[col]]
            
        df_channel['has_channel'] = df_channel[channel] != 0
        df_channel = df_channel.query('has_channel')
        df_channel.pop('has_channel')
        
        
        list_of_nodes = list(df_channel['id1']) + list(df_channel['id2'])
        list_of_nodes = np.unique(list_of_nodes)
        
        
        num_edges = len(df_channel)
        print(infile,channel,len(list_of_nodes), num_edges)
        
        
        if not len(list_of_nodes): #if it is zero
            continue
        
        
        
        G = nx.Graph()
        for i in range(len(df_channel)):
            line = df_channel.iloc[i]
            n1, n2 = line['id1'], line['id2']
            weight = line[channel]

            G.add_edge(n1, n2, weight=weight)
        
        if ( len(list_of_nodes) == len(G.nodes) ) and ( num_edges == len(G.edges) ):
            print(infile,channel,'Everything ok!')
        else:
            print('################Something is wrong.')
        
        
        all_string_ids = list(df_channel['id1']) + list(df_channel['id2'])
        all_string_ids = np.unique(all_string_ids)
        num_nodes = len(all_string_ids)
        
        out_name =  f'../data/{infile}/channels/{num_nodes}_{channel}_'
        #Getting the STRING network image
        
        
        if 'LEISH' in infile.upper():
            org_id = 5671
            # print('\t...Leishmania infantum PPI')
        elif 'TRYP' in infile.upper():
            org_id = 353153
            # print('\t...Trypanosoma cruzi PPI')
        else:
            org_id = 9606
            # print('\t...Searching Human PPI')
            
            
            
        string_api.Network_Image(search_ids = all_string_ids,
                                 call_id = 'caio_lqmc',
                                 species_id = org_id,
                                 folder_img = out_name)
        
        df_channel.to_csv(out_name+'table.csv')
        
        nx.write_graphml_lxml(G, out_name + 'graph.graphml')
        
    
    
#To do list:
#-Colormap regions to color edges
#-Ver diferença exclusive e outras regiões
# 	tipo, mesmas coordenadas do chagas_compound e as conexões do exclusive chagas comp
#-Ver as que ficaram pra trás
#-Desenhar um modo que dê pra ver cada canal, que nem eu fiz antes com as coordenadas
#Plotar as redes com o Nx de uma maneira legal
    # "STRING uses a spring model to generate the network images.
    # Nodes are modeled as masses and edges as springs;
    # the final position of the nodes in the image
    # is computed by minimizing the 'energy' of the system.
    # We give high confidence edges a higher 'spring strength'
    # so that they will reach an optimal position before lower confidence edges.
    # The user also can optionally reduce the 'natural length' of
    # a high confidence edge - this forces them closer together and
    # sometimes results in a clearer picture of high confidence interactions.
    # We set the high confidence edge length to 80% of the normal length by default.!
    # https://string-db.org/help/getting_started/






