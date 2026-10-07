#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 17:33:56 2025

@author: caio
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm


df_servers = pd.read_csv('../data/support/disease_servers.csv')
print(len(df_servers))

diseases = ['leish','chagas']

for disease in diseases:
    uniprot_intersection_matrix = np.empty(shape=(len(df_servers),len(df_servers))).astype(str)
    len_intersection_matrix = -1*np.ones(shape=(len(df_servers),len(df_servers)))
    norm_len_intersection_matrix = len_intersection_matrix.copy()
    
    # i,j = 3,3
    print(len(df_servers))
    for j in range(len(df_servers)):
        for i in range(j+1):
        
            
            try:
                uni_i = df_servers['UniProts_' + disease][i].split('|')
                uni_j = df_servers['UniProts_' + disease][j].split('|')
                
                inter = np.intersect1d(uni_i, uni_j)
                
                uniprot_inter = '|'.join(inter)
                
                len_inter = len(inter)
                
                norm_len_inter = len_inter / len(uni_i) 
                #quantos % de i tão em j, pq o i é menor q j
                
            
            
            except:
                len_inter = 0
                uniprot_inter = ''
                norm_len_inter = -1
    
            uniprot_intersection_matrix[j][i] = uniprot_inter
            len_intersection_matrix[j][i] = len_inter
            norm_len_intersection_matrix[j][i] = norm_len_inter*100
    
            # print(len_inter)
            
    
    
    
    plot = plt.imshow(norm_len_intersection_matrix[::-1], cmap = 'inferno')
    plt.colorbar(plot)
    
    
    yticks = []
    for i in range(len(df_servers)):
        line = df_servers.iloc[i]
        
        yticks.append(f"{line['Name']} ({line['Size_'+disease]})")
        
    plt.yticks(np.arange(len(df_servers)),
                yticks[::-1])
    plt.xticks(np.arange(len(df_servers)),
                list(df_servers['Name']), rotation = 'vertical')
    
    plt.title(disease)
    plt.tight_layout()
    plt.savefig('../data/support/disease_'+disease+'_intersections.svg', bbox_inches='tight')
    plt.clf()
    
    
    df_uniprots = df_servers.dropna( subset = ['UniProts_'+disease])#, inplace = True)
    all_uniprots = '|'.join(df_uniprots['UniProts_'+disease])
    all_uniprots = all_uniprots.split('|')
    all_uniprots = np.unique(all_uniprots)
    
    all_uniprots_to_save = '\n'.join(all_uniprots)
    
    print(disease)
    with open('../data/support/all_uniprots_'+disease+'.txt','w') as file:
        file.write(all_uniprots_to_save)