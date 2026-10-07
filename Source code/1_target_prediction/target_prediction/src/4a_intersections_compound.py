#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:35:05 2025

@author: caio
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm


df_servers = pd.read_csv('../data/support/compound_servers.csv')



uniprot_intersection_matrix = np.empty(shape=(len(df_servers),len(df_servers))).astype(str)
len_intersection_matrix = -1*np.ones(shape=(len(df_servers),len(df_servers)))
norm_len_intersection_matrix = len_intersection_matrix.copy()

i,j = 3,3

for j in range(len(df_servers)):
    for i in range(j+1):
            
        try:
            uni_i = df_servers['UniProts'][i].split('|')
            uni_j = df_servers['UniProts'][j].split('|')
            
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


yticks = []
for i in range(len(df_servers)):
    line = df_servers.iloc[i]
    
    yticks.append(f"{line['Name']} ({line['Size']})")


plot = plt.imshow(norm_len_intersection_matrix[::-1], cmap = 'inferno')
plt.grid(False) #pq tava dando erro
plt.colorbar(plot)

plt.tight_layout()
plt.yticks(np.arange(len(df_servers)),
           yticks[::-1])
plt.xticks(np.arange(len(df_servers)),
           list(df_servers['Name']), rotation = 'vertical')

plt.savefig('../data/support/compound_intersections.svg', bbox_inches='tight')

df_servers = df_servers.dropna(subset = ['UniProts'])
all_uniprots = '|'.join(df_servers['UniProts'])
all_uniprots = all_uniprots.split('|')
all_uniprots = np.unique(all_uniprots)



all_uniprots_to_save = '\n'.join(all_uniprots)

with open('../data/support/all_uniprots_compound.txt','w') as file:
    file.write(all_uniprots_to_save)
    
norm_len_intersection_matrix = norm_len_intersection_matrix[::-1]  
len_intersection_matrix = len_intersection_matrix[::-1]
