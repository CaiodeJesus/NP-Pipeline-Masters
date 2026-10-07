#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:30:10 2025

@author: caio
"""

import pandas as pd
import numpy as np
import os




metrics_dict = {'genecards':'',
                'malacards':'',
                'omim':'',
                'disgenet':'_ScoreGDA',
                'ncbi':'',}



threshold_values_dict = {'genecards':0.1,
                         'malacards':-1,
                         'omim':-1,
                         'disgenet':0.06,
                         'ncbi':-1}
    

if threshold_values_dict.keys() != metrics_dict.keys():
    print('Warning, the dictionaries are not equal!')


df_servers = pd.DataFrame()
df_servers['Name'] = threshold_values_dict.keys()
df_servers['Chosen_Metric'] = metrics_dict.values()
df_servers['Threshold_value'] = threshold_values_dict.values()



diseases = ['leish', 'chagas']

uniprot_all = []
size_all = []
for i in range(len(df_servers)):

    # i = 2

    s = df_servers['Name'][i]
    print(s)
    
    size_diseases = []
    uniprot_diseases = []
    for disease in diseases:
    

        print('\t',disease)
        try:
            df = pd.read_csv(f'../data/clean/disease/{disease}/{s}.csv')
        
        except:
            print(f'../data/clean/disease/{disease}/{s}.csv not found')
            uniprot_diseases.append('')
            size_diseases.append(0)
            continue
            
            
        initial_len = len(df)
        df.dropna(subset = ['UniProt'], inplace = True)
        print('\t\t%d without UniProt'%(initial_len - len(df)))
            
            
        p = df_servers['Threshold_value'][i]
            
        metric =  'Metric' + metrics_dict[s]
        print('\t\t'+metric)
        print('\t\tminimum value was %.4f'%df[metric].min())
        print('\t\ta median is',np.median(df[metric]))
            
        print('\t\tcorte de',p)
        df['does_pass'] = df[metric] > p
        df = df.query('does_pass')
        df.pop('does_pass')
            
        print('\tminimal value post filter %.4f'%df[metric].min())
        
        
        unique_unis = np.unique(df['UniProt'])
        uni = '|'.join(unique_unis)
        uniprot_diseases.append(uni)
        
        size_diseases.append(len(unique_unis))

    size_all.append(size_diseases)
    uniprot_all.append(uniprot_diseases)
    
 
size_all_transpose = np.array(size_all).T
uniprot_all_transpose = np.array(uniprot_all).T


for i in range(len(diseases)):
    
    disease = diseases[i]
    df_servers['UniProts_'+disease] = uniprot_all_transpose[i]
    df_servers['Size_'+disease] = size_all_transpose[i]
    

df_servers.to_csv('../data/support/disease_servers.csv')





