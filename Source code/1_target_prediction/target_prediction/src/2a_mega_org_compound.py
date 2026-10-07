#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:35:05 2025

@author: caio
"""

import pandas as pd
import numpy as np
import os

metrics_dict = {'stitch':'',
                'swiss':'',
                'pharmapper':'_Norm Fit',
                'target_net':'',
                'superpred':'_Probability',
                'sea':'_neg_logPv',
                'chemmapper':'_score'}



threshold_values_dict = {'stitch':0.16,
                         'swiss':0.0,
                         'pharmapper':0.9,
                         'target_net':0.5,
                         'superpred':50.0,
                         'sea':5,
                         'chemmapper':0.0}


if threshold_values_dict.keys() != metrics_dict.keys():
    print('Warning, dictionaries are not equal!')


df_servers = pd.DataFrame()
df_servers['Name'] = threshold_values_dict.keys()
df_servers['Chosen_Metric'] = metrics_dict.values()
df_servers['Threshold_value'] = threshold_values_dict.values()


sizes_of_result = []
uniprots = []

for i in range(len(df_servers)):
    # i = 6
    
    s = df_servers['Name'][i]
    print(s)
    df = pd.read_csv('../data/clean/compound/'+s+'.csv')
    
    
    initial_len = len(df)
    df.dropna(subset = ['UniProt'], inplace = True)
    print('\t%d without UniProt'%(initial_len - len(df)))
    
    
    p = df_servers['Threshold_value'][i]
    
    metric =  'Metric' + metrics_dict[s]
    # print(metric)
    print('\tminimum value was %.4f'%df[metric].min())
    
    print('\tcorte de',p)
    df['does_pass'] = df[metric] > p
    df = df.query('does_pass')
    df.pop('does_pass')
    
    print('\tminimum value post threshold is %.4f'%df[metric].min())


    unique_unis = np.unique(df['UniProt'])
    uni = '|'.join(unique_unis)
    uniprots.append(uni)

    sizes_of_result.append(len(unique_unis))


df_servers['Size'] = sizes_of_result
df_servers['UniProts'] = uniprots


df_servers.sort_values(by=['Size'], ascending = True, inplace=True)
df_servers.reset_index(inplace = True)

df_servers.to_csv('../data/support/compound_servers.csv')













