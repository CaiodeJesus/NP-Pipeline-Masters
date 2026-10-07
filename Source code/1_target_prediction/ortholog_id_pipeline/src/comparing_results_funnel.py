#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:29:54 2026

@author: caio
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


do_reunion = False

with open('../support/list_of_ids_tritrypdb','r') as f:
    list_of_ids = [linf.strip('\n') for linf in f.readlines()]
 
     
if do_reunion:
    
    dict_steps = {
    'Step1.csv': '1.OrthoMDL',
    'Step2.csv': '2.Reciprocal BLAST',
    'Step3.csv': '3.Synteny'   
    }
    
    results_folder = '../results/common/'
    
    
    print('Total TP->Linf')
    print(len(np.unique(list_of_ids)))
    print()
    
    leish_lists = []
    in_lists = []
    
    df_save = pd.DataFrame()
    df_save['Linf_Targets'] = list_of_ids
    
    for file in sorted(os.listdir(results_folder)):
    
        common = pd.read_csv(results_folder+file)        
        
        print(dict_steps[file])
        print('...reading',file)
        print(len(common),'Linf-Tcruzi pairs')
        for col in common.columns:
            if common.dtypes[col] == object:
                print('...',len(np.unique(common[col])),'unique',col)
        
        
        
        leish_lists.append(np.unique(common['Linf']))
        
        is_in_list = np.empty(len(list_of_ids))
        for i in range(len(list_of_ids)):
            
            # if dict_steps[file] == '2.Reciprocal BLAST':
            #     is_in_list[i] = list_of_ids[i]+'-T1' in np.unique(common['Linf'])
            # else:
            is_in_list[i] = list_of_ids[i] in np.unique(common['Linf'])
        
        in_lists.append(is_in_list)
        df_save[dict_steps[file]] = is_in_list
        print('Common = ',int(is_in_list.sum()),'in the total list')
        print()
         
     
    df_save.to_csv('../results/3_Steps_Together.csv')
    
else:
    
    #Comparing with the original (independent) pipeline
    funnel = pd.read_csv('../results/3_Steps_Together.csv')
    funnel['All'] = funnel['1.OrthoMDL']*funnel['2.Reciprocal BLAST']*funnel['3.Synteny']
    funnel.set_index('Linf_Targets', inplace = True)

    original = pd.read_csv('../../original/results/3_Steps_Together.csv')
    original['All'] = original['1.OrthoMDL']*original['2.Reciprocal BLAST']*original['3.Synteny']
    original.set_index('Linf_Targets', inplace = True)
    
    join = pd.DataFrame()
    join['Linf_Targets'] = funnel.index
    join.set_index('Linf_Targets', inplace = True)
    join['All_funnel'] = funnel['All']
    join['All_original'] = original['All']
    
    join['Both'] = 1 - np.abs(join['All_funnel']-join['All_original'])
    
    
    #Saving results
    funnel.reset_index(inplace = True)
    funnel['All_bool'] = funnel['All'].astype(bool)
    passing_all = funnel.query('All_bool').copy()
    
    list_of_pass_leish = np.unique(passing_all['Linf_Targets'])
    
    print('Passing in all:', len(list_of_pass_leish))
    print()
    
    
    dict_out = {
        'common': [],
        'exclusive': []
                }
    
    for leish in list_of_ids:
        if leish in list_of_pass_leish:
            dict_out['common'].append(leish)
        else:
            dict_out['exclusive'].append(leish)
    
    print('Total number of targets:',len(list_of_ids))
    print('...in common with Trypanosoma cruzi:',len(dict_out['common']))
    print('...exclusive to Leishmania infantum:',len(dict_out['exclusive']))
    
    print()
    for key in dict_out.keys():
        print(key)
        out_file = f'../results/final_leish_{key}_targets_ttdb'
        
        with open(out_file, 'w') as f:
            f.write('\n'.join(dict_out[key]))
            
        print(f'...Saving {key} targets in',out_file) 
