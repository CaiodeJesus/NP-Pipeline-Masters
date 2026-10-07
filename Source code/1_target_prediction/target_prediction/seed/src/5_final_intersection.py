#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:29:25 2025

@author: caio
"""

def inter(arr1, arr2):
    if len(arr1)>len(arr2):
        big_arr = arr1
        small_arr = arr2
    else:
        big_arr = arr2
        small_arr = arr1
        
    interseccao = []
    for s in small_arr:
        if s in big_arr:
            interseccao.append(s)

    # print(len(interseccao))
    return np.array(interseccao)
    
    
    
    
    
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn as venn






targets_set_list = [] 

names_sets = ['leish','chagas','compound']

for name in names_sets:
    targets_set = np.loadtxt('../data/support/all_uniprots_'+name+'.txt', dtype = str)
    targets_set = np.unique(targets_set) 
    print('Targets %s:'%name,len(targets_set))
    
    targets_set_list.append(targets_set)



print()
names_plot = []

dict_names={'leish':'Leishmaniasis\nTargets',
            'chagas':'Chagas Disease\nTargets',
            'compound':'Hit Targets'}



center_intersection = inter( inter(targets_set_list[0],
                                  targets_set_list[1]),
                             targets_set_list[2] )

print('Middle overlap', len(center_intersection) )
with open(f'../out/overlaps/middle_{len(center_intersection)}.txt','w') as file:
    file.write('\n'.join(center_intersection))

for target in center_intersection:
    print(target)


for j in range(len(targets_set_list)):
    
    label = dict_names[names_sets[j]]
    names_plot.append(label + '\n(%d)'%len(targets_set_list[j]) )

    for i in range(j):
        
        print(names_sets[i],names_sets[j])
        
        overlap = inter(targets_set_list[i],
                                  targets_set_list[j])
        
        print('\tthe intersection has',len(overlap))
        
        with open(f'../out/overlaps/{names_sets[i]}_{names_sets[j]}_{len(overlap)}.txt','w') as file:
            file.write('\n'.join(overlap))

        exclusive = []
        for target in overlap:
            if target not in center_intersection:
                exclusive.append(target)
        
        exclusive = np.array(exclusive)
       
        print('\twith',len(exclusive),'exclusive targets')
        
        with open(f'../out/overlaps/exclusive_{names_sets[i]}_{names_sets[j]}_{len(exclusive)}.txt','w') as file:
            file.write('\n'.join(exclusive))

   
    



        
print()



fs = 40
L=5
fig = plt.figure(figsize = (4*L,3*L))

plt.rcParams.update({'font.size': fs})


base_color = 'gray'

v = venn.venn3(subsets=[set(targets_set) for targets_set in targets_set_list],
               set_labels = names_plot,
               set_colors = (base_color, base_color, base_color))

for label in ['A','B','C']:
    v.get_label_by_id(label).set_fontweight('bold')
    v.get_label_by_id(label).set_fontsize(fs+7)


keys = list(v.id2idx.keys())

region_dict = {'E_L' : keys[0], #exclusive leish
               'E_C': keys[1], #exclusive chagas
               'E_H': keys[6], #exclusive hit
               'E_LC': keys[2], #exclusive leish chagas
               'E_LH': keys[7], #exclusive leish hit
               'E_CH': keys[8], #exclusive chagas hit
               'M': keys[9], #exclusive leish chagas
               }

green = '#00841dff'
pink = '#ee1787ff'
v.get_patch_by_id(region_dict['E_LH']).set_color(pink)
v.get_patch_by_id(region_dict['E_LH']).set_alpha(1)
v.get_patch_by_id(region_dict['E_LH']).set_edgecolor('black')


v.get_patch_by_id(region_dict['M']).set_color(green)
v.get_patch_by_id(region_dict['M']).set_alpha(1)

v.get_patch_by_id(region_dict['M']).set_edgecolor('black')

      
plt.savefig('../out/venn_diagram.svg')



