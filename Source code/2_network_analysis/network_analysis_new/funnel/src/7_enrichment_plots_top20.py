#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:32:15 2026

@author: caio
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_colors(inp, colormap, vmin=None, vmax=None):
    
#    https://stackoverflow.com/questions/28144142/how-can-i-generate-a-colormap-array-from-a-simple-array-in-matplotlib
 
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp)), norm

def size_marker(array):
    
    x = np.array(array)
    a = 20
    b = 100
    mask = x<0
    x[mask] = -b
    y = a*x + b
    return y


region = 'exclusive'

src_all_file = f'../enrichment/string/enrichment_{region}.all.tsv'
src_top_file = f'../enrichment/string/enrichment_{region}top.all.tsv'


src_all = pd.read_csv(src_all_file, sep = '\t')
src_top = pd.read_csv(src_top_file, sep = '\t')

dict_gsets = {}
gsets = [
    'GO Component',
    'GO Function',
    'GO Process',
    'KEGG'
    ]





for gset in gsets[:]:
    src_all['is_gset'] = src_all['#category'] == gset
    df_all = src_all.query('is_gset').copy()
    # number_of_ids = dict_gsets[gset+'_all']['term ID']
    # print(len(number_of_ids), len(np.unique(number_of_ids)))
    df_all.set_index('term ID', inplace = True)
    

    src_top['is_gset'] = src_top['#category'] == gset
    df_top = src_top.query('is_gset').copy()
    # number_of_ids = dict_gsets[gset+'_top']['term ID']
    # print(len(number_of_ids), len(np.unique(number_of_ids)))
    df_top.set_index('term ID', inplace = True)
    # dict_gsets[gset+'_top'] = df_top.copy()

    
    for col in ['observed gene count','signal','false discovery rate']:
        df_all[f'{col}_top'] = df_top[col]
    
    df_all[r'$-log10(FDR)$'] = -np.log10(df_all['false discovery rate'])
    df_all['minus_log_fdr_top'] = -np.log10(df_all['false discovery rate_top'])
    
    
    df_all.fillna(-1, inplace = True)
    # dict_gsets[gset] = df_all
    dict_gsets[gset] = df_all.copy()

# [1.8, 6, 0.1, 6.5]
coord_cbar_dict = {
    'all_GO Component': [1.5, 4, 0.1, 6.5], 
    'all_GO Function': [5., 4, 0.2, 6.5],
    'all_GO Process': [5., 4, 0.2, 6.5],
    'all_KEGG': [2, 4, 0.05, 6.5],
    
    'common_GO Component': [1.7,4,0.07,6.5],
    'common_GO Function': [3.5,4,0.15,6.5],
    'common_GO Process': [3.2,4,0.15,6.5],
    'common_KEGG': [1.9,4,0.1,6.5],
    
    
    'exclusive_GO Component': [1.15,0,0.05,2.5],
    'exclusive_GO Function': [3.4,5,0.1,5],
    'exclusive_GO Process': [3.4,5,0.1,5],
    'exclusive_KEGG': [0.6,-0.3,0.025,0.3],
        }

# (24, 25)
figsize_dict = {
    'all_GO Component': (24, 25),
    'all_GO Function':  (24, 25 ),
    'all_GO Process':   ( 24, 25 ),
    'all_KEGG':         ( 24, 25 ),
    'common_GO Component': ( 24, 25 ),
    'common_GO Function':  ( 24, 25 ),
    'common_GO Process':   ( 24, 25 ),
    'common_KEGG':         ( 24, 25 ),
    'exclusive_GO Component': ( 24, 15 ),
    'exclusive_GO Function':  ( 24, 25 ),
    'exclusive_GO Process':   ( 24, 25 ),
    'exclusive_KEGG':         ( 10, 10 ),
        }

titles_dict = {
    'all_GO Component': 'Componentes Celulares\nAlvos Totais', 
    'all_GO Function': 'Funções Moleculares\nAlvos Totais',
    'all_GO Process': 'Processos Biológicos\nAlvos Totais',
    'all_KEGG': 'Vias KEGG\nAlvos Totais',
    
    'common_GO Component': 'Componentes Celulares\nAlvos Comuns',
    'common_GO Function': 'Funções Moleculares\nAlvos Comuns',
    'common_GO Process': 'Processos Biológicos\nAlvos Comuns',
    'common_KEGG': 'Vias KEGG\nAlvos Comuns',
    
    
    'exclusive_GO Component': 'Componentes Celulares\nAlvos Exclusivos',
    'exclusive_GO Function': 'Funções Moleculares\nAlvos Exclusivos',
    'exclusive_GO Process': 'Processos Biológicos\nAlvos Exclusivos',
    'exclusive_KEGG': 'Vias KEGG\nAlvos Exclusivos',
        }

do_save = True
N = 3

####################
# gset = gsets[3]



# for gset in gsets[:]:

for gset in gsets[N:N+1]:
# if True: #just a placeholder for the identation
    df = dict_gsets[gset]
    
    #x, y, width, height

    
    variables = {
        'xaxis':'signal',
        'marker_size':'observed gene count',
        'color':r'$-log10(FDR)$'
        }
    
    
    
    width, height = figsize_dict[region+'_'+gset]
    fs = 30
    print(region+'_'+gset)
    
    
    fig = plt.figure(figsize = (width, height) )    
    ax = fig.add_subplot(1,1,1)
    
    
    
    df['is_significant'] = df['false discovery rate'] < 10**-1
    df = df.query('is_significant')
    df.sort_values(by=[variables['xaxis']],
                   inplace = True, ascending = False)
    
    df = df.iloc[:20]
    
    
    df.sort_values(by=[variables['xaxis']],
                   inplace = True, ascending = True)

    #Ylabels
    ylabels = df['term description']
    ypos = np.arange(len(ylabels))
    
    cmap = plt.cm.cool
    
    array_to_colors = df[variables['color']]
    colors, norm = get_colors(array_to_colors, cmap)
    
    color_top = 'black'
    
    
    #Lines
    
    if region+'_'+gset == 'exclusive_GO Component':
        bar_plot = ax.barh(ypos,
                           df[variables['xaxis']],
                           height = 0.08, #grossura das barras
                           color = colors,
                           )
    elif region+'_'+gset == 'exclusive_KEGG':
        bar_plot = ax.barh(ypos,
                           df[variables['xaxis']],
                           height = 0.02, #grossura das barras
                           color = colors,
                           )
    else:
        bar_plot = ax.barh(ypos,
                           df[variables['xaxis']],
                           height = 0.2, #grossura das barras
                           color = colors,
                           )
    
    
    #Markers
    array_to_size = df[variables['marker_size']]
    
    size = size_marker(array_to_size)
    
    xplot = np.array(df[variables['xaxis']])
    ax.scatter(xplot, ypos,
               s = size , color = colors)
    
    
    #Mostrando os que se conservaram nos top
    size_top = size_marker(df[variables['marker_size']+'_top']) 
    mask = size_top > 0
    
    
    ax.scatter(xplot[mask], ypos[mask],
               s = size_top[mask],
               color = color_top)
    
    bar_plot = ax.barh(ypos[mask], xplot[mask],
            height = 0.1, #grossura das barras
            color = color_top,
            )
    
    
    yticks = []
    text = ''
    for ylabel in ylabels:
        
        # print(ylabel)
        # ylabel_res = ylabel[:40]
        if len(ylabel) > 80:
            print(ylabel,'é maior q 80')
            if ylabel[len(ylabel)//2 - 1] != ' ':
                ylabel_res = ylabel[:len(ylabel)//2] + '\n' + ylabel[len(ylabel)//2:]
            else:
                ylabel_res = ylabel[:len(ylabel)//2] + '-\n' + ylabel[len(ylabel)//2:]  
        elif len(ylabel) > 40:
            # print('maior q 40')
            words = ylabel.split(' ')
            # print(words)
            ylabel_res = ''
            is_done = False
            for n in range(len(words)):
                # print()
                word = words[n]
                if not is_done:
                    if len(ylabel_res+word)+1 < 40:
                        ylabel_res = ylabel_res + word + ' '
                        # print(ylabel_res)
                    else:
                        ylabel_res = ylabel_res + '\n' + ' '.join(words[n:])
                        # print(ylabel_res)
                        is_done = True
        else:
            ylabel_res = ylabel
        
        # print()
        yticks.append(ylabel_res.strip(' '))
        
    
    
    ax.set_yticks(ypos, yticks, size = fs)
    ax.tick_params(axis = 'x', labelsize=fs)
    
    if region+'_'+gset == 'exclusive_GO Component' :
        ax.set_title(titles_dict[region+'_'+gset].split('\n')[1], fontsize = fs+10, y = 1.0)
    elif region+'_'+gset == 'exclusive_KEGG' :
        ax.set_title(titles_dict[region+'_'+gset].split('\n')[1], fontsize = fs+10, y = 0.97)
      
    else:
        ax.set_title(titles_dict[region+'_'+gset].split('\n')[1], fontsize = fs+10, y = 1.05)
        
    fig.suptitle(titles_dict[region+'_'+gset].split('\n')[0], fontsize = fs+15, fontweight = 'bold')


    ax.set_xlabel('Sinal de Enriquecimento', fontsize = fs+5)
    ax.grid(axis='x')
    
    
    
    
    # n_legend = 4
    

    array_size_legend = np.linspace(10,
                                    array_to_size.max(),
                                    4)

    array_size_legend = np.array([np.ceil(i/10)*10 for i in array_size_legend])
    array_size_legend = np.unique(array_size_legend)
    
    if len(array_size_legend) < 4:
        array_size_legend = np.linspace(5,
                                        array_to_size.max(),
                                        4)

        array_size_legend = np.array([np.ceil(i/5)*5 for i in array_size_legend])
        array_size_legend = np.unique(array_size_legend)
        
    
    
    size_legend = size_marker(array_size_legend)
    
    for i in range(len(array_size_legend)):
        ax.scatter(-i-2,-i-2,
                   s = size_legend[i],
                   label = int(array_size_legend[i]),
                   color = 'k')
    
    ax.set_ylim(0-0.5,len(ylabels)-1+0.5)
    ax.set_xlim(0)
    ax.legend(loc = 'upper right',  
              bbox_to_anchor=(1.35, 1),
              #https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
              title = 'Contagem de genes',
              frameon = True,
              title_fontsize = fs+5,
              labelspacing = 2,
              fontsize = fs)
    
    
    coord_axes_cbar = coord_cbar_dict[region+'_'+gset]
    
    cax = ax.inset_axes(coord_axes_cbar,
                        transform=ax.transData)
    
    cb = mpl.colorbar.ColorbarBase(cax,
                                   orientation='vertical', 
                                   cmap = cmap,
                                   norm=norm,  # vmax and vmin
                                   extend='both',
                                   label='This is a label',
                                   )
    # #https://stackoverflow.com/questions/16595138/standalone-colorbar
    
    cax.set_ylabel(variables['color'],
                   size = fs+10,
                   labelpad = 10,
                   loc = 'center'
                   )
    
    cax.tick_params(labelsize=fs)
    # plt.suptitle(region, fontsize = fs)
    print(gset)
    # print(ylabels)
    
    
    
    text = ''
    indices = list(ylabels.index)
    termos = list(ylabels)
    
    for i in range(len(ylabels)):
        text += termos[i]+'('+indices[i]+')' + '; '
    
    text += '\n Se conservaram nos top20 \n'


    indices = list(ylabels.iloc[ypos[mask]].index)
    termos = list(ylabels.iloc[ypos[mask]])
    for i in range(len(ylabels.iloc[ypos[mask]])):
        text += termos[i]+'('+indices[i]+')' + '; '
    
    
    
    print(text)
    print()    
    
    # print('-Width:',width)
    # print('-Height:',height)
    # print('Axes cbar:', coord_axes_cbar)
    
    if do_save:
        plt.savefig(f'../enrichment/plots/top20/svgs/{region}_{gset}.svg',
                    bbox_inches = 'tight')
        plt.savefig(f'../enrichment/plots/top20/{region}_{gset}.png',
                    bbox_inches = "tight")

    # plt.clf()


#########################################################
# ax.stem(df[variables['xaxis']],
#         orientation = 'horizontal',
#         basefmt = '')

#To do:
#-colocar alguns ifs no caso de alguns values do dicionário variables seja "", pra entender q n é pra plotar aquela parte
#-e rotacionar as variáveis nas chaves do dicionário variables tbm
#-formula pra posição da barra de acordo com o width e height
###por exemplo:
# # h_cbar = np.ceil(height/2)
# # coord_axes_cbar = [5, #x
# #                    0, #y
# #                    h_cbar/40, #width
# #                    h_cbar] #height
