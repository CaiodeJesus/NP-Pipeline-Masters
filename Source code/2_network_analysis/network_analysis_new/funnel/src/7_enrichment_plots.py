#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:54:35 2026

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
    a = 5
    b = 5
    mask = x<0
    x[mask] = -b
    y = a*x +b
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





for gset in gsets:
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



coord_cbar_dict = {
    'all_GO Component': [1.7, 2, 0.15, 6.5],
    'all_GO Function': [5, 38, 0.2, 16],
    'all_GO Process': [4.8, 60, 0.15, 16],
    'all_KEGG': [2.2, 0, 0.1, 12],
    'common_GO Component': [1.8, 7, 0.1, 6.5],
    'common_GO Function': [3.3, 27, 0.2, 16],
    'common_GO Process': [3, 50, 0.15, 16],
    'common_KEGG': [2, 1.5, 0.1, 12],
    'exclusive_GO Component': [1.2, 1.1, 0.1, 4.5],
    'exclusive_GO Function': [3.2, 14, 0.2, 7],
    'exclusive_GO Process': [3.3, 17, 0.15, 7],
    'exclusive_KEGG': [0.8, -0.6, 0.07, 1.5],
        }


figsize_dict = {
    'all_GO Component': ( 10, 16 ),
    'all_GO Function':  ( 24, 38 ),
    'all_GO Process':   ( 33, 55 ),
    'all_KEGG':         ( 12, 16 ),
    'common_GO Component': ( 10, 16 ),
    'common_GO Function':  ( 24, 38 ),
    'common_GO Process':   ( 33, 55 ),
    'common_KEGG':         ( 12, 16 ),
    'exclusive_GO Component': ( 10, 14 ),
    'exclusive_GO Function':  ( 24, 38 ),
    'exclusive_GO Process':   ( 33, 55 ),
    'exclusive_KEGG':         ( 12, 16 ),
        }


####################
# gset = gsets[3]




for gset in gsets[:]:
# if True: #just a placeholder for the identation
    df = dict_gsets[gset]
    
    #x, y, width, height

    
    variables = {
        'xaxis':'signal',
        'marker_size':'observed gene count',
        'color':r'$-log10(FDR)$'
        }
    
    
    
    width, height = figsize_dict[region+'_'+gset]
    fs = 20
    print(region+'_'+gset)
    
    
    fig = plt.figure(figsize = (width, height) )   
    ax = fig.add_subplot(1,1,1)
    
    
    
    df['is_significant'] = df['false discovery rate'] < 10**-1
    df = df.query('is_significant')
    df.sort_values(by=[variables['xaxis']],
                   inplace = True, ascending = True)
    
    #Ylabels
    ylabels = df['term description']
    ypos = np.arange(len(ylabels))
    
    
    #Mapping color
    # ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    # ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
    # 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
    # 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    # ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
    # 'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
    # 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper']
    # ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
    # 'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
    # 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper']
    
    
    cmap = plt.cm.cool
    
    array_to_colors = df[variables['color']]
    colors, norm = get_colors(array_to_colors, cmap)
    
    color_top = 'black'
    
    
    #Lines
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
            height = 0.05, #grossura das barras
            color = color_top,
            )
    
    
    
    ax.set_yticks(ypos, ylabels, size = fs)
    ax.tick_params(axis = 'x', labelsize=fs)
    ax.set_title(gset, fontsize = fs+10, fontweight = 'bold')
    ax.set_xlabel(variables['xaxis'].capitalize(), fontsize = fs+5)
    ax.grid(axis='x')
    
    
    
    
    n_legend = 4
    
    array_size_legend = np.linspace(10,
                                    array_to_size.max(),
                                    n_legend)
    
    array_size_legend = np.array([(i//10)*10 for i in array_size_legend])
    size_legend = size_marker(array_size_legend)
    
    for i in range(n_legend):
        ax.scatter(-i-2,-i-2,
                   s = size_legend[i],
                   label = int(array_size_legend[i]),
                   color = 'k')
    
    ax.set_ylim(-1,len(ylabels)+1)
    ax.set_xlim(0)
    ax.legend(loc = 'upper left',  
              bbox_to_anchor=(1.04, 1),
              #https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
              title = variables['marker_size'].capitalize(),
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
    
    print(gset)
    print('-Width:',width)
    print('-Height:',height)
    print('Axes cbar:', coord_axes_cbar)
    
    plt.savefig(f'../enrichment/plots/svgs/{region}_{gset}.svg',
                bbox_inches = 'tight')
    plt.savefig(f'../enrichment/plots/{region}_{gset}.png',
                bbox_inches = "tight")

    plt.clf()


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
