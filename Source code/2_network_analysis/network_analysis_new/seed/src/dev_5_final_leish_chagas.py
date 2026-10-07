#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 17:27:03 2025

@author: caio
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from string_api import Network_Table, Network_Image
import os
import networkx as nx
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def Clean_Rank(overlap):
    raw_rank = pd.read_csv(f'../out/overlaps/{overlap}.txt',
                      sep = '\t',
                      skiprows = 2,
                      header = None)

    uniprots_clean = [up.split(':')[1] for up in raw_rank[2]]

    rank_clean = pd.DataFrame()
    rank_clean['Rank'] = raw_rank[0]
    rank_clean['Name'] = raw_rank[1]
    rank_clean['Uniprot'] = uniprots_clean
    
    return rank_clean

# #Gathering the three regions
# leish_comp = Clean_Rank('')
# middle = Clean_Rank('middle_22')
# exclusive = Clean_Rank('exclusive_leish_compound_26')


# middle['Region'] = 'in Chagas, middle'
# exclusive['Region'] = 'not in Chagas, exclusive'

# regions = pd.concat([middle, exclusive])
# regions.set_index('Uniprot', inplace = True)


# text = ''

# region_list = np.empty(len(leish_comp)).astype(str)
# pos_region_list = np.empty(len(leish_comp)).astype(str)

# #Writing information to Output
# for i in range(len(leish_comp)):
#     line = leish_comp.iloc[i]

#     text += '\t'.join(list(line)[:-1])
    
#     uniprot = line['Uniprot']
        
    
#     region = regions.loc[uniprot]['Region']
#     pos_region = regions.loc[uniprot]['Rank']
    
#     region_list[i] = region.split(',')[1].strip()
#     pos_region_list[i] = pos_region
    
#     text += f"\t{region}'s {pos_region}"
    
#     text += f"\tuniprot.org/uniprotkb/{uniprot}/entry\n"


# leish_comp['Region'] = region_list
# leish_comp['Rank_region'] = pos_region_list


# with open('../input') as f:
#     input_lines = f.readlines()


# today = datetime.datetime.now()
# header = '##################################\n'
# header += '#INPUT PARAMETERS\n'
# header += '#'+today.strftime("%c") +'\n'
# header += '#'+'#'.join(input_lines).replace('-','#') #already has \n
# header += '##################################'
# columns = 'Overall_rank\tName\tVenn_position\tUniProt_Link'

# print(header)
# print(columns)
# print(text)
    

# with open(f"../out/out_leish_chagas_{today.strftime('%y_%m_%d')}", 'w') as f:
#     f.write(header+'\n')
#     f.write(columns+'\n')
#     f.write(text)
    
    
# my_name = 'caio_lqmc'
# name_set = 'test'
# human_id = 9606 #Homo Sapiens


# #Searching STRING API for final network
# aliases = pd.read_csv('../data/leish_compound_48/aliases.csv')
# aliases = aliases.set_index(['Uniprot'])
# leish_comp = leish_comp.set_index(['Uniprot'])
# leish_comp['string_id'] = aliases['Unnamed: 0']


# is_table = 'table.csv' in os.listdir('../out/final_net/')
# is_image = 'network.svg' in os.listdir('../out/final_net/')
# are_both_in_folder = is_table and is_image


# if not are_both_in_folder:

#     print('Getting Network Table')
#     df = Network_Table(search_ids = list(leish_comp['string_id']),
#                        call_id = my_name,
#                        species_id = human_id,
#                        add_nodes = 0)
    
#     df.to_csv('../out/final_net/table.csv')
    
#     all_string_ids = list(df['stringId_A']) + list(df['stringId_B'])
#     all_string_ids = np.unique(all_string_ids)
    
#     #Getting the STRING network image
#     print('Getting Network Image')
#     Network_Image(search_ids = all_string_ids,
#                    call_id = my_name,
#                    species_id = human_id,
#                    folder_img = '../out/final_net/')


# #Creating Graph type from final network
# df_channel = pd.read_csv('../out/final_net/table.csv')
# G = nx.Graph()
# for i in range(len(df_channel)):
#     line = df_channel.iloc[i]
#     n1, n2 = line['preferredName_A'], line['preferredName_B']
#     weight = line['score']

#     G.add_edge(n1, n2, weight=weight)





# #Drawing Network

# fig = plt.figure(figsize = (40,15))   
   
# ax1 = fig.add_subplot(1,1, 1)

# pos = nx.spring_layout(G, k = 2.5, seed = 80,
#                        weight = 'weight' )

# edge_width = 2
# node_size = 1000
# fs = 7

# big_node, small_node = 4000, 700



# nodes = np.array((leish_comp['Name']))

# node_ranks = [float(r.strip('º')) for r in leish_comp['Rank']]
# node_ranks = np.array(node_ranks)

# region_ranks = [float(r.strip('º')) for r in leish_comp['Rank_region']]
# region_ranks = np.array(region_ranks)
# region_nodes = np.array(leish_comp['Region'])

# N = node_ranks.max()
# node_sizes = 3*(N-node_ranks)*(big_node-small_node)/(N-1)+ small_node



# #Drawing nodes from the middle region (with Chagas)
# cmap_middle = mpl.colormaps['Greens']
# colors_middle = cmap_middle(np.linspace(0,1,20))
# node_colors = [colors_middle[int(20-r)-1] for r in region_ranks[region_nodes == 'middle']]

# nx.draw_networkx_nodes(G, pos, 
#                        nodelist = nodes[region_nodes == 'middle'],
#                        node_size = node_sizes[region_nodes == 'middle'],
#                        alpha = 1,
#                        node_color = node_colors,
#                        ax = ax1, linewidths=0)

# for j in range(len(node_colors)):
#     nc = node_colors[j]
#     nx.draw_networkx_labels(G.subgraph(nodes[region_nodes == 'middle'][j:j+1]),
#                             font_size=20, font_weight = 'bold',
#                             font_color = 'white',
#                             pos=pos, 
#                             bbox={"ec": "None", "fc": nc , "alpha": 1},
#                             ax  = ax1)
#     #Iterating for eah node so that the label box is the same color as the node

# #Drawing nodes from the exclusive region (not in Chagas)
# cmap_exclusive = mpl.colormaps['RdPu']
# colors_exclusive = cmap_exclusive(np.linspace(0,1,20))
# node_colors = [colors_exclusive[int(20-r)-1] for r in region_ranks[region_nodes == 'exclusive']]
# nx.draw_networkx_nodes(G, pos, 
#                        nodelist = nodes[region_nodes == 'exclusive'],
#                        node_size = node_sizes[region_nodes == 'exclusive'],
#                        alpha = 1,
#                        node_color = node_colors,
#                        ax = ax1, linewidths=0,
#                        )


# for j in range(len(node_colors)): 
#     nx.draw_networkx_labels(G.subgraph(nodes[region_nodes == 'exclusive'][j:j+1]),
#                             font_size=20, font_weight = 'bold',
#                             pos=pos, 
#                             font_color = 'white',
#                             bbox={"ec": "None", "fc": node_colors[j] , "alpha": 1},
#                             ax  = ax1)


# # Drawing edges in 3 groups (low, medium and high confidence according to STRING classification)
# lc, hc = 0.4, 0.7 #low and high confidence limits
# edge_low = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] < lc]
# edge_medium = [(u, v) for (u, v, d) in G.edges(data=True) if (d["weight"] >= lc and d["weight"] < hc)]
# edge_high = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] >= hc]

# # .: ‘-’, ‘–’, ‘-.’, ‘:’ or w
# nx.draw_networkx_edges(G, pos, edgelist = edge_low,
#                        width=edge_width, style = ':',
#                        ax = ax1, alpha = 0.5, label = 'Low')
# nx.draw_networkx_edges(G, pos, edgelist = edge_medium,
#                        width=edge_width, style = '--',
#                        ax = ax1, alpha = 0.5, label = 'Medium')
# nx.draw_networkx_edges(G, pos, edgelist = edge_high,
#                        width=1.5*edge_width, style = '-',
#                        ax = ax1, alpha = 1, label = 'High')


# #Drawing legend
# ax1.legend(loc = 'lower center',
#           ncols = 3,
#           fontsize = 20,
#           title = 'Interaction Confidence Score',
#           title_fontsize = 23,
#           bbox_to_anchor=(0.5, -0.1),
#           fancybox = True)

# ax1.axis('off')
# ax1.set_title('Protein-Protein Interaction Network\nof Leishmania-Hit Targets',
#           weight = 'bold', size = 30)


# #Drawing colormaps for other region ranks
# ranks_image = []

# l = 5

# for i in range(len(region_ranks)):
#     r = region_ranks[i]

#     region = region_nodes[i]
    
#     if region == 'middle':
#         pixel = colors_middle[int(20-r)-1]
#         ranks_image.append(l*[pixel])
        
#     elif region == 'exclusive':
#         pixel = colors_exclusive[int(20-r)-1]
#         ranks_image.append(l*[pixel])
        
#     else:
#         print('Something is wrong')
#         pixel = (0,0,0,0)
#         ranks_image.append(l*[pixel])


# divider = make_axes_locatable(ax1)
# ax2 = divider.append_axes("left", size="15%", pad=0.5)
# imshow = ax2.imshow(ranks_image)
# ax2.set_title('Leishmania-Hit\nRanking\n',
#               fontsize = 23)#, weight = 'italic')


# big_rank_ticks = [f"{int(r)}º {n}" for r,n in zip(node_ranks,nodes)]

# ax2.set_yticks(np.arange(len(nodes)), big_rank_ticks,
#                fontsize = 20)
# ax2.set_xticks([])



# cax1 = divider.append_axes("right", size="2%", pad=1)
# cbar1 = mpl.colorbar.ColorbarBase(cax1, orientation='vertical', 
#                                   cmap='Greens')
# rank_yticks = [f"{r}º" for r in np.arange(20)+1][::-1]
# cbar1.ax.set_yticks(np.linspace(0,1,20),rank_yticks)
# cbar1.ax.tick_params(labelsize=20)
# cax1.set_title('Middle\nRanking\n',
#               fontsize = 23)#, weight = 'italic')


# cax2 = divider.append_axes("right", size="2%", pad=2)
# cbar2 = mpl.colorbar.ColorbarBase(cax2, orientation='vertical', 
#                                cmap='RdPu')
# cbar2.ax.set_yticks(np.linspace(0,1,20),rank_yticks)
# cbar2.ax.tick_params(labelsize=20)
# cax2.set_title('Exclusive LH\nRanking\n',
#               fontsize = 23)#, weight = 'italic')


# plt.savefig('../out/Ranking_plot.svg')
# nx.write_graphml_lxml(G, '../out/Ranking_graph.graphml')

# pd.DataFrame(pos).to_csv('../out/positions.csv')
