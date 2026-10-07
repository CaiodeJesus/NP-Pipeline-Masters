#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 17:04:45 2025

@author: caio
"""

       
import string_api
import os
import requests 
import numpy as np
import pandas as pd
import time

input_folder = '../../../ortholog_id_pipeline/funnel/results/'

regions = ['common', 'exclusive', 'all']
dict_regions = {}


for region in regions:
    
    infile_src = f'Query_{region}_UniProt.csv'
    print(f'{region.capitalize()} Targets')
    print(f'...Searching {infile_src}')
    
    dict_regions[region] = pd.read_csv(input_folder + infile_src)

    
all_links = ''






#Global Parameters
my_name = 'caio_lqmc'
name_set = 'test'
linf_org_id = 5671

    
print('\t...May take about 3 min each')
print()

start = time.asctime()
print('...Starting at',start)
print()
for region in regions:
    
    target_list = list(dict_regions[region]['accession'])
    
    print('...Doing',region)
    
    
    #Getting STRING table
    df = string_api.Network_Table(search_ids = target_list,
                                   call_id = my_name,
                                   species_id = linf_org_id)
    
    df.to_csv(f"../data/{region}/interactions.csv")
    
    
    
    df = pd.read_csv(f"../data/{region}/interactions.csv")
    
    gene_as_value_A = {v: k for v, k in zip(df['stringId_A'], df['preferredName_A']) }
    gene_as_value_B = {v: k for v, k in zip(df['stringId_B'], df['preferredName_B']) }
    
    
    aliases = pd.DataFrame()
    # aliases
    
    
    
    aliases['NameA'] = gene_as_value_A
    aliases['NameB'] = gene_as_value_B
    aliases['Name'] = aliases['NameA'].astype(str) + '_' + aliases['NameB'].astype(str)
    aliases.pop('NameA')
    aliases.pop('NameB')
    
    names = [ ''.join( np.unique(n.replace('nan','').split('_') ) ) for n in aliases['Name'] ]
    names = np.array(names)
    names[names == ''] = 'no_ppi'
    aliases['Name'] = names    
    aliases.to_csv(f"../data/{region}/aliases.csv")
    
    
    
    
    all_string_ids = list(aliases.index)
    
    
    #Getting the STRING network image
    string_api.Network_Image(search_ids = all_string_ids,
                               call_id = my_name,
                               species_id = linf_org_id,
                               folder_img = f"../data/{region}/")
    
    
    #Getting the link for STRING page   
    url = 'https://string-db.org/api/tsv/get_link?'
    url += f"identifiers={'%0d'.join(all_string_ids)}&"
    url += f'species={linf_org_id}&add_white_nodes=0'
    
    link_to_string = requests.post(url).text.split('\n')[1]
    print('Link:\n',link_to_string)
    
    all_links += region+'\t'+link_to_string+'\n'
    
    
    
    
    
    print()
    
end = time.asctime()
print('...Finishing at',end)   
with open('../data/links.txt','w') as f:
    f.write(all_links)
