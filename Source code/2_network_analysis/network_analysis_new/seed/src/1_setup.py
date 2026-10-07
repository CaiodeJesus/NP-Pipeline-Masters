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



input_folder = '../input/'
infile_src = 'Leish_targets_pass.txt'

all_links = ''
    
    
infile = infile_src.split('.')[0] #tirando o .txt
print(f'Searching {infile}')
    
#Making directories for each overlap
if infile not in os.listdir('../data/'):
    os.mkdir('../data/' + infile)
    # print('../data/'+area.split('.')[0])
    

#Reading Target Files
with open(f"{input_folder}/{infile_src}") as f:
    target_list = f.readlines()        

target_list = [p.strip('\n') for p in target_list if p[0] != '#']
    

#Global Parameters
my_name = 'caio_lqmc'
name_set = 'test'

if 'LEISH' in infile.upper():
    org_id = 5671
    print('\t...Leishmania infantum PPI')
elif 'TRYP' in infile.upper():
    org_id = 353153
    print('\t...Trypanosoma cruzi PPI')
else:
    org_id = 9606
    print('\t...Searching Human PPI')
    
print('\t...May take about 3 min')

        


# #Getting STRING table
# df = string_api.Network_Table(search_ids = target_list,
#                                call_id = my_name,
#                                species_id = org_id)

# df.to_csv(f"../data/{infile}/interactions.csv")



df = pd.read_csv(f"../data/{infile}/interactions.csv")

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
aliases.to_csv(f"../data/{infile}/aliases.csv")




all_string_ids = list(aliases.index)


#Getting the STRING network image
string_api.Network_Image(search_ids = all_string_ids,
                           call_id = my_name,
                           species_id = org_id,
                           folder_img = f"../data/{infile}/")


#Getting the link for STRING page   
url = 'https://string-db.org/api/tsv/get_link?'
url += f"identifiers={'%0d'.join(all_string_ids)}&"
url += f'species={org_id}&add_white_nodes=0'

link_to_string = requests.post(url).text.split('\n')[1]
print('Link:\n',link_to_string)

all_links += infile+'\t'+link_to_string+'\n'





print()

with open('../data/links.txt','w') as f:
    f.write(all_links)
