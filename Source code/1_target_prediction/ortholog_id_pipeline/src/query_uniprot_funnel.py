#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:21:14 2026

@author: caio
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import wget


with open('../support/list_of_ids_tritrypdb','r') as f:
    total_list_of_ids = [leish.strip('\n') for leish in f.readlines()]
    
# id_list = []
do_query = False
queries_folder = '../query/'



with open('../results/final_leish_common_targets_ttdb','r') as f:
    common_targets = [leish.strip('\n') for leish in f.readlines()]
    
with open('../results/final_leish_exclusive_targets_ttdb','r') as f:
    exclusive_targets = [leish.strip('\n') for leish in f.readlines()]

dict_regions = {
    'common': common_targets,
    'exclusive': exclusive_targets,
    'all':common_targets + exclusive_targets
                }

if do_query:

    for ttdb_id in total_list_of_ids[:]:
        # print(ttdb_id)
        if ttdb_id+'.txt' not in os.listdir(queries_folder):
    
            url = 'https://rest.uniprot.org/uniprotkb/search?query='
            url +='xref:veupathdb-'+ttdb_id
            url += "&fields=accession,gene_names,organism_name,xref_veupathdb"
            # print(url)
            wget.download(url, out = queries_folder+ttdb_id+'.txt')
    
    print('Every ID in Query!')

else:

    
    for key in dict_regions.keys():
        
        print(key.capitalize(),'targets.')
        cols = {
                'ttdb_id':[],
                'accession':[],
                'organism_name':[],
                'xref_veupathdb':[]
                }
        
        something_wrong = {}
        
        for ttdb_id in dict_regions[key][:]:
            
                if ttdb_id+'.txt' in os.listdir(queries_folder):
                    with open(queries_folder+ttdb_id+'.txt', 'r') as file:
                        query = json.load(file)
                else:
                    print('...Missing',ttdb_id)
                    continue
            
                
                try:
                    accession = query['results'][0]['primaryAccession']
                    
                    org_name = query['results'][0]['organism']['scientificName']
                    
                    x_ref_list = query['results'][0]['uniProtKBCrossReferences']
                    
                    xref_ttdb = ''
                    for x_ref in x_ref_list:
                        if x_ref['id'] == 'TriTrypDB:' + ttdb_id:
                            xref_ttdb = x_ref['id']
    
                except:
                    
                    accession = ''
                    org_name = ''
                    xref_ttdb = ''
                                                                                           
                    print('...',ttdb_id,'Deu ruim')
        
        
                cols['ttdb_id'].append(ttdb_id)
                cols['accession'].append(accession)
                cols['organism_name'].append(org_name)
                cols['xref_veupathdb'].append(xref_ttdb)
                
        
        df = pd.DataFrame(cols)
        #conferindo
        
        print('...Organismos:')
        
        diff_org = 0
        for org in np.unique(df['organism_name']):
            print('.....',org)
            diff_org +=1
        print('...Deveria ser só leishmania')
        
        if diff_org == 1:
            print('.....','E é!')
            df.pop('organism_name')
        else:
            print('...','..','tem',diff_org-1,'a mais')
            
        print()
        print('...Qual pegou uma ID diferente?:')
        
        count = 0
        for i in range(len(df)):
            line = df.iloc[i]
            if line['xref_veupathdb'] != 'TriTrypDB:'+line['ttdb_id']:
                print('',line['ttdb_id'],line['xref_veupathdb'])
                count+=1
                
        if count == 0:
            print('...','..','Nenhuma')
            df.pop('xref_veupathdb')
        else:
            print('...','..',count, 'no total')
        
        
        df.to_csv(f'../results/Query_{key}_UniProt.csv')
        print()
    
    
    
    
    