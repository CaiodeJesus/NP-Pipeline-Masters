#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:56:13 2026

@author: caio
"""

import pandas as pd
import os
import wget
import json



do_query = False


dict_regions = {}

folder = '../out/infiles/'
queries_folder = '../query/'

for file in os.listdir(folder):
    
    region = file.split('.')[0]
    
    df = pd.read_csv(folder+file,
                     skiprows=2, header = None, sep = '\t')
    
    dict_regions[region] = [uni.split(':')[1] for uni in df[2]]
    
    
if do_query:

    for region in dict_regions:
        for uniprot_id in dict_regions[region][:]:
            # print(uniprot_id)
            if uniprot_id+'.txt' not in os.listdir(queries_folder):
        
                url = 'https://rest.uniprot.org/uniprotkb/search?query='
                url += uniprot_id
                # print(url)
                wget.download(url, out = queries_folder+uniprot_id+'.json')
        
    print('Every ID in Query!')
    
    
else:

    
    for key in dict_regions.keys():

        print(key.capitalize(),'targets.')
        cols = {
            'Query_Acess': [],
            'Organism': [],
            'Found_Acess': [],
            'Existence': [],
            'Sequence': [],
            'Name': [],
            'E.C.': [],
            'PDB': [],
            'GO_ids': [],
            'GO_terms': [],
            'AlphaFoldDB': []
                }
        
        something_wrong = {}
        
        for uniprot_id in dict_regions[key][:]:
            
                do_print_erro = True
                
                with open(queries_folder+uniprot_id+'.json', 'r') as file:
                    query = json.load(file)

                try:
                    org = query['results'][0]['organism']['scientificName']
                except:
                    org = '-'
                    # print(uniprot_id,'Deu ruim')
                    
                    
                try:
                    unipros_acess = query['results'][0]['primaryAccession']
                except:
                    unipros_acess = '-'
                    # if do_print_erro:
                    #     print(uniprot_id,'Deu ruim')
                    #     do_print_erro = False

                try:
                    existence = query['results'][0]['proteinExistence']
                except:
                    existence = '-'
                    # if do_print_erro:
                    #     print(uniprot_id,'Deu ruim')
                    #     do_print_erro = False


                try:
                    seq = query['results'][0]['sequence']['value']
                except:
                    seq = '-'
                    # if do_print_erro:
                    #     print(uniprot_id,'Deu ruim')
                    #     do_print_erro = False

                try:
                    name = query['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                except:
                    name = '-'
                    # if do_print_erro:
                    #     print(uniprot_id,'Deu ruim')
                    #     do_print_erro = False
                        
                        
                try:
                    eclist = query['results'][0]['proteinDescription']['recommendedName']['ecNumbers']
                    
                    ecs = [ec['value'] for ec in eclist]
                    ecs = '|'.join(ecs)
                    if ecs == '':
                        ecs = 'no_data'
                    
                except:
                    ecs = '-'
                    # if do_print_erro:
                    #     print(uniprot_id,'Deu ruim')
                    #     do_print_erro = False
                
                try:
                    cross_ref = query['results'][0]['uniProtKBCrossReferences']
                    
                    pdbs = '|'.join([db['id'] for db in cross_ref if db['database'] == 'PDB'])
                    go_ids = '|'.join([db['id'] for db in cross_ref if db['database'] == 'GO'])
                    go_terms = '|'.join([db['properties'][0]['value'] for db in cross_ref if db['database'] == 'GO'])
                    afs = '|'.join([db['id'] for db in cross_ref if db['database'] == 'AlphaFoldDB'])
                    
                    if pdbs == '':
                        pdbs = 'no_data'
                    if go_ids == '':
                        go_ids = 'no_data'
                    if go_terms == '':
                        go_terms = 'no_data'
                    if afs == '':
                        afs = 'no_data'
                        
                except:
                    
                    pdbs = '-'
                    go_ids = '-' 
                    go_terms = '-' 
                    afs = '-'
                
                cols['Query_Acess'].append(uniprot_id)
                cols['Organism'].append(org)
                cols['Found_Acess'].append(unipros_acess)
                cols['Existence'].append(existence)
                cols['Sequence'].append(seq)
                cols['Name'].append(name)
                cols['E.C.'].append(ecs)
                cols['PDB'].append(pdbs)
                cols['GO_ids'].append(go_ids) 
                cols['GO_terms'].append(go_terms) 
                cols['AlphaFoldDB'].append(afs)

                
        
        out = pd.DataFrame(cols)
        out.to_csv(f'../out/data/{key}.csv')