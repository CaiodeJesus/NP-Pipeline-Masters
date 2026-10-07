#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 13:20:59 2025

@author: caio
"""

import requests 
import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import pandas as pd
import os

def Request(call_id,
            output_format,
            method,
            params,
            species_id):
    

    string_api_url = "https://version-12-0.string-db.org/api"
    request_url = "/".join([string_api_url, output_format, method])

    params["species"] = species_id
    params["caller_identity"] = call_id

    response = requests.post(request_url, data=params)
    
    return response




def Get_ids(search_proteins,
                   call_id,
                   species_id,
                   do_print = False):

    output_format = "tsv-no-header"
    method = "get_string_ids"

    params = {
    
        "identifiers" : "\r".join(search_proteins), # your protein list
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
    }
    
    
    #Accessing API
    response = Request(call_id, output_format,
                       method, params, species_id)
        
    #Response processing

    string_dict = {}
    for line in response.text.strip().split("\n"):
        l = line.split("\t")
        # print(l)
        input_identifier, string_identifier = l[0], l[2]
        #header :['queryItem', 0
                # 'queryIndex', 
                # 'stringId', 2
                # 'ncbiTaxonId',
                # 'taxonName',
                # 'preferredName',
                # 'annotation']
    
        string_dict[input_identifier] = string_identifier
        if do_print:
            print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
         
    return string_dict

def Network_Image(search_ids, call_id, species_id,
                  folder_img = './'):

    #Params for task
    output_format = "svg"
    method = "network"
    
    params = {

        "identifiers" : "\r".join(search_ids), # your protein list
        "add_white_nodes": 0, # add 15 white nodes to my protein 
        "network_type": "functional", # show confidence links
        "flat_node_design" : 1,
        "required_score":0
    }

    #Accessing API
    response = Request(call_id, output_format,
                       method, params, species_id)

    #Response processing
    file_name = folder_img + "network.svg" 

    with open(file_name, 'wb') as fh:
        fh.write(response.content)

    # image = Image.open(file_name)

    # return image




def Homology(search_ids, call_id, species_id):

    
    #Params for task
    output_format = "tsv"
    method = "homology"
    
    params = {

        "identifiers" : "%0d".join(search_ids), # your protein list
    }

    #Accessing API
    response = Request(call_id, output_format,
                       method, params, species_id)

    matrix = []
    for line in response.text.strip().split("\n"):
    
        l = line.strip().split("\t")
        matrix.append(l)
        # print(l)
    
    df = pd.DataFrame(matrix[1:], columns = matrix[0])
    for col in matrix[0]:
        if 'score' in col:
            df[col] = df[col].astype(float)
            
            
    df['diff_id'] = df['stringId_A'] != df['stringId_B']
    df = df.query('diff_id')
    df.pop('diff_id')
    
    return df




def Network_Table(search_ids, call_id, species_id,
                  add_nodes = 0):

    
    #Params for task
    output_format = "tsv"
    method = "network"
    
    params = {

        "identifiers" : "%0d".join(search_ids), # your protein list
        "add_nodes":add_nodes,
    }

    #Accessing API
    response = Request(call_id, output_format,
                       method, params, species_id)


    
    matrix = []
    for line in response.text.strip().split("\n"):
    
        l = line.strip().split("\t")
        matrix.append(l)
        # print(l)
    
    df = pd.DataFrame(matrix[1:], columns = matrix[0])
    for col in matrix[0]:
        if 'score' in col:
            df[col] = df[col].astype(float)
    
    
    
    df['joined_col'] = df['stringId_A'] + '_' + df['stringId_B']
    df.set_index('joined_col', inplace = True)
    
        
    all_string_ids = list(df['stringId_A']) + list(df['stringId_B'])
    all_string_ids = np.unique(all_string_ids)
    
 
    
    #Getting homology scores in a separate query
    homology = Homology(search_ids = all_string_ids,
                        call_id = call_id,
                        species_id = species_id)
    homology['joined_col'] = homology['stringId_A'] + '_' + homology['stringId_B']
    homology.set_index('joined_col', inplace = True)

    df['hscore'] = homology['bitscore']
    df = df.fillna(0)
    df = df.reset_index()
    df.pop('joined_col')
      
    
    return df









if __name__ == "__main__":

    

    #Global Parameters
    my_name = 'caio_lqmc'
    name_set = 'test'
    human_id = 9606 #Homo Sapiens
    search_proteins = ["p53", "BRCA1", "cdk2", "Q99835", 'egfr']

    try:
        os.mkdir('../data/test/')
    except:
        print('Output will be saved on ../data/test/')
        
    #Getting STRING ids
    string_ids = Get_ids(search_proteins,
                         call_id = my_name,
                         species_id = human_id,
                         do_print = False)

    #Getting STRING table
    df = Network_Table(search_ids = list(string_ids.values()),
                       call_id = my_name,
                       species_id = human_id,
                       add_nodes = 10)
    
    df.to_csv(f'../data/test/{name_set}.csv')

    all_string_ids = list(df['stringId_A']) + list(df['stringId_B'])
    all_string_ids = np.unique(all_string_ids)

    #Getting the STRING network image
    Network_Image(search_ids = all_string_ids,
                   call_id = my_name,
                   species_id = human_id,
                   folder_img = '../data/test/')
    
    
    #Getting the link for STRING page   
    url = 'https://string-db.org/api/tsv/get_link?'
    url += f"identifiers={'%0d'.join(all_string_ids)}&"
    url += f'species={human_id}&add_white_nodes=0'

    link_to_string = requests.post(url).text.split('\n')[1]
    print('Link:\n',link_to_string)
    
      
