#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:30:27 2025

@author: caio
"""



def query(search, out):
    queries_folder = '../data/query/'
    query_out_name = out
    
    if query_out_name in os.listdir(queries_folder):
        os.remove(queries_folder+query_out_name)
        
    url = "https://rest.uniprot.org/uniprotkb/search?query="
    url += search + "&fields=id,protein_name"
    print(url)
    
    try:
        wget.download(url, out = queries_folder+query_out_name)
        with open(queries_folder+query_out_name, 'r') as file:
            query = json.load(file)

    except:
        query = {'results': []}
        print('Deu ruim a pesquisa')
        

    
    
    return query





import numpy as np
import pandas as pd
import wget
import os
import json


#clean é quando eu consigo
#abrir em csv o arquivo, pegar o UniProt e o Metric de um modo padronizado



raw_folder = 'raw/disease/'
flag = 'disgenet'

diseases_list = ['leish', 'chagas']
raw_folder = '../data/raw/disease/'

servers = ['genecards',
           'omim',
           'disgenet',
           'ncbi',
           'malacards']



for flag in servers:
    print(flag)
    ##GENECARDS
    
    if flag == 'genecards':
    
        for disease in diseases_list[:]:
    
            gene_cards_disease_raw = pd.read_csv(raw_folder + 'GeneCards-SearchResults_%s.csv'%disease)
            
            
            gene_cards_disease_clean = pd.DataFrame()
            gene_cards_disease_clean['UniProt'] = gene_cards_disease_raw['Uniprot ID']
            gene_cards_disease_clean['Metric'] = gene_cards_disease_raw['Relevance score']
            gene_cards_disease_clean.dropna(inplace = True)
            gene_cards_disease_clean.to_csv('../data/clean/disease/%s/genecards.csv'%disease)
    
    
    ##OMIM
    
    if flag == 'omim':
        
        for disease in diseases_list[:]:
            print(disease)
    
            omim_disease_raw = pd.read_csv(raw_folder + 'OMIM-Entry-Retrieval_%s.tsv'%disease,
                                           sep = '\t', header = 3)
            
            
            omim_disease_raw.dropna(subset = ['Cytogenetic Location'], inplace = True)
            #filtranso por aqui pq as vezes ele lista doença e não alvo,
            #e conferindo esses que dão problema não tem nan nessa coluna
            omim_disease_raw.reset_index(inplace = True)
            
        
            
            
            mim_accession_list = np.empty(len(omim_disease_raw)).astype(str)
            uniprot_list = np.empty(len(omim_disease_raw)).astype(str)
            found_names_list = np.empty(len(omim_disease_raw)).astype(str)
        
        
            for i in range(len(omim_disease_raw)):
        
                print(omim_disease_raw['Title'][i])
                
                mim_number_i = omim_disease_raw['MIM Number'][i]
                
                mim_accession_list[i] = mim_number_i
                mim_number_i = mim_number_i.strip('%').strip('#').strip('*').strip('+')
                
                
                
                try:
                    title_i, gene_name_i = omim_disease_raw['Title'][i].split(';')
                    title_i = title_i.replace(',','').replace(' ','+')
                    gene_name_i = gene_name_i.strip(' ')
        
                except:
                    title_i = omim_disease_raw['Title'][i]
        
                
                print('\t', title_i)
                print('\t', mim_number_i)
                print('\t', gene_name_i)
                
                search = title_i + '+AND+gene:'
                search += gene_name_i + '+AND+xref:mim-' + mim_number_i 
                
                print('\t'+search)
                data = query(search,title_i)
                print()
                
                print(len(data['results']), 'resultados foram encontrados')
                
                try:
                    uniprot_id = data['results'][0]['primaryAccession']
                    try:
                        name_found = data['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                    except:
                        name_found = data['results'][0]['proteinDescription']['submissionNames'][0]['fullName']['value']
                        print('except')
                except:
                    uniprot_id = ''
                    name_found = ''
            
                print('Uniprot',uniprot_id)
                uniprot_list[i] = uniprot_id
                found_names_list[i] = name_found
            
            
            
            omim_disease_clean  = pd.DataFrame()
            omim_disease_clean['UniProt'] = uniprot_list
            omim_disease_clean['Metric'] = np.zeros(len(omim_disease_clean))
            omim_disease_clean['MIM_number'] = omim_disease_raw['MIM Number']
            omim_disease_clean['Names'] = omim_disease_raw['Title']
            omim_disease_clean['Found_Names'] = found_names_list
           
            omim_disease_clean.to_csv('../data/clean/disease/%s/omim.csv'%disease)
            
        
    
            
    ##DISGENET
    if flag == 'disgenet':
        print('')
        
        # disease = diseases_list[0]
        
        for disease in diseases_list[:]:
            print(disease)
            
            disgenet_disease_raw = pd.read_csv(raw_folder + 'disgenet_%s.tsv'%disease,
                                               sep = '\t')
     
    
            disgenet_disease_clean = pd.DataFrame()
            disgenet_disease_clean['UniProt'] = disgenet_disease_raw['UnitProt']
            
            metric_column_names = ['EvidenceIndexGDA',
                                   'GeneDPI',
                                   'GeneDSI',
                                   'GenePLI',
                                   'NumDiseasesAssociatedToGene',
                                   'ScoreGDA']
            
            for column in metric_column_names:
                
                
                metric = disgenet_disease_raw[column]                
                disgenet_disease_clean['Metric_'+column] = metric
            
                print('\tcom um range de',
                      metric.min(), 'a', metric.max(),'para', column)
                print('\túltimo elemtento = ',list(metric)[-1]) #pra mostrar que não é igual
    
    
            print()
            
            
            keep_line = np.zeros(len(disgenet_disease_clean))
            new_lines = []
            for i in range(len(disgenet_disease_clean['UniProt'])):
                uniprot = disgenet_disease_clean['UniProt'][i]
                if ',' in uniprot:
                    # print(uniprot)
                    line_metrics = list(disgenet_disease_clean.iloc[i])[1:]
                    # print(metrics)
                    keep_line[i] = False
                    for uni in uniprot.split(','):
                        new_line = [uni] + line_metrics
                        # print(new_line)
                        new_lines.append(new_line)
                else:
                    keep_line[i] = True
                
            keep_line = np.append(keep_line,
                                  np.ones(len(new_lines)))
            
            new_lines_df = pd.DataFrame(new_lines, columns = disgenet_disease_clean.columns)
            disgenet_disease_clean = pd.concat([disgenet_disease_clean, new_lines_df])
            
            disgenet_disease_clean['keep_line'] = keep_line.astype(bool)
            disgenet_disease_clean = disgenet_disease_clean.query('keep_line')
            disgenet_disease_clean.reset_index( inplace = True )
            disgenet_disease_clean.pop('keep_line')
            
            disgenet_disease_clean.to_csv('../data/clean/disease/%s/disgenet.csv'%disease)
                
    
    
    ##NCBI
    
    if flag == 'ncbi':
        
        for disease in ['leish', 'chagas']:
            
            ncbi_disease_raw = pd.read_csv(raw_folder+'NCBI_gene_result_'+disease+'.txt',
                                           sep='\t')
            
            
            #selecting human
            ncbi_disease_raw['is_human'] = ncbi_disease_raw['Org_name'] =='Homo sapiens'
            ncbi_disease_raw = ncbi_disease_raw.query('is_human')
            ncbi_disease_raw.pop('is_human')
            ncbi_disease_raw.reset_index(inplace = True)
            print()
            
            uniprot_ids = []
            found_names = []
            geneid_list = []
            for i in range(len(ncbi_disease_raw)):
                geneid = str(ncbi_disease_raw['GeneID'][i])
                search = 'xref:geneid-' + geneid
                data = query(search, out = 'ncbi_'+search)
                
                try:
                    uniprot_id = data['results'][0]['primaryAccession']
                    try:
                        name_found = data['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                    except:
                        name_found = data['results'][0]['proteinDescription']['submissionNames'][0]['fullName']['value']
                        print('except')
                except:
                    uniprot_id = ''
                    name_found = ''
                
                uniprot_ids.append(uniprot_id)
                found_names.append(name_found)
                geneid_list.append(geneid)
                print(geneid,uniprot_id,name_found)
                
                
            ncbi_disease_clean = pd.DataFrame()
            ncbi_disease_clean['UniProt'] = uniprot_ids
            ncbi_disease_clean['Metric'] = np.zeros(len(ncbi_disease_clean))
            ncbi_disease_clean['GeneID'] = geneid_list
            
            ncbi_disease_clean.to_csv('../data/clean/disease/%s/ncbi.csv'%disease)
        
    ##MALACARDS
    
    if flag == 'malacards':
        
        for disease in ['leish', 'chagas']:
            malacards_disease_raw = pd.read_csv(raw_folder+'malacards_copypaste_' + disease + '.tsv',
                                                sep = '\t', index_col = False)
            
            columns = list(malacards_disease_raw)
        
            malacards_disease_raw.dropna( subset = [columns[1]], inplace = True)
            malacards_disease_raw.reset_index( inplace = True )
            
            gene_list = []
            uniprot_ids = []
            found_names = []
            
            for i in range(len(malacards_disease_raw)):
                gene = str(malacards_disease_raw[columns[1]][i])
                
                search = 'gene:' + gene
                data = query(search, out = 'ncbi_'+search)
                
                try:
                    uniprot_id = data['results'][0]['primaryAccession']
                    try:
                        name_found = data['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                    except:
                        name_found = data['results'][0]['proteinDescription']['submissionNames'][0]['fullName']['value']
                        print('except')
                except:
                    uniprot_id = ''
                    name_found = ''
                
                uniprot_ids.append(uniprot_id)
                found_names.append(name_found)
                gene_list.append(gene)
                print(gene,uniprot_id,name_found)
    
    
    
            malacards_disease_clean = pd.DataFrame()
            malacards_disease_clean['UniProt'] = uniprot_ids
            malacards_disease_clean['Metric'] = np.zeros(len(malacards_disease_clean))
            malacards_disease_clean['GeneID'] = gene_list
            
            malacards_disease_clean.to_csv('../data/clean/disease/%s/malacards.csv'%disease)        
    
        
        