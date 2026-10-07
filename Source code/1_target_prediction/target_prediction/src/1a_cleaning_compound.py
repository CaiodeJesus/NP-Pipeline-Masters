#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:03:50 2025

@author: caio
"""




def query(protacxn):
    queries_folder = '../data/query/'
    query_out_name = 'pubchem_%s_protacxn'%(protacxn)
    
    if query_out_name in os.listdir(queries_folder):
        os.remove(queries_folder+query_out_name)
        
    url = "https://rest.uniprot.org/uniprotkb/search?query="
    url += protacxn+"&fields=id,protein_name"
    
    wget.download(url, out = queries_folder+query_out_name)
    
    with open(queries_folder+query_out_name, 'r') as file:
        query = json.load(file)
    
    
    return query










import numpy as np
import pandas as pd
import wget
import os
import json






raw_folder = '../data/raw/compound/'
servers = ['swiss',
           'pharmmapper',
           'stitch',
           'targetnet',
           'superpred',
           'sea',
           'chemmapper']
print('Demora um tempin pra rodar tudo, uns 4 min')


# flag = servers[2]
for flag in servers[:]:
    
    ##SWISS
    
    if flag == 'swiss':
        swiss_raw = pd.read_csv(raw_folder+'SwissTargetPrediction.csv')
        
        
        #separando os nomes de gene que vem juntos
        genes = []
        uniprots = []
        metrics = []
        
        common_name = swiss_raw['Common name']
        uniprots_raw = swiss_raw['Uniprot ID']
        metrics_raw = swiss_raw['Probability*']
    
        
        for i in range(len(swiss_raw)):
            genes += common_name[i].split(' ')
    
            uniprots += uniprots_raw[i].split(' ')
    
            metrics += [metrics_raw[i]]*len(uniprots_raw[i].split(' '))
        
        swiss_clean = pd.DataFrame()
        swiss_clean['Gene'] = genes
        swiss_clean['UniProt'] = uniprots
        swiss_clean['Metric'] = metrics
        
        
        #cortando por um valor de probabilidade
    
        print('swiss', len(swiss_clean))
    
        print('\tcom um range de',
              swiss_clean['Metric'].min(), 'a', swiss_clean['Metric'].max())
        
    
    
        swiss_clean.to_csv('../data/clean/compound/swiss.csv')
    
    
    
    ##PHARMMAPPER
    
    if flag == 'pharmmapper':
        pharmapper_raw = pd.read_csv(raw_folder+'pharmmapper.csv', skiprows=1)
        
        
        pharmapper_clean  = pd.DataFrame()
        pharmapper_clean['UniProt'] = pharmapper_raw['Uniplot']
        
    
        
        
        print('pharmapper', len(pharmapper_clean))
        
        metric_column_names = ['Fit', 'Norm Fit', 'zscore']
        for column in metric_column_names:
            
            metric = pharmapper_raw[column]
            pharmapper_clean['Metric_'+column] = metric
        
            print('\tcom um range de',
                  metric.min(), 'a', metric.max(),'para', column)
        
        
        
        pharmapper_clean.to_csv('../data/clean/compound/pharmapper.csv')
    
    
    
    ##STITCH
    
    if flag == 'stitch':
        print('stitch')
        
        stitch_raw = pd.read_csv(raw_folder+'stitch_interactions (1).tsv', sep='\t')
    
        
        stitch_raw['node2_is_protein'] = [node2 in list(stitch_raw['#node1']) for node2 in stitch_raw['node2']]
        stitch_raw['node2_is_protein'] = stitch_raw['node2_is_protein'].astype(float)
    
        scores = []
        targets = np.unique(stitch_raw['#node1'])
        
        uniprot_ids = np.empty(len(targets)).astype(str)
        found_names = np.empty(len(targets)).astype(str)
        
        
        for i in range(len(targets)):
        #aparentemente esses são só os que aparece como resultado,
        #não os compostos que se coloca de input
    
            node1 = targets[i]                          
    
            # print(node1)
            stitch_raw['has_node'] = stitch_raw['#node1'] == node1
            stitch_raw['has_node'] = stitch_raw['has_node'].astype(float)
    
            stitch_raw['is_viable'] = stitch_raw['has_node']*(1 - stitch_raw['node2_is_protein'])
            stitch_raw['is_viable'] = stitch_raw['is_viable'].astype(bool)
            stitch_sel = stitch_raw.query('is_viable')
            
            score = stitch_sel['combined_score'].max()
            # print(score)
            scores.append(score)
            # print()
            
            #pesquisando na API do UniProt
            data = query(node1+'+human')
        
            uniprot_id = data['results'][0]['primaryAccession']
            name_found = data['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
            # print('\t', name_found, uniprot_id)
        
            uniprot_ids[i] = uniprot_id
            found_names[i] = name_found
            
            
        stitch_clean = pd.DataFrame()
        stitch_clean['Names'] = targets
        stitch_clean['Metric'] = scores
    
        stitch_clean['UniProt'] = uniprot_ids
        stitch_clean['Found_Names'] = found_names
        
        
        print('stitch', len(stitch_clean))
    
        print('\tcom um range de',
              stitch_clean['Metric'].min(), 'a', stitch_clean['Metric'].max())
        
            
        stitch_clean.to_csv('../data/clean/compound/stitch.csv')
    
    
    ##TARGETNET
    
    if flag == 'targetnet':
        target_net_raw = pd.read_csv(raw_folder+'target_net.tsv',sep='\t')
        metric = target_net_raw['Prob']
        
        
        
        target_net_clean = pd.DataFrame()
        target_net_clean['UniProt'] = target_net_raw['Uniprot_ID']
        target_net_clean['Metric'] = metric
        
        print('target net', len(target_net_clean))
        
        print('\tcom um range de',
              metric.min(), 'a', metric.max())
        
        target_net_clean.to_csv('../data/clean/compound/target_net.csv')
    
    
    ##SUPERPRED
    
    if flag == 'superpred':
        superpred_raw = pd.read_csv(raw_folder+'superpred.csv')
        
    
        superpred_clean = pd.DataFrame()
        superpred_clean['UniProt'] = superpred_raw['UniProt ID']
    
        print('superpred', len(superpred_clean))
        
    
        metric_column_names = ['Probability', 'Model accuracy']
        for column in metric_column_names:
            
            
            metric = np.zeros(len(superpred_raw))
        
            for i in range(len(superpred_raw)):
                metric[i] = float(superpred_raw[column][i].replace('%',''))
            
            
            superpred_clean['Metric_'+column] = metric
        
            print('\tcom um range de',
                  metric.min(), 'a', metric.max(),'para', column)
            
    
        superpred_clean.to_csv('../data/clean/compound/superpred.csv')
    
    
    ##SEA
    
    if flag == 'sea':
        print(flag)
        sea_raw = pd.read_csv(raw_folder+'sea-results.csv')
        names_raw = sea_raw['Target ID']
       
        names = sea_raw['Target ID']
        uniprots = []
        
        for i in range(len(sea_raw[:])):
            # print(i)
            name = names[i]
            data = query(name.replace(' ','+'))
        
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
    
    
            # print(name)
            # print('\t', name_found, uniprot_id)
                
            uniprots.append(uniprot_id)
        
        
        
        
        sea_clean = pd.DataFrame()
        
        sea_clean['ID'] = names
        sea_clean['UniProt'] = uniprots
        
        
        print('sea', len(sea_clean))
        metric_column_names = ['P-Value', 'Max Tc', 'Cut Sum', 'Z-Score']
    
        for column in metric_column_names[:]:
            
            if sea_raw[column].dtypes == 'object':
                metric = np.empty(len(sea_raw))
                
                i=0
                for m in sea_raw[column]:
                    metric[i] = float(m.replace(',','.'))
                    i+=1
            else:
                metric = sea_raw[column]
            
            sea_clean['Metric_'+column] = metric
        
            print('\tcom um range de',
                  metric.min(), 'a', metric.max(),'para', column)
        
        
        column = 'neg_logPv'
        metric = -np.log10(sea_clean['Metric_P-Value'])
        sea_clean['Metric_'+column] = metric
    
        print('\tcom um range de',
              metric.min(), 'a', metric.max(),'para', column)
        
        
        sea_clean.to_csv('../data/clean/compound/sea.csv')
        
        
        
    ##CHEMMAPPER
    
    if flag == 'chemmapper':
        print(flag)
        
        chemmapper_raw = pd.read_csv(raw_folder+'chemmapper_new.csv')
        
        chemmapper_raw.dropna(subset=['uniprot'], inplace = True)
           
        
        chemmapper_clean = pd.DataFrame()
        chemmapper_clean['UniProt'] = chemmapper_raw['uniprot']
        
        
        print('chemmapper', len(chemmapper_clean))
        
        metric_column_names = ['score', 'simiscore']
    
        for column in metric_column_names[:]:
            
    
            metric = chemmapper_raw[column]
            
            chemmapper_clean['Metric_'+column] = metric
        
            print('\tcom um range de',
                  metric.min(), 'a', metric.max(),'para', column)
            
            
    
        
        chemmapper_clean.to_csv('../data/clean/compound/chemmapper.csv')
        
    print()