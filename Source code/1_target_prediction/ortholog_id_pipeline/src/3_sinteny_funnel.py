#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:52:50 2026

@author: caio
"""

import pandas as pd
import numpy as np

synteny_raw = pd.read_csv('../data/GeneByLocusTag_OrthologsLite.csv')

trypsp_ref = 'Trypanosoma cruzi CL Brener Esmeraldo-like'

synteny_raw['is_ref_tcruzi'] = synteny_raw['Organism'] == trypsp_ref
synteny_raw['is_syntenic'] = synteny_raw['is syntenic'] == 'yes'
synteny_raw['do_pass'] = synteny_raw['is_ref_tcruzi']*synteny_raw['is_syntenic']
synteny_raw['do_pass'] = synteny_raw['do_pass']

common = synteny_raw.query('do_pass').copy()

cols_to_pop = common.columns

common['Linf'] = common['Gene ID']
common['Tcruzi'] = common['Ortholog']

for col in cols_to_pop:
    common.pop(col)

out_file = '../results/common/Step3.csv'
common.to_csv(out_file)


print('3. Synteny')
print(len(common),'Leish-Tryp pairs')
for col in common.columns:
    print('...',len(np.unique(common[col])),'unique',col)
    
    


        
        
    