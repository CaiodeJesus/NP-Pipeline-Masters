#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:36:29 2026

@author: caio
"""

import numpy as np
import pandas as pd
from Bio import SeqIO

with open('../support/list_of_ids_tritrypdb','r') as f:
    list_of_ids_raw = f.readlines()    
    
df_ids = pd.DataFrame()
# list_of_ids_raw.append('LINF_010005100')
df_ids['ID'] = ['linf|'+ttdb_id.strip('\n') for ttdb_id in list_of_ids_raw]
df_ids.set_index('ID', inplace = True)
df_ids['is_in_tp_leish'] = 1


shared_og = pd.read_csv('../data/Query_leish_BySharedOrtholog.csv')


shared_og.set_index('Accession', inplace = True)
shared_og['is_in_tp_leish'] = df_ids['is_in_tp_leish']
shared_og.fillna(0, inplace = True)

shared_og['is_in_tp_leish'] = shared_og['is_in_tp_leish'].astype(bool)


common_area = shared_og.query('is_in_tp_leish').copy()
common_area.reset_index(inplace = True)


org_dict = {
    'tcba' : 'Trypanosoma cruzi Brazil A4', #TcBrA4_0030010 tcba
    'tcbe' : 'Trypanosoma cruzi Berenice', #ECC02_000231 tcbe
    'tccd' : 'Trypanosoma cruzi Dm28c 2017', #BCY84_00916 tccd
    'tccl' : 'TcCL_ESM04623_t1', #TcCL_ESM04623 tccl
    'tcrg' : 'Trypanosoma cruzi strain G', #TcG_00364 tcrg
    'tcrk' : 'Trypanosoma cruzi TCC', #C3747_41g236 tcrk
    'tcrl' : 'Trypanosoma cruzi CL Brener Esmeraldo-like', #TcCLB.511577.83 tcrl
    'tcrm' : 'Trypanosoma cruzi marinkellei strain B7', #Tc_MARK_1872 tcrm
    'tcrn' : 'Trypanosoma cruzi CL Brener Non-Esmeraldo-like', #TcCLB.504797.40 tcrn
    'tcrs' : 'Trypanosoma cruzi Sylvio X10/1', #TcSYL_0109540.t1 tcrs
    'tcru' : 'Trypanosoma cruzi Dm28c 2018', #C4B63_42g75 tcru
    'tcsx' : 'Trypanosoma cruzi Sylvio X10/1-2012', #TCSYLVIO_003162 tcsx
    'tcyc' : 'Trypanosoma cruzi Y C6', #TcYC6_0095990 tcyc    
}  


total_tcruzi = list(common_area['Target ID'])
total_leish = list(common_area['Accession'])

tcruzis_single = []
tcruzi_sp = []
tcruzis_leish = []
tcruzis_og = []

for i in range(len(common_area))[:]:
    
    line = common_area.iloc[i]
    
    tcruzi = line['Target ID']
    leish = line['Accession'].replace('linf|','')
    
    ortho_group = line['Group ID']
    
    
    tcruzi_list = tcruzi.split(',')
    
    for t in tcruzi_list:
        if '|' in t:
            tcruzi_sp.append(org_dict[t.split('|')[0]])
            tcruzis_single.append(t.split('|')[1])
            tcruzis_leish.append(leish)
            tcruzis_og.append(ortho_group)




clean_common = pd.DataFrame()

clean_common['Linf'] = tcruzis_leish
clean_common['Tcruzi'] = tcruzis_single
clean_common['Organism'] = tcruzi_sp
clean_common['OrthoGroup'] = tcruzis_og



tcruzisp_ref = 'Trypanosoma cruzi CL Brener Esmeraldo-like'

clean_common['is_tryp_ref'] = clean_common['Organism'] == tcruzisp_ref
clean_common_ref = clean_common.query('is_tryp_ref').copy()
clean_common_ref.pop('is_tryp_ref')
clean_common_ref.pop('Organism')



clean_common_ref.to_csv('../results/common/Step1.csv')

print('1. OrthoMDL')
print(len(clean_common_ref),'Leish-Tcruzi pairs')
print(tcruzisp_ref)
for col in clean_common_ref.columns:
    print('...',len(np.unique(clean_common_ref[col])),'unique',col)



after_step1 = []
search_fasta_file = '../support/blast/ProteinsLinf.fasta'

print()
print('Writing FASTA')
#https://biopython.org/wiki/SeqIO
for record in SeqIO.parse(search_fasta_file, "fasta"):
    # sequence_dict[record.id] = str(record.seq)
    record_woT1 = record.id.replace('-T1','')
    if record_woT1 in np.unique(clean_common_ref['Linf']):       
        after_step1.append(record)
    
print('...Going to the next step with', len(after_step1), 'targets.')    
SeqIO.write(after_step1,"../support/blast/search_proteins.fasta","fasta")    
    
    
# SeqIO.write(after_step1,"./bla","fasta")    



