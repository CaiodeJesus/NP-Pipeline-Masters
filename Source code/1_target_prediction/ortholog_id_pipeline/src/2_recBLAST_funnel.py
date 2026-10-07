#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:45:33 2026

@author: caio
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import timeit

def filter_blast(raw, name, metric, threshold):
    print(name)
    print('Before filter')
    print('\tQuery (search) size:', len(np.unique(raw['Query id'])))
    print('\tSubject (DB) size:', len(np.unique(raw['Subject id'])))
    print()

    raw['do_pass'] = raw[metric] > threshold
    filtered = raw.query('do_pass').copy()
    
    print('After filter')
    print('\tQuery (search) size:', len(np.unique(filtered['Query id'])))
    print('\tSubject (DB) size:', len(np.unique(filtered['Subject id'])))
    print()

    
    homology_pairs = pd.DataFrame()
    homology_pairs['Query'] = filtered['Query id']
    homology_pairs['Subject'] = filtered['Subject id']

    return homology_pairs



"""
(base) caio@lqmc09:~/2026/ortholog_id_pipeline/funnel/support/blast$ pwd
/home/caio/2026/ortholog_id_pipeline/funnel/support/blast
(base) caio@lqmc09:~/2026/ortholog_id_pipeline/funnel/support/blast$ makeblastdb -in ProteinsLinf.fasta -title Linf -dbtype prot -out LinfDB/Linf -parse_seqids


Building a new DB, current time: 03/19/2026 15:35:36
New DB name:   /home/caio/2026/ortholog_id_pipeline/funnel/support/blast/LinfDB/Linf
New DB title:  Linf
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 8527 sequences in 0.23746 seconds.


(base) caio@lqmc09:~/2026/ortholog_id_pipeline/funnel/support/blast$ makeblastdb -in ProteinsTcruzi.fasta -title Tcruzi -dbtype prot -out TcruziDB/Tcruzi -parse_seqids


Building a new DB, current time: 03/19/2026 15:35:46
New DB name:   /home/caio/2026/ortholog_id_pipeline/funnel/support/blast/TcruziDB/Tcruzi
New DB title:  Tcruzi
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 9039 sequences in 0.296173 seconds.


(base) caio@lqmc09:~/2026/ortholog_id_pipeline/funnel/support/blast$ blastp -query search_proteins.fasta -db TcruziDB/Tcruzi -out query_after1.txt -outfmt 6
contra o proteoma de Tcruzi
"""


with open('../support/blast/input_blast_filter', 'r') as f:
    metric, threshold = f.readlines()[1].split('\t')
threshold = int(threshold)

columns = ["Query id",
           "Subject id",
           "% identity",
           "alignment length",
           "mismatches",
           "gap openings",
           "query start",
           "query end",
           "subject start",
           "subject end",
           "e-value",
           "bit score"]
#https://conmeehan.github.io/blast+tutorial.html

raw_query_after1 = pd.read_csv('../support/blast/query_after1.txt',
                          sep = '\t', header = None,
                          names = columns,
                          # index_col = 'Query id'
                          )
after1_in_tcruzi = filter_blast(raw_query_after1,
                                'Step1 -> Tryp', metric, threshold)['Subject']
after1_in_tcruzi = np.unique(after1_in_tcruzi)


fasta_after1_in_tcruzi = []
search_fasta_file = '../support/blast/ProteinsTcruzi.fasta'
#https://biopython.org/wiki/SeqIO
for record in SeqIO.parse(search_fasta_file, "fasta"):
    # sequence_dict[record.id] = str(record.seq)
    
    if record.id in after1_in_tcruzi:       
        fasta_after1_in_tcruzi.append(record)
    
    
SeqIO.write(fasta_after1_in_tcruzi,"../support/blast/After1toTcruzi.fasta","fasta")    


"""
(base) caio@lqmc09:~/2026/ortholog_id_pipeline/funnel/support/blast$ blastp -query After1toTcruzi.fasta -db LinfDB/Linf -out query_tryp.txt -outfmt 6
contra o proteoma de Linf
"""




raw_query_tryp = pd.read_csv('../support/blast/query_tryp.txt',
                          sep = '\t', header = None,
                          names = columns,
                          # index_col = 'Query id'
                          )
reciprocal = filter_blast(raw_query_tryp,
                          'Tcruzi back to Linf', metric, threshold)['Subject']
reciprocal = np.unique(reciprocal)


common_last_step = pd.read_csv('../results/common/Step1.csv')

is_in_next = np.empty(len(common_last_step))
i=0
for linf in common_last_step['Linf']:
    is_in_next[i] = linf+'-T1' in reciprocal
    i+=1


common_last_step['is_in_next'] = is_in_next.astype(bool)

next_step = common_last_step.query('is_in_next').copy()
next_step.pop('is_in_next')


print('2. Reciprocal BLAST')
print(len(next_step),'Linf-Tryp pairs in original TP->Linf set')
for col in next_step.columns:
    print('...',len(np.unique(next_step[col])),'unique',col)

next_step.to_csv('../results/common/Step2.csv')




id_list = '\n'.join(np.unique(next_step['Linf']))
out_file = '../results/commonIDListAfterStep2.txt'
with open(out_file,'w') as f:
    f.write(id_list)
    
print()
print('Saving list of ids to',out_file)
