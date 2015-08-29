#!/usr/bin/env python

import pandas as pd

cols = ['target_name', 'target_accession', 'tlen', 'query_name', 
        'query_accession', 'query_len', 'full_evalue', 'full_score', 
        'full_bias', 'domain_num', 'domain_total', 'domain_c_value', 
        'domain_i_evalue', 'domain_score', 'domain_bias', 'hmm_coord_from', 
        'hmm_coord_to', 'ali_coord_from', 'ali_coord_to', 'env_coord_from', 
        'env_coord_to', 'accuracy', 'description']

def hmmscan_to_df(fn):

    data = []
    with open('_work/lamp10.fasta.pfam-A.out') as fp:
        for ln in fp:
            tokens = ln.split()
            data.append(tokens[:len(cols)-1] + [' '.join(tokens[len(cols)-1:])])
    return pd.DataFrame(data[3:], columns=cols)
