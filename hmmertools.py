#!/usr/bin/env python

import pandas as pd

hmmscan_cols = ['target_name', 'target_accession', 'tlen', 'query_name', 
        'query_accession', 'query_len', 'full_evalue', 'full_score', 
        'full_bias', 'domain_num', 'domain_total', 'domain_c_value', 
        'domain_i_evalue', 'domain_score', 'domain_bias', 'hmm_coord_from', 
        'hmm_coord_to', 'ali_coord_from', 'ali_coord_to', 'env_coord_from', 
        'env_coord_to', 'accuracy', 'description']

gff3_transdecoder_cols = ['seqid', 'feature_type', 'start', 'end', 'strand']

def hmmscan_to_df(fn):
    def split_query(item):
        q, _, _ = item.partition('|')
        return q
    data = []
    with open(fn) as fp:
        for ln in fp:
            if ln.startswith('#'):
                continue
            tokens = ln.split()
            data.append(tokens[:len(hmmscan_cols)-1] + [' '.join(tokens[len(hmmscan_cols)-1:])])
    df = pd.DataFrame(data, columns=hmmscan_cols)
    df.query_name = df.query_name.apply(split_query)
    df.set_index('query_name', inplace=True)
    return df

def gff3_transdecoder_to_df(fn):
    data = []
    with open(fn) as fp:
        for ln in fp:
            if ln == '\n':
                continue
            tokens = ln.split('\t')
            try:
                data.append([tokens[0]] + tokens[2:5] + [tokens[6]])
            except IndexError as e:
                print e
                print tokens
                break
    df = pd.DataFrame(data, columns=gff3_transdecoder_cols)
    df.set_index('seqid', inplace=True)
    return df

