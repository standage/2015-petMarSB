#!/usr/bin/env python

import csv
import sys
import pandas as pd

outfmt6 = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

def blast_to_df(fn, names=outfmt6, delimiter='\t', index_col=0):
    raw_df = pd.read_table(fn, header=None, index_col=index_col, names=names,
                             delimiter=delimiter, skipinitialspace=True)
    
    def qstrand_fn(row):
        if row.qstart > row.qend:
            return '-'
        return '+'

    def sstrand_fn(row):
        if row.sstart > row.send:
            return '-'
        return '+'

    def sstart_fn(row):
        if row.sstart < row.send:
            return row.sstart - 1
        else:
            return row.send

    def send_fn(row):
        if row.sstart < row.send:
            return row.send
        else:
            return row.sstart + 1
 
    def qstart_fn(row):       
        if row.qstart < row.qend:
            return row.qstart - 1
        else:
            return row.qend

    def qend_fn(row):
        if row.qstart < row.qend:
            return row.qend
        else:
            return row.qstart + 1

    raw_df['qstrand'] = raw_df.apply(qstrand_fn, axis=1)
    raw_df['sstrand'] = raw_df.apply(sstrand_fn, axis=1)

    qstart = raw_df.apply(qstart_fn, axis=1)
    qend = raw_df.apply(qend_fn, axis=1)
    sstart = raw_df.apply(sstart_fn, axis=1)
    send = raw_df.apply(send_fn, axis=1)

    raw_df['qstart'] = qstart
    raw_df['qend'] = qend
    raw_df['sstart'] = sstart
    raw_df['send'] = send

    return raw_df

def best_hits(df, comp_col='evalue'):
    assert type(df) is pd.DataFrame
    df['ind'] = df.index
    df.sort(columns=['ind', comp_col], inplace=True)
    df.drop_duplicates(subset='ind', inplace=True)
    del df['ind']

def get_orthologies(A, B, idx, comp_col='evalue'):

    assert type(A) is pd.DataFrame
    assert type(B) is pd.DataFrame
    assert type(idx) is pd.Index

    best_hits(A)
    best_hits(B)

    X = pd.merge(A, B, how='inner', left_on='sseqid', right_index=True)
    X = pd.merge(pd.DataFrame(index=idx),
                      X[(X.index == X.sseqid_y)], left_index=True, right_index=True, how='left')
    del X['sseqid_x']
    del X['sseqid_y']

    return X
