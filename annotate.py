#!/usr/bin/env python
from __future__ import print_function
import pyximport; pyximport.install()
import remap_blast

from itertools import izip
import json
import os
from pprint import pprint
import sys

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task
import jinja2
import pandas as pd
import screed

from tasks import *
import blasttools
from blasttools import outfmt6
import hmmertools

@create_task_object
def aggregate_annotations_task(assembly_fn, blast_targets, transdecoder_fn, hmmscan_fn,
                               sample_df, annotation_fn):

    def get_store():
        return pd.HDFStore(annotation_fn, complib='zlib', complevel=5)

    def get_assembly_df():
        print('--- Building initial assembly DataFrame...')
        store = get_store()
        idx = []
        lens = []
        for record in screed.open(assembly_fn):
            idx.append(record.name)
            lens.append(len(record.sequence))
        annot_df = pd.DataFrame({'length': lens}, index=idx)
        print(annot_df.head())
        store['annot_df'] = annot_df
        store.close()

    best_hits_files = []
    for db_fn in blast_targets:
        target = '{}.x.{}.db.tsv'.format(assembly_fn, db_fn)
        name = 'x_{}'.format(strip_seq_extension(db_fn))
        best_hits_files.append((name, target))

    def get_best_hits():
        print('--- Calculating best hits...')
        store = get_store()
        for i, (name, target) in enumerate(best_hits_files):
            print('\t...working on {}'.format(target))
            
            df = blasttools.blast_to_df(target)
            remap_blast.remap_blast_coords_df(df)

            #store[target] = df

            tmp = pd.merge(pd.DataFrame(index=tpm_df.index), df,
                              left_index=True, right_index=True, how='left')
            blasttools.best_hits(tmp)

            if i == 0:
                best_hits_panel = pd.Panel({name: tmp})
            else:
                best_hits_panel[name] = tmp
        store['best_hits_panel'] = best_hits_panel
        best_hits = best_hits_panel.minor_xs('sseqid')
        store['best_hits_df'] = best_hits
        #store['blast_filter_df'] = best_hits.minor_xs('evalue') >= 0
        store.close()

    ortho_files = []
    for db_fn in blast_targets:
        A_fn = '{}.x.{}.db.tsv'.format(assembly_fn, db_fn)
        B_fn = '{}.db.x.{}.tsv'.format(db_fn, assembly_fn)
        name = 'x_'.format(strip_seq_extension(db_fn))
        ortho_files.append((name, A_fn, B_fn))

    def get_orthologies():
        print('--- Calculating recipricol best hits...')
        store = get_store()
        for i, (name, A_fn, B_fn) in enumerate(ortho_files):
            print('\t...working on {}, {}'.format(A_fn, B_fn))
            
            A = blasttools.blast_to_df(A_fn)
            B = blasttools.blast_to_df(B_fn)

            remap_blast.remap_blast_coords_df(A)
            remap_blast.remap_blast_coords_df(B)
            
            X = blasttools.get_orthologies(A, B, tpm_df.index)
            
            if i == 0:
                ortho_panel = pd.Panel({name: X})
            else:
                ortho_panel[name] = X
        store['ortho_panel'] = ortho_panel
        orthos = ortho_panel.minor_xs('sseqid')
        store['ortho_df'] = orthos
        #store['ortho_filter_df'] = orthos.minor_xs('evalue') >= 0
        store.close()

    def get_transdecoder():
        print('--- Parsing TransDecoder results...')
        store = get_store()
        transdecoder_df = hmmertools.gff3_transdecoder_to_df(transdecoder_fn)
        store['transdecoder_dr'] = transdecoder_df
        store.close()
    

    def get_pfam():
        print('--- Parsing HMMER P-fam results...')
        store = get_store()
        pfam_df = hmmertools.hmmscan_to_df(hmmscan_fn)
        store['pfam_df'] = pfam_df
        store.close()

    def get_tpm():
        print('--- Getting expression data...')
        store = get_store()
        tpm_df = pd.read_csv(tpm_fn, delimiter='\t', index_col=0)
        labels = dict(zip(sample_df.filename, sample_df.label))
        tpm_df.rename(columns=labels, inplace=True)
        tpm_df.sort(axis=1, inplace=True)
        store['tpm_df'] = tpm_df
        store.close()

    actions = [(get_assembly_df, []), (get_best_hits, []), (get_orthologies, []),
               (get_transdecoder, []), (get_pfam, []), (get_tpm, [])]

    blast_file_dep = []
    for _, fn in best_hits_files:
        blast_file_dep.append(fn)
    for _, A_fn, B_fn in ortho_files:
        blast_file_dep.extend([A_fn, B_fn])
    file_dep = [assembly_fn, transdecoder_fn, hmmscan_fn] + blast_file_dep

    return {'title': title_with_actions,
            'name': 'aggregate_annotations_' + assembly_fn,
            'actions': actions,
            'file_dep': file_dep,
            'targets': [annotation_fn],
            'clean': [clean_targets]}
