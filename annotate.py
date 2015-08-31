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
import hmmertools

def annotate_task(assembly_fn, blast_targets, transdecoder_fn, hmmscan_fn,
                  busco_fn, annotation_fn):

    def get_store():
        return pd.HDFStore(annotation_fn, complib='zlib', complevel=5)

    def annotate():
        store = get_store()
        idx = []
        lens = []
        for record in screed.open(assembly_fn):
            idx.append(record.name)
            lens.append(len(record.sequence))
        annot_df = pd.DataFrame({'length': lens}, index=idx)
        print(annot_df.head())

        print('--- Calculating best hits...')
        for i, db_fn in enumerate(blast_targets):
            target = '{}.x.{}.db.tsv'.format(assembly_fn, db_fn)
            name = 'x_{}'.format(strip_seq_extension(db_fn))
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

        print('--- Calculating recipricol best hits...')
        for i, db_fn in enumerate(blast_targets):
            A_fn = '{}.x.{}.db.tsv'.format(assembly_fn, db_fn)
            B_fn = '{}.db.x.{}.tsv'.format(db_fn, assembly_fn)
            name = 'x_'.format(strip_seq_extension(db_fn))
            print('\t...working on {}, {}'.format(A_fn, B_fn))
            
            A = pd.read_table(A_fn, header=None, index_col=0, names=outfmt6)
            B = pd.read_table(B_fn, header=None, index_col=0, names=outfmt6)

            fix_blast_coords_df(A)
            fix_blast_coords_df(B)
            
            X = blasttools.get_orthologies(A, B, tpm_df.index)
            
            if i == 0:
                ortho_panel = pd.Panel({name: X})
            else:
                ortho_panel[name] = X

        print('--- Parsing TransDecoder results...')
        transdecoder_df = hmmertools.gff3_transdecoder_to_df(transdecoder_fn)

        print('--- Parsing HMMER P-fam results...')
        pfam_df = hmmertools.hmmscan_to_df(hmmscan_fn)

        print('--- Aggregating...')

