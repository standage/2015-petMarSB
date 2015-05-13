#!/usr/bin/env python

'''
Libraries, functions, settings common to all petmar notebooks.
'''

###########
# Imports #
###########

from itertools import izip
import os
   
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles

import seaborn as sns

import screed

from peasoup import configtools
from peasoup.plot import FigManager
from peasoup.plot import plot_dendro

from collections import defaultdict as ddict

import pandas as pd
from pandas import Categorical


import numpy as np
from numpy.random import rand

from scipy import stats 
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch

from sklearn.decomposition import PCA, KernelPCA
from sklearn.cluster import KMeans

from IPython.display import SVG
from IPython.display import FileLink
from IPython.display import HTML

import mygene
from bioservices import UniProt

from pprint import pprint

import warnings

#######################
# Functions, settings #
#######################

pd.set_option('display.max_rows', 200)

def wdir(fn=''):
    return os.path.join('../../work/', fn)

outfmt6 = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', \
           'qend', 'sstart', 'send', 'evalue', 'bitscore']

metadata = configtools.get_cfg('../../metadata.ini', '../../metadata.spec.ini')

warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings("ignore",category=UserWarning)

def uniprot_str(uid):
    tmp = uid.split('|')
    return tmp[1]

def filter_species(gene_set, mg, species='7757'):
    if mg is None:
        mg = mygene.MyGeneInfo()
    results = set()
    for result in mg.querymany(list(gene_set), scopes='symbol,name,ensemblgene,entrezgene,uniprot', species=species):
        if 'notfound' in result:
            results.add(result['query'])
    return list(results)

def gene_from_symbol(gene_list, mg):
    if mg is None:
        mg = mygene.MyGeneInfo()
    results = mg.querymany(gene_list, scopes='symbol,name,alias', 
                            fields='name,alias,summary,go,refseq',
                            species='all', as_dataframe=True, returnall=False)
    return results

def dedupe_gene_func(X):
    pos = np.array(X.isnull().sum(axis=1)).argmin()
    return X.ix[pos]

def transcript_write_func(X, gene_map, db, outfp):
    gene_tag = X.tag
    tr_names = gene_map[gene_tag]
    gene_alias = X.alias_str
    gene_name = X.name
    
    if type(gene_tag) == float and np.isnan(gene_tag):
        gene_tag = ''
    else:
        gene_tag = 'tag={}'.format(gene_tag)
    if type(gene_alias) == float and np.isnan(gene_alias):
        gene_alias = ''
    else:
        gene_alias = 'alias={}'.format(gene_alias)
    if type(gene_name) == float and np.isnan(gene_name):
        gene_name = ''
    else:
        gene_name = 'name={}'.format(gene_name)

    for tr_name in tr_names:
        h_line = '{0} {1} {2} {3}'.format(tr_name, gene_tag, gene_alias, gene_name)
        tr = db[tr_name].sequence
        outfp.write('>{}\n{}\n'.format(h_line, tr))


