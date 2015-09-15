
# coding: utf-8

## Lamprey Transcriptome Analysis: Gene Model Overlap Notebook

# ```
# Camille Scott [camille dot scott dot w @gmail.com] [@camille_codon]
# 
# camillescott.github.io
# 
# Lab for Genomics, Evolution, and Development
# Michigan State University
# ```

### About

# Explores consensus between assembled transcripts and existing lamprey gene models, using the `gtf_to_genes` library.

#     assembly version: lamp03
# 
#     assembly program: Trinity
#     
#     gtf model versin: 7.0.78, at ftp://ftp.ensembl.org/pub/release-78/gtf/petromyzon_marinus/Petromyzon_marinus.Pmarinus_7.0.78.gtf.gz

# In[1]:

get_ipython().system(u'echo $PWD')


### Contents

# 1. [Libraries](#Libraries) 
# * [Metadata](#Metadata) 
#      1. [Databases](#Databases)
#      * [Samples](#Samples)
# * [Data Loading](#Data-Loading)
# * [Matplotlib RC Settings](#Matplotlib-RC-Settings)

### Libraries

# In[1]:

get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')
from petmar_common import *


# In[2]:

from IPython.html.widgets import interactive, RadioButtonsWidget
from IPython.display import display
from IPython.html import widgets


# In[35]:

import gffutils


# In[22]:

def wdir(fn):
    return os.path.join('../../work/', fn)


# In[4]:

pd.set_option('display.max_rows', 200)


# In[5]:

store = pd.HDFStore(wdir('{}.store.h5'.format(metadata['prefix'])), complib='zlib', complevel=9)


# In[6]:

prefix = metadata['prefix'] + '_'


# In[7]:

db_metadata = store['db_metadata']


# In[8]:

sample_metadata = store['sample_metadata']


# In[23]:

gtf_file = wdir('Petromyzon_marinus.Pmarinus_7.0.78.gtf')


### Data Loading

# In[9]:

assem_db = screed.ScreedDB(wdir(db_metadata.fn.assembly))


# In[10]:

blast_panel = store['blast_panel']
blast_df = store['blast_df']


# In[16]:

get_ipython().magic(u'config InlineBackend.close_figures = False')


# In[38]:

#gtf_db = gffutils.create_db(gtf_file, dbfn=prefix+'gtf.db', force=True, keep_order=False,
#                            merge_strategy='merge', sort_attribute_values=False)
print 'load db'
gtf_db = gffutils.FeatureDB(prefix+'gtf.db')


# In[97]:

pM_df = blast_panel['petMar2']


# In[107]:

pM_df['gff'] = None


# In[151]:

def query_func(row):
    if type(row.sseqid) is str:
        return [f.astuple() for f in gtf_db.region('{}:{}-{}'.format(row.sseqid, int(row.sstart), int(row.send)), 
                                         completely_within=False)]
    return None


# In[ ]:
print 'do apply'
pM_df['gff'] = pM_df.apply(query_func, axis=1)

print 'done with apply, store in hdf5'
store['pM_df'] = pM_df
store.close()

