#!/usr/bin/env python

import os
import sys
import pandas as pd

def parse_busco(fn):
    res = {}
    with open(fn) as fp:
        for ln in fp:
            if ln.strip().startswith('C:'):
                tokens = ln.split(',')
                for token in tokens:
                    key, _, val = token.partition(':')
                    key = key.strip()
                    val = val.strip().strip('%')
                    if key == 'C':
                        valc, _, vald = val.partition('%')
                        valc = valc.strip()
                        vald = vald.strip('D:][%')
                        res['C(%)'] = valc
                        res['D(%)'] = vald
                    else:
                        if key != 'n':
                           key += '(%)'
                        res[key] = val.strip().strip('%')
    return res

def busco_to_df(fn_list, dbs=['metazoa', 'vertebrata']):

    data = []
    for fn in fn_list:
        data.append(parse_busco(fn))

    df = pd.DataFrame(data)
    df['fn'] = [os.path.basename(fn)[14:-14].strip('.') for fn in fn_list]
    df['db'] = None
    for db in dbs:
        idx = df.fn.str.contains(db)
        df.loc[idx,'db'] = db
        df.loc[idx,'fn'] = df.loc[idx, 'fn'].apply(lambda fn: fn[:fn.find(db)].strip('. '))
    return df

def formatted(df, dbs=['metazoa', 'vertebrata']):

    df.set_index('fn', inplace=True)
    columns = pd.MultiIndex.from_product([dbs, ['C(%)', 'D(%)', 'F(%)', 'M(%)', 'n']])
    separated = []
    for db in dbs:
        separated.append(df[df.db == db])
        del separated[-1]['db']
    combined = pd.concat(separated, axis=1)
    combined.columns = columns
    return combined

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--busco-files', nargs='+')
    parser.add_argument('-o', dest='output', type=argparse.FileType('wb'), default=sys.stdout)
    args = parser.parse_args()

    df = busco_to_df(args.busco_files)
    df = df.pivot('fn', 'db').stack()
    df.reset_index().to_json(args.output)
