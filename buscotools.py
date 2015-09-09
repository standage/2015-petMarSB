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
                    res[key.strip()] = val.strip()
    return res

def busco_to_df(fn_list):

    data = []
    for fn in fn_list:
        data.append(parse_busco(fn))

    return pd.DataFrame(data, index=[os.path.basename(fn)[14:-14] for fn in fn_list])

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--busco-files', nargs='+')
    parser.add_argument('-o', dest='output', type=argparse.FileType('wb'), default=sys.stdout)
    args = parser.parse_args()

    df = busco_to_df(args.busco_files)
    df.to_latex(args.output)
