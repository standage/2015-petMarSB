#!/usr/bin/env python

from joshua.intervaltree import Interval, IntervalTree
import pandas as pd
import pyprind

class DataFrameInterval(Interval):

    def __init__(self, idx, df, *args, **kwargs):
        Interval.__init__(self, *args, **kwargs)
        self.idx = idx
        self.df = df

    @property
    def value(self):
        return self.df.loc[self.idx]

    def __str__(self):
        return 'DataFrameInterval({s}, {e}, idx={idx}, '\
               'chrom={chrom}, strand={strand})'.format(s=self.start, e=self.end, idx=self.idx,
                                                                          chrom=self.chrom, strand=self.strand)


def get_gtf_coords(row):
    '''
    Extract the coordinates and Interval information from a GTF series
    '''

    return row.start, row.end, row.contig_id, row.strand

def get_blast_subject_coords(row):
    '''
    Extract the coordinates and Interval information for the subject in a BLAST hit
    '''

    return row.sstart, row.send, row.qseqid, row.sstrand

def get_blast_query_coords(row):
    '''
    Extract the coordinates and Interval information for the query in a BLAST hit
    '''

    return row.qstart, row.qend, row.sseqid, row.qstrand

def build_tree_from_group(group, ref_dataframe, coord_func, bar=None):
    '''
    Build an IntervalTree from a sub-DataFrame, for groupby() operations
    '''
    tree = IntervalTree()

    for idx, row in group.iterrows():
        if bar is not None:
            bar.update()
        start, end, chrom, strand = coord_func(row)
        tree.insert_interval(DataFrameInterval(idx, ref_dataframe, start, end,
                                               chrom=chrom, strand=strand))

    return tree

def merge_overlapping(iv_list):
    '''
    Merge all overlapping Intervals from the given list, and return a new list
    with the collapsed Intervals
    '''
    if len(iv_list) == 1:
        return iv_list

    iv_list = sorted(iv_list, key=lambda iv: iv.start)
    merge_stack = [iv_list[0]]
    for iv in iv_list[1:]:
        top = merge_stack[-1]
        if top.end < iv.start:
            merge_stack.append(iv)
        elif top.end < iv.end:
            top.end = iv.end
            merge_stack.pop()
            merge_stack.append(top)
    return merge_stack

def calc_bases_overlapped(iv, overlap_ivs):
    '''
    Given a target Interval and a list of *overlapping* Intervals,
    calculate how many bases in the target are overlapped
    '''
    merged = merge_overlapping(overlap_ivs)

    covered = 0
    for overlap_iv in merged:
        if overlap_iv.start <= iv.start:
            if overlap_iv.end <= iv.end:
                covered += overlap_iv.end - iv.start
            else:
                covered = iv.end - iv.start
                break
        else:
            if overlap_iv.end <= iv.end:
                covered += overlap_iv.end - overlap_iv.start
            else:
                covered = iv.end - overlap_iv.start
                break
    assert covered <= len(iv)
    return covered


def tree_intersect(tree_A, tree_B, cutoff=0.9):
    '''
    Find all overlaps of Intervals in A by B and return a `dict` of the results,
    where keys are (annotation_id, alignment_id) tuples (with the id being the
    interger index in the original DataFrames) and the values are the number
    of bases overlapped in A.
    '''

    if type(tree_A) is not IntervalTree:
        raise TypeError('tree_A must be a valid IntervalTree (got {t})'.format(t=type(tree_A)))
    if type(tree_B) is not IntervalTree:
        return {}

    overlaps = {}
    def overlap_fn(node):
        iv = node.interval
        for ov in tree_B.find(iv.start, iv.end):
            ov_len = calc_bases_overlapped(iv, [ov])
            assert ov_len <= (iv.end - iv.start)
            if (float(ov_len) / len(iv)) >= cutoff:
                overlaps[(iv.idx, ov.idx)] = ov_len
    tree_A.traverse(overlap_fn)
    return overlaps


def tree_coverage_intersect(tree_A, tree_B, cutoff=0.9):
    '''
    Like `tree_intersect`, but merges overlapping interavals in `tree_B` to
    calculate the overlap length.    
    '''

    if type(tree_A) is not IntervalTree:
        raise TypeError('tree_A must be a valid IntervalTree (got {t})'.format(t=type(tree_A)))
    if type(tree_B) is not IntervalTree:
        return {}

    overlaps = {}
    def overlap_fn(node):
        iv = node.interval
        ov_list = tree_B.find(iv.start, iv.end)
        if ov_list:
            ov_len = calc_bases_overlapped(iv, ov_list)
            assert ov_len <= (iv.end - iv.start)
            if (float(ov_len) / len(iv)) >= cutoff:
                overlaps[iv.idx] = ov_len
    tree_A.traverse(overlap_fn)
    return overlaps


def get_ann_aln_overlap_df(tree_df, cutoff=0.9, merge=False, bar=None):
    '''
    Perform the tree intersect between all pairs of (annotation, alignment) IntervalTrees
    in the given DataFrame
    '''

    data = []
    if merge:
        for contig_id, ann_tree, aln_tree in tree_df.itertuples():
            if type(ann_tree) is IntervalTree: # as opposed to NaN
                d = tree_coverage_intersect(ann_tree, aln_tree, cutoff=cutoff)
                if d:
                    data.append(pd.DataFrame({'overlap_len': d.values()}, index=d.keys()))
            if bar:
                bar.update()
    else:
        for contig_id, ann_tree, aln_tree in tree_df.itertuples():
            if type(ann_tree) is IntervalTree:
                d = tree_intersect(ann_tree, aln_tree, cutoff=cutoff)
                if d:
                    index = pd.MultiIndex.from_tuples(d.keys(), names=['ann_id', 'aln_id'])
                    data.append(pd.DataFrame({'overlap_len': d.values()}, index=index))
            if bar is not None:
                bar.update()

    return pd.concat(data, axis=0)


def get_aln_ann_overlap_df(tree_df, cutoff=0.9, merge=False, bar=None):
    '''
    Perform the tree intersect between all pairs of (alignment, annotation) IntervalTrees
    in the given DataFrame
    '''

    data = []
    if merge:
        for contig_id, ann_tree, aln_tree in tree_df.itertuples():
            if type(aln_tree) is IntervalTree: # as opposed to NaN
                d = tree_coverage_intersect(aln_tree, ann_tree, cutoff=cutoff)
                if d:
                    data.append(pd.DataFrame({'overlap_len': d.values()}, index=d.keys()))
            if bar:
                bar.update()
    else:
        for contig_id, ann_tree, aln_tree in tree_df.itertuples():
            if type(aln_tree) is IntervalTree:
                d = tree_intersect(aln_tree, ann_tree, cutoff=cutoff)
                if d:
                    index = pd.MultiIndex.from_tuples(d.keys(), names=['aln_id', 'ann_id'])
                    data.append(pd.DataFrame({'overlap_len': d.values()}, index=index))
            if bar is not None:
                bar.update()

    return pd.concat(data, axis=0)

def check_ann_covered_single(ann_df, aln_tree_df, cutoff=0.90):

    def overlap_fn(row):
        iv = Interval(row.start, row.end)
        try:
            overlaps = aln_tree_df.loc[row.contig_id].find(iv.start, iv.end)
        except KeyError:
            return False
        for ov in overlaps:
            covered = calc_bases_overlapped(iv, [ov])
            if float(covered) / (row.end - row.start) >= cutoff:
                return True
        return False

    return ann_df.apply(overlap_fn, axis=1)


def check_ann_covered(ann_df, aln_tree_df, cutoff=0.90):

    def overlap_fn(group):
        contig_id = group.contig_id.iloc[0]
        try:
            tree = aln_tree_df.loc[contig_id]
        except KeyError:
            return pd.Series([False] * len(group), index=group.index)

        results = []
        for key, row in group.iterrows():
            result = False
            iv = Interval(row.start, row.end)
            overlaps = tree.find(iv.start, iv.end)
            if overlaps:
                covered = calc_bases_overlapped(iv, overlaps)
                if float(covered) / (row.end - row.start) >= cutoff:
                    result = True
            results.append(result)
        
        return pd.Series(results, index=group.index)

    result_df = ann_df.groupby('contig_id').apply(overlap_fn)
    result_df.index = result_df.index.droplevel(0)
    return result_df
