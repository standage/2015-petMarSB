#!/usr/bin/env python

from bx.intervals.intersection import Interval, IntervalTree
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

    def __len__(self):
        return self.end - self.start

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

    if row.sstart < row.send:
        start = row.sstart - 1
        end = row.send
    else:
        start = row.send
        end = row.sstart + 1
    strand = '+'
    if row.sstart > row.send:
        strand = '-'
    return start, end, row.qseqid, strand

def get_blast_query_coords(row):
    '''
    Extract the coordinates and Interval information for the query in a BLAST hit
    '''

    if row.qstart < row.qend:
        start = row.qstart - 1
        end = row.qend
    else:
        start = row.qend
        end = row.qstart + 1
    strand = '+'
    if row.qstart > row.qend:
        strand = '-'
    return start, end, row.sseqid, strand

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
    return covered

def coverage_intersect(tree_A, tree_B):
    '''
    Perform a true intersection A^B: for each interval in A, find its length covered by intervals in B.
    '''

    if type(tree_A) is not IntervalTree:
        raise TypeError('tree_A must be a valid IntervalTree')

    overlaps = []
    def overlap_fn(node):
        iv = node.interval
        if type(tree_B) is IntervalTree:
            overlaps.append( (iv, tree_B.find(node.interval.start, node.interval.end)) )
        else:
            overlaps.append( (iv, []) )

    covered_data = {}
    tree_A.traverse(overlap_fn)
    for iv, overlap_list in overlaps:
        covered = 0
        if overlap_list:
            covered = calc_bases_overlapped(iv, overlap_list)
        assert covered <= len(iv)
        covered_data[iv.idx] = pd.Series({'covered': covered, 'length': len(iv)})
    return pd.DataFrame(covered_data).T

def get_gtf_aln_coverage(tree_df):
    '''
    Get the coverage intersect all pairs of (annotation, alignment) IntervalTrees
    in the given DataFrame; return as a new, larger DataFrame
    '''

    data = []
    for contig_id, ann_tree, aln_tree in tree_df.itertuples():
        if type(ann_tree) is IntervalTree:
            d = coverage_intersect(ann_tree, aln_tree)
            data.append(d)
    return pd.concat(data, axis=0)

def tree_intersect(tree_A, tree_B):
    '''
    Find all overlaps of Intervals in A by B and return a `dict` of the results,
    where keys are (annotation_id, alignment_id) tuples (with the id being the
    interger index in the original DataFrames) and the values are the number
    of bases overlapped in A.
    '''

    if type(tree_A) is not IntervalTree:
        raise TypeError('tree_A must be a valid IntervalTree')

    overlaps = {}
    def overlap_fn(node):
        iv = node.interval
        if type(tree_B) is IntervalTree:
            for ov in tree_B.find(node.interval.start, node.interval.end):
                ov_len = calc_bases_overlapped(iv, [ov])
                assert ov_len <= (node.interval.end - node.interval.start)
                overlaps[(iv.idx, ov.idx)] = ov_len
    tree_A.traverse(overlap_fn)
    return overlaps

def get_gtf_aln_overlap_df(tree_df, bar=None):
    '''
    Perform the tree intersect between all pairs of (annotation, alignment) IntervalTrees
    in the given DataFrame; return as a new, larger DataFrame
    '''

    data = []
    for contig_id, ann_tree, aln_tree in tree_df.itertuples():
        if type(ann_tree) is IntervalTree:
            d = tree_intersect(ann_tree, aln_tree)
            if d:
                index = pd.MultiIndex.from_tuples(d.keys(), names=['ann_id', 'aln_id'])
                data.append(pd.DataFrame({'overlap_len': d.values()}, index=index))
        if bar is not None:
            bar.update()

    return pd.concat(data, axis=0)

def check_ann_covered(ann_df, aln_tree_df, cutoff=0.90):

    def overlap_fn(row):
        overlaps = aln_tree_df.loc[row.
