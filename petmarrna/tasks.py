#!/usr/bin/env python
from __future__ import print_function

from itertools import izip
import json
import os
import pprint
import re
from shutil import rmtree
import shutil
import sys

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task

from bioservices import UniProt
import jinja2
import pandas as pd
import screed

#from peasoup.tasks import BlastTask

def clean_folder(target):
    try:
        rmtree(target)
    except OSError:
        pass

seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]

def create_task_object(task_dict_func):
    '''Wrapper to decorate functions returning pydoit
    Task dictionaries and have them return pydoit Task
    objects
    '''
    def d_to_t(*args, **kwargs):
        ret_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(ret_dict)
    return d_to_t

@create_task_object
def diginorm_task(input_files, dg_cfg, label, ct_outfn=None):

    ksize = dg_cfg['ksize']
    table_size = dg_cfg['table_size']
    n_tables = dg_cfg['n_tables']
    coverage = dg_cfg['coverage']

    name = 'normalize_by_median_' + label
    report_fn = label + '.report.txt'

    inputs = ' '.join(input_files)
    ct_out_str = ''
    if ct_outfn is not None:
        ct_out_str = '-s ' + ct_outfn

    cmd = 'normalize-by-median.py -f -k {ksize} -x {table_size} -N {n_tables} '\
          '-C {coverage} -R {report_fn} {ct_out_str} {inputs}'.format(**locals())

    targets = [fn + '.keep' for fn in input_files]
    targets.append(report_fn)
    if ct_out_str:
        targets.append(ct_outfn)

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'file_dep': input_files,
            'targets': targets,
            'clean': [clean_targets]}

@create_task_object
def filter_abund_task(input_files, ct_file, fab_cfg, label):

    name = 'filter_abund_' + label

    min_abund = fab_cfg['min_abund']
    coverage = fab_cfg['coverage']

    inputs = ' '.join(input_files)
    cmd = 'filter-abund.py -C {min_abund} -V -Z {coverage} '\
          '{ct_file} {inputs}'.format(**locals())

    targets = [fn + '.abundfilt' for fn in input_files]

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'file_dep': input_files + [ct_file],
            'targets': targets,
            'clean': [clean_targets]}

@create_task_object
def download_and_gunzip_task(url, target_fn, label=''):
    cmd = 'curl {url} | gunzip -c > {target_fn}'.format(**locals())

    name = '_'.join(['download_gunzip', target_fn, label])

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [run_once]}

@create_task_object
def download_and_untar_task(url, target_dir, label=''):

    cmd1 = 'mkdir -p {target_dir}; curl {url} | tar -xz -C {target_dir}'.format(**locals())
    name = '_'.join(['download_untar', target_dir.strip('/'), label])
    done = name + '.done'
    cmd2 = 'touch {name}'.format(name=name)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd1, cmd2],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [run_once]}

@create_task_object
def uniprot_query_task(query, target_fn, fmt='fasta', label=''):

    def func():
        u = UniProt()
        res = u.search(query, frmt=fmt)
        with open(target_fn, 'wb') as fp:
            fp.write(res)

    if not label:
        label = query
    name = 'uniprot_query_' + label

    return {'name': name,
            'title': title_with_actions,
            'actions': [(func, [])],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [run_once] }

@create_task_object
def truncate_fasta_header_task(fasta_fn, length=500):

    def func():
        tmp_fn = fasta_fn + '.tmp'
        with open(tmp_fn, 'wb') as fp:
            for r in screed.open(fasta_fn):
                try:
                    name = r.name.split()[0]
                except IndexError:
                    name = r.name[:length]
                fp.write('>{}\n{}\n'.format(name, r.sequence))

        shutil.move(tmp_fn, fasta_fn)

    name = 'truncate_fasta_header_{fasta_fn}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [(func, [])],
            'file_dep': [fasta_fn],
            'uptodate': [run_once]}

@create_task_object
def create_folder_task(folder, label=''):

    if not label:
        label = 'create_folder_{folder}'.format(**locals())

    return {'title': title_with_actions,
            'name': label,
            'actions': [(create_folder, [folder])],
            'targets': [folder],
            'uptodate': [run_once],
            'clean': [clean_targets] }

@create_task_object
def blast_format_task(db_fn, db_out_fn, db_type):
    assert db_type in ['nucl', 'prot']

    cmd = 'makeblastdb -in {db_fn} -dbtype {db_type} -out {db_out_fn}'.format(**locals())

    target_fn = ''
    if db_type == 'nucl':
        target_fn = db_out_fn + '.nhr'
    else:
        target_fn = db_out_fn + '.phr'

    name = 'makeblastdb_{db_out_fn}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd, 'touch '+db_out_fn],
            'targets': [target_fn, db_out_fn],
            'file_dep': [db_fn],
            'clean': [clean_targets, 'rm -f {target_fn}.*'.format(**locals())] }

#def blast_task(row, config, assembly):
#    ''' TODO: Remove peasoup usage
#    '''
#    blast_threads = config['pipeline']['blast']['threads']
#    blast_params = config['pipeline']['blast']['params']
#    blast_evalue = config['pipeline']['blast']['evalue']
#
#    db_name = row.filename + '.db'
#
#    t1 = '{0}.x.{1}.tsv'.format(assembly, db_name)
#    t2 = '{0}.x.{1}.tsv'.format(db_name, assembly)
#
#
#    if row.db_type == 'prot':
#        yield BlastTask('blastx', assembly, db_name, t1,
#                        num_threads=blast_threads, evalue=blast_evalue,
#                        params=blast_params).tasks().next()
#        yield BlastTask('tblastn', row.filename, '{}.db'.format(assembly),
#                        t2, num_threads=blast_threads, evalue=blast_evalue,
#                        params=blast_params).tasks().next()
#    else:
#        yield BlastTask('blastn', assembly, db_name, t1,
#                        num_threads=blast_threads, evalue=blast_evalue,
#                        params=blast_params).tasks().next()
#        yield BlastTask('blastn', row.filename, '{}.db'.format(assembly),
#                        t2, num_threads=blast_threads, evalue=blast_evalue,
#                        params=blast_params).tasks().next()

@create_task_object
def link_file_task(src):
    ''' Soft-link file to the current directory
    '''
    cmd = 'ln -fs {src}'.format(src=src)
    return {'title': title_with_actions,
            'name': 'ln_' + os.path.basename(src),
            'actions': [cmd],
            'targets': [os.path.basename(src)],
            'uptodate': [run_once],
            'clean': [clean_targets]}


@create_task_object
def split_pairs_task(pe_reads_fn):

    left_fn = pe_reads_fn + '.1'
    right_fn = pe_reads_fn + '.2'
    orphan_fn = pe_reads_fn + '.0'
    cmd = 'split-paired-reads.py {pe_reads_fn} -1 {left_fn} -2 {right_fn} -0 {orphan_fn}'.format(**locals())

    name = 'split_paired_' + pe_reads_fn

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [left_fn, right_fn, orphan_fn],
            'file_dep': [pe_reads_fn],
            'clean': [clean_targets] }

@create_task_object
def bowtie2_build_task(input_fn, db_basename, bowtie2_cfg):

    extra_args = bowtie2_cfg['extra_args']
    cmd = 'bowtie2-build {extra_args} {input_fn} {db_basename}'.format(**locals())

    targets = [db_basename+ext for ext in \
                ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']]
    targets.append(db_basename)

    name = 'bowtie2_build_{db_basename}'.format(**locals())
    return {'title': title_with_actions,
            'name': db_basename,
            'actions': [cmd, 'touch {db_basename}'.format(**locals())],
            'targets': targets,
            'file_dep': [input_fn],
            'clean': [clean_targets] }

@create_task_object
def bowtie2_align_task(db_basename, target_fn, bowtie2_cfg, left_fn='', right_fn='', singleton_fn='',
                        read_fmt='-q', samtools_convert=True,
                        encoding='phred33'):

    assert read_fmt in ['-q', '-f']
    assert encoding in ['phred33', 'phred64']
    encoding = '--' + encoding
    if (left_fn or right_fn):
        assert (left_fn and right_fn)
    n_threads = bowtie2_cfg['n_threads']
    extra_args = bowtie2_cfg['extra_args']
    cmd = 'bowtie2 -p {n_threads} {extra_args} {encoding} {read_fmt} -x {db_basename} '.format(**locals())

    file_dep = [db_basename]
    targets = []

    name = 'bowtie2_align' + ''.join('+' + fn if fn else fn for fn in [left_fn, right_fn, singleton_fn, db_basename])

    if left_fn:
        file_dep.extend([left_fn, right_fn])
        left_fn = '-1 ' + left_fn
        right_fn = '-2 ' + right_fn
    if singleton_fn:
        file_dep.append(singleton_fn)
        singleton_fn = '-U ' + singleton_fn
    if samtools_convert:
        targets.append(target_fn + '.bam')
        target_fn = ' | samtools view -Sb - > {target_fn}.bam'.format(**locals())
    else:
        targets.append(target_fn)
        target_fn = '-S ' + target_fn

    cmd = cmd + '{left_fn} {right_fn} {singleton_fn} {target_fn}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': targets,
            'file_dep': file_dep,
            'clean': [clean_targets] }


@create_task_object
def eXpress_task(transcripts_fn, hits_fn, results_folder):

    cmd = 'express -o {results_folder} {transcripts_fn} {hits_fn}'.format(**locals())

    name = 'eXpress_{transcripts_fn}_{results_folder}'.format(**locals())

    targets = [os.path.join(results_folder, 'params.xprs'),
                os.path.join(results_folder, 'results.xprs'),
                results_folder]

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': targets,
            'file_dep': [transcripts_fn, hits_fn],
            'clean': [clean_targets] }

@create_task_object
def samtools_sort_task(bam_fn):

    cmd = 'samtools sort -n {bam_fn} {bam_fn}.sorted'.format(**locals())

    name = 'samtools_sort_{bam_fn}'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [bam_fn + '.sorted.bam'],
            'file_dep': [bam_fn],
            'clean': [clean_targets] }

@create_task_object
def aggregate_express_task(results_files, tpm_target_fn, eff_target_fn, tot_target_fn, label=''):

    import csv

    def recursive_open(file_list, fp_dict, call, *args, **kwargs):
        if not file_list:
            call(fp_dict, *args, **kwargs)
        else:
            fn = file_list.pop()
            with open(fn, 'rb') as fp:
                fp_dict[fn] = fp
                recursive_open(file_list, fp_dict, call, *args, **kwargs)

    def cmd(results_files, tpm_target, eff_target, tot_target):
        with open(tpm_target, 'wb') as tpm_fp, \
             open(eff_target, 'wb') as eff_fp, \
             open(tot_target, 'wb') as tot_fp:

            def agg(fp_dict, tpm_fp, eff_fp, tot_fp):
                readers = [csv.DictReader(fp, delimiter='\t') for fp in fp_dict.values()]
                names = [fn + sep for fn, sep, _ in map(lambda fn: fn.partition('.fq.gz'), fp_dict.keys())]
                print(names)

                tpm_writer = csv.DictWriter(tpm_fp, ['target_id'] + names, delimiter='\t')
                tpm_writer.writeheader()

                eff_writer = csv.DictWriter(eff_fp, ['target_id'] + names, delimiter='\t')
                eff_writer.writeheader()

                tot_writer = csv.DictWriter(tot_fp, ['target_id'] + names, delimiter='\t')
                tot_writer.writeheader()

                for lines in izip(*readers):
                    data = {fn:d['tpm'] for fn, d in zip(names, lines)}
                    data['target_id'] = lines[0]['target_id']
                    tpm_writer.writerow(data)

                    data = {fn:d['eff_counts'] for fn, d in zip(names, lines)}
                    data['target_id'] = lines[0]['target_id']
                    eff_writer.writerow(data)

                    data = {fn:d['tot_counts'] for fn, d in zip(names, lines)}
                    data['target_id'] = lines[0]['target_id']
                    tot_writer.writerow(data)

            fpd = {}
            recursive_open(results_files, fpd, agg, tpm_fp, eff_fp, tot_fp)

    if not label:
        label = 'aggregate_express_task_' + os.path.basename(tpm_target_fn)

    return {'name': label,
            'title': title_with_actions,
            'actions': [(cmd, [results_files, tpm_target_fn, eff_target_fn, tot_target_fn])],
            'targets': [eff_target_fn, tpm_target_fn, tot_target_fn],
            'clean': [clean_targets],
            'file_dep': results_files}

@create_task_object
def trimmomatic_pe_task(left_in, right_in, left_paired_out, left_unpaired_out,
                     right_paired_out, right_unpaired_out, encoding, trim_cfg):

    assert encoding in ['phred33', 'phred64']
    name = 'TrimmomaticPE_' + left_in + '_' + right_in

    params = trim_cfg['params']
    n_threads = trim_cfg['n_threads']
    cmd = 'TrimmomaticPE -{encoding} -threads {n_threads} {left_in} {right_in} {left_paired_out} {left_unpaired_out} '\
          '{right_paired_out} {right_unpaired_out} {params}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'file_dep': [left_in, right_in],
            'targets': [left_paired_out, left_unpaired_out, right_paired_out, right_unpaired_out],
            'clean': [clean_targets]}

@create_task_object
def trimmomatic_se_task(sample_fn, output_fn, encoding, trim_cfg):
    assert encoding in ['phred33', 'phred64']
    name = 'TrimmomaticSE_' + sample_fn

    params = trim_cfg['params']
    n_threads = trim_cfg['n_threads']
    cmd = 'TrimmomaticSE -{encoding} -threads {n_threads} {sample_fn} '\
          '{output_fn} {params}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'file_dep': [sample_fn],
            'targets': [output_fn],
            'clean': [clean_targets]}

@create_task_object
def interleave_task(left_in, right_in, out_fn, label=''):

    if not label:
        label = 'interleave_' + out_fn

    cmd = 'interleave-reads.py {left_in} {right_in} -o {out_fn}'.format(**locals())

    return {'title': title_with_actions,
            'name': label,
            'actions': [cmd],
            'file_dep': [left_in, right_in],
            'targets': [out_fn],
            'clean': [clean_targets]}

@create_task_object
def cat_task(file_list, target_fn):

    cmd = 'cat {files} > {t}'.format(files=' '.join(file_list), t=target_fn)

    return {'title': title_with_actions,
            'name': 'cat_' + target_fn,
            'actions': [cmd],
            'file_dep': file_list,
            'targets': [target_fn],
            'clean': [clean_targets]}

@create_task_object
def group_task(group_name, task_names):
    return {'name': group_name,
            'actions': None,
            'task_dep': task_names}

# python3 BUSCO_v1.1b1/BUSCO_v1.1b1.py -in petMar2.cdna.fa -o petMar2.cdna.busco.test -l vertebrata/ -m trans -c 4
@create_task_object
def busco_task(input_filename, output_dir, busco_db_dir, input_type, busco_cfg):

    name = '_'.join(['busco', input_filename, os.path.basename(busco_db_dir)])

    assert input_type in ['genome', 'OGS', 'trans']
    n_threads = busco_cfg['n_threads']
    busco_path = busco_cfg['path']

    cmd = 'python3 {busco_path} -in {in_fn} -o {out_dir} -l {db_dir} '\
            '-m {in_type} -c {n_threads}'.format(busco_path=busco_path,
            in_fn=input_filename, out_dir=output_dir, db_dir=busco_db_dir,
            in_type=input_type, n_threads=n_threads)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': ['run_' + output_dir,
                        os.path.join('run_' + output_dir, 'short_summary_' + output_dir.rstrip('/'))],
            'file_dep': [input_filename],
            'uptodate': [run_once],
            'clean': [(clean_folder, ['run_' + output_dir])]}

@create_task_object
def cmpress_task(db_fileame):

    cmd = 'cmpress ' + db_filename

    return {'name': 'cmpress:' + os.path.basename(db_filename),
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.i1f', '.i1i', '.i1m', '.i1p']],
            'file_dep': [db_filename],
            'clean': [clean_targets]}

@create_task_object
def cmscan_task(input_filename, output_filename, db_filename, cmscan_cfg, label=''):

    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)

    n_threads = cmscan_cfg['n_threads']
    cmd = 'cmscan --cpu {n_threads} --cut_ga --rfam --nohmmonly --tblout {output_filename}.tbl'\
          ' {db_filename} {input_filename} > {output_filename}.cmscan'.format(**locals())

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename + 'i3p'],
            'targets': [output_filename + '.tbl', output_filename + '.cmscan'],
            'clean': [clean_targets]}

@create_task_object
def hmmpress_task(db_filename, label=''):

    if not label:
        label = 'hmmpress_' + os.path.basename(db_filename)

    cmd = 'hmmpress ' + db_filename

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.h3f', '.h3i', '.h3m', '.h3p']],
            'file_dep': [db_filename],
            'clean': [clean_targets]}

@create_task_object
def hmmscan_task(input_filename, output_filename, db_filename, hmmscan_cfg, label=''):

    if not label:
        label = 'hmmscan_' + os.path.basename(input_filename) + '.x.' + \
                os.path.basename(db_filename)

    n_threads = hmmscan_cfg['n_threads']
    cmd = 'hmmscan --cpu {n_threads} --domtblout {output_filename} \
          {db_filename} {input_filename}'.format(**locals())

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename+'.h3p'],
            'targets': [output_filename],
            'clean': [clean_targets]}

@create_task_object
def transdecoder_orf_task(input_filename, transdecoder_cfg, label=''):

    if not label:
        label = 'TransDecoder.LongOrfs_' + os.path.basename(input_filename)

    min_prot_len = transdecoder_cfg['min_prot_len']
    cmd = 'TransDecoder.LongOrfs -t {input_filename} -m {min_prot_len}'.format(**locals())

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'targets': [input_filename + '.transdecoder_dir/longest_orfs.pep'],
            'clean': [(clean_folder, [input_filename + '.transdecoder_dir'])]}

# TransDecoder.Predict -t lamp10.fasta --retain_pfam_hits lamp10.fasta.pfam-A.out
@create_task_object
def transdecoder_predict_task(input_filename, db_filename, transdecoder_cfg, label=''):

    if not label:
        label = 'TransDecoder.Predict_' + os.path.basename(input_filename)

    orf_cutoff = transdecoder_cfg['orf_cutoff']
    n_threads = transdecoder_cfg['n_threads']

    cmd = 'TransDecoder.Predict -t {input_filename} --retain_pfam_hits {db_filename} \
            --retain_long_orfs {orf_cutoff} --cpu {n_threads}'.format(**locals())

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, input_filename + '.transdecoder_dir/longest_orfs.pep', db_filename],
            'targets': [input_filename + ext for ext in ['.bed', '.cds', '.pep', '.gff3', '.mRNA']],
            'clean': [clean_targets, (clean_folder, [input_filename + '.transdecoder_dir'])]}
