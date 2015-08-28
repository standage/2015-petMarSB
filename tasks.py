#/usr/bin/env python
from __future__ import print_function

from itertools import izip
import json
import os
import pprint
import re
from shutil import rmtree
import sys

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task
import jinja2
import pandas as pd

from peasoup import task_funcs
from peasoup.tasks import BlastTask, BlastFormatTask, CurlTask, GunzipTask, \
                            UniProtQueryTask, TruncateFastaNameTask

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
def diginorm_task(input_files, ksize, htsize, cutoff,
                  n_tables=4, ht_outfn=None, label=None):

    if not label:
        label = 'normalize_by_median_'
        label += '_'.join(input_files[:5])
        if len(input_files) > 5:
            label += '_and_' + str(len(input_files) - 5) + '_more'

    if len(input_files) == 1:
        report_fn = input_files[0] + '.keep.report.txt'
    else:
        report_fn = label + '.report.txt'

    inputs = ' '.join(input_files)
    ht_out_str = ''
    if ht_outfn is not None:
        ht_out_str = '-s ' + ht_outfn
    cmd = 'normalize-by-median.py -k {ksize} -x {htsize} -N {n_tables} '\
          '-C {cutoff} -R {report_fn} {ht_out_str} {inputs}'.format(**locals())

    targets = [fn + '.keep' for fn in input_files]
    targets.append(report_fn)
    if ht_out_str:
        targets.append(ht_outfn)

    return {'title': title_with_actions,
            'name': label,
            'actions': [cmd],
            'file_dep': input_files,
            'targets': targets,
            'clean': [clean_targets]}

@create_task_object
def filter_abund_task(input_files, ct_file, minabund, coverage, label=None):

    if not label:
        label = 'filter_abund_'
        label += '_'.join(input_files[:5])
        if len(input_files) > 5:
            label += '_and_' + str(len(input_files) - 5) + '_more'

    inputs = ' '.join(input_files)
    cmd = 'filter-abund.py -C {minabund} -V -Z {coverage} '\
          '{ct_file} {inputs}'.format(**locals())

    targets = [fn + '.abundfilt' for fn in input_files]

    return {'title': title_with_actions,
            'name': label,
            'actions': [cmd],
            'file_dep': input_files + [ct_file],
            'targets': targets,
            'clean': [clean_targets]}

@create_task_object
def curl_task(row):
    cmd = 'curl -o {target_fn} {url}'.format(target_fn=row.filename+'.gz', url=row.url)
    tsk = task_funcs.general_cmdline_task([], [row.filename+'.gz'], cmd)
    tsk['uptodate'] = [run_once]
    return tsk

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

def blast_task(row, config, assembly):

    blast_threads = config['pipeline']['blast']['threads']
    blast_params = config['pipeline']['blast']['params']
    blast_evalue = config['pipeline']['blast']['evalue']

    db_name = row.filename + '.db'

    t1 = '{0}.x.{1}.tsv'.format(assembly, db_name)
    t2 = '{0}.x.{1}.tsv'.format(db_name, assembly)


    if row.db_type == 'prot':
        yield BlastTask('blastx', assembly, db_name, t1,
                        num_threads=blast_threads, evalue=blast_evalue,
                        params=blast_params).tasks().next()
        yield BlastTask('tblastn', row.filename, '{}.db'.format(assembly),
                        t2, num_threads=blast_threads, evalue=blast_evalue,
                        params=blast_params).tasks().next()
    else:
        yield BlastTask('blastn', assembly, db_name, t1,
                        num_threads=blast_threads, evalue=blast_evalue,
                        params=blast_params).tasks().next()
        yield BlastTask('blastn', row.filename, '{}.db'.format(assembly),
                        t2, num_threads=blast_threads, evalue=blast_evalue,
                        params=blast_params).tasks().next()

@create_task_object
def link_files_task(src):
    cmd = 'ln -fs {src}'.format(src=src)
    return task_funcs.general_cmdline_task([src], [os.path.basename(src)], cmd)

@create_task_object
def bowtie2_build_task(row):
    return task_funcs.bowtie2_build_task(row.filename, row.filename.split('.fasta')[0])

@create_task_object
def split_pairs_task(row):
    return task_funcs.split_pairs_task(row.filename)

@create_task_object
def bowtie2_align_task(row, idx_fn, n_threads):
    target = '{sample}.x.{idx}'.format(sample=row.filename, idx=idx_fn)
    encoding = '--' + row.phred
    if row.paired == False:
        return task_funcs.bowtie2_align_task(idx_fn, target, singleton_fn=row.filename,
                                             num_threads=n_threads, encoding=encoding)
    else:
        return task_funcs.bowtie2_align_task(idx_fn, target, left_fn=row.filename+'.1',
                                                right_fn=row.filename + '.2', singleton_fn=row.filename + '.0',
                                                num_threads=n_threads, encoding=encoding)

@create_task_object
def express_task(hits_fn, transcripts_fn):
    folder = hits_fn.split('.bam')[0]
    return task_funcs.eXpress_task(transcripts_fn, hits_fn, folder)

@create_task_object
def samtools_sort_task(hits_fn):
    return task_funcs.samtools_sort_task(hits_fn)

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
def group_task(group_name, task_names):
    return {'name': group_name,
            'actions': None,
            'task_dep': task_names}

@create_task_object
def download_and_untar_task(url, target_dir, label=''):

    cmd = 'mkdir {target_dir}; curl {url} | tar -xz -C {target_dir}'.format(**locals())

    if not label:
        label = 'dl_untar_' + target_dir.strip('/')

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [target_dir],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [run_once]}

# python3 BUSCO_v1.1b1/BUSCO_v1.1b1.py -in petMar2.cdna.fa -o petMar2.cdna.busco.test -l vertebrata/ -m trans -c 4
@create_task_object
def busco_task(input_filename, output_dir, busco_db_dir, input_type, busco_cfg, label=''):
    
    if not label:
        label = 'busco_' + input_filename

    assert input_type in ['genome', 'OGS', 'trans']
    n_threads = busco_cfg['n_threads']
    busco_path = busco_cfg['path']

    cmd = 'python3 {busco_path} -in {in_fn} -o {out_dir} -l {db_dir} -m {in_type} -c {n_threads}'.format(
            busco_path=busco_path, in_fn=input_filename, out_dir=output_dir, db_dir=busco_db_dir, 
            in_type=input_type, n_threads=n_threads)

    return {'name': label,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output_dir],
            'file_dep': [input_filename],
            'uptodate': [run_once],
            'clean': [(clean_folder, [output_dir])]}

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
            'targets': [input_filename + '.transdecoder_dir'],
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
            'file_dep': [input_filename, input_filename + '.transdecoder_dir', db_filename],
            'targets': [input_filename + ext for ext in ['.bed', '.cds', '.pep', '.gff3', '.mRNA']],
            'clean': [clean_targets, (clean_folder, [input_filename + '.transdecoder_dir'])]}
