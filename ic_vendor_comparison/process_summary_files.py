import json
import gzip
import glob
import os
from collections import defaultdict
import pandas as pd


def recursive_dict():
    return defaultdict(recursive_dict)


summary_dir = '/home/jmontgomery/mnt/ic_summary_files/'
all_files = glob.glob(os.path.join(summary_dir, '*'))
path_ls = [path for path in all_files if path.endswith('.gz')]
batch_paths = [path for path in all_files if path.endswith('.json')]
cutoff = 0.7
batch_dict = {}
for batch_path in batch_paths:
    with open(batch_path) as batch_file:
        batch_json = json.load(batch_file)
    for lib in batch_json['libraries']:
        seqsple = lib['seqSple']
        batch_dict.update({seqsple: batch_json})
    batch_id = batch_json['batch']['libBatchId']
    batch_dict.update({batch_id: batch_json})
count_dict = recursive_dict()
for path in path_ls:
    seqsple_path = os.path.basename(path).split('.')[0]
    if seqsple_path.endswith('-d') or seqsple_path.endswith('-r'):
        seqsple_path = seqsple_path[:-2]
    # batch_id = '-'.join(seqsple.split('-')[:2])
    if not seqsple_path in batch_dict:
        print(f"{seqsple_path} not found in batch file")
        continue
    batch = batch_dict[seqsple_path]
    print(f"Analyzing {seqsple_path}")
    total_reads = 0
    for lib in batch['libraries']:
        if seqsple_path in lib['seqSple']:
            libtype = lib['libType']
            total_reads = lib['readInfo']['totalReads']
            accession = lib['bioSple']
    if total_reads == 0:
        print(f"{seqsple_path} not found in batch file")
        continue
    with gzip.open(path, 'rt') as file:
        if '-d-' in path and 'bacterial' in path:
            for line in file:
                obj = json.loads(line)
                if obj['taxid'] == 562:
                    count_dict[accession].update({'genome_coverage': obj['coverage']})
                    count_dict[accession].update({'normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
                    break
            for field in ['genome_coverage', 'normalized_reads']:
                if field not in count_dict[accession]:
                    count_dict[accession].update({field: None})
        elif '-d-' in path and 'viral' in path:
            for line in file:
                obj = json.loads(line)
                if obj['taxid'] == 10665:
                    count_dict[accession].update({'t4_coverage': obj['coverage']})
                    count_dict[accession].update({'t4_normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
                if obj['taxid'] == 10760:
                    count_dict[accession].update({'t7_coverage': obj['coverage']})
                    count_dict[accession].update({'t7_normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
            for field in ['t4_coverage', 't4_normalized_reads', 't7_coverage', 't7_normalized_reads']:
                if field not in count_dict[accession]:
                    count_dict[accession].update({field: None})
        elif '-r-' in path and 'bacterial' in path:
            for line in file:
                obj = json.loads(line)
                if obj['taxid'] == 562:
                    count_dict[accession].update({'normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
                    for gene in obj['gene_info']:
                        if gene['geneid'] == 0:
                            count_dict[accession].update({'16s_coverage': gene['coverage']})
                            break
            for field in ['16s_coverage']:
                if field not in count_dict[accession]:
                    count_dict[accession].update({field: None})
        elif '-r-' in path and 'viral' in path:
            for line in file:
                obj = json.loads(line)
                if obj['taxid'] == 12022:
                    count_dict[accession].update({'ms2_coverage': obj['coverage']})
                    count_dict[accession].update({'ms2_normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
                if obj['taxid'] == 39803:
                    count_dict[accession].update({'qbeta_coverage': obj['coverage']})
                    count_dict[accession].update({'qbeta_normalized_reads': 1e7 * (obj['read_count'] / total_reads)})
            for field in ['ms2_coverage', 'ms2_normalized_reads', 'qbeta_coverage', 'qbeta_normalized_reads']:
                if field not in count_dict[accession]:
                    count_dict[accession].update({field: None})
        for field in ['ms2_coverage', 'ms2_normalized_reads', 'qbeta_coverage',
                      'qbeta_normalized_reads', '16s_coverage',
                      't4_coverage', 't4_normalized_reads', 't7_coverage', 't7_normalized_reads'
                      'genome_coverage', 'normalized_reads']:
            if field not in count_dict[accession]:
                count_dict[accession].update({field: None})

with open('ecoli_results_2.json', 'w') as outfile:
    json.dump(count_dict, outfile)

df = pd.DataFrame.from_dict(count_dict, orient='index').to_csv('ecoli_results_2.csv')