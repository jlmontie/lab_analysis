import json
from bin_sample_composition import bin_reads
import glob
import os
from ncbi_taxonomy_utils import ncbi_taxonomy
import gzip
import pandas as pd
import numpy as np

ncbi = ncbi_taxonomy()

results_dir = '/data/analysis_group1/idbd_rnd/results/200119_NB551702_0107_AH2N3HAFX2/'
batch_file_paths = glob.glob(os.path.join(results_dir, 'batch', '*'))

results = {}
for batch_file in batch_file_paths:
    with open(batch_file) as input:
        batch_obj = json.load(input)
    for lib in batch_obj['libraries']:
        seq_sple = lib['seqSple']
        accession = lib['bioSple']
        for batch_lib in batch_obj['batch']['readsDist']['DNA']:
            if batch_lib['bioSple'] == accession:
                total_reads = batch_lib['postQualityReads']
        summary_paths = lib['diagnosticOutput']
        composition_path_dna = glob.glob(os.path.join(results_dir, 'tax', seq_sple + '*dna.sample_composition.out'))
        composition_path_rna = glob.glob(os.path.join(results_dir, 'tax', seq_sple + '*rna.sample_composition.out'))
        if composition_path_dna:
            sample_comp = bin_reads(composition_path_dna[0], ncbi_class=ncbi, quantification='absolute', ctrl_taxids=[10760])
            for key in sample_comp:
                sample_comp[key] = sample_comp[key] * 10e6 / total_reads
            libtype = 'dna'
        else:
            sample_comp = {
                "Human": None,
                "Bacteria": None,
                "Virus": None,
                "Parasite": None,
                "Fungus": None,
                "Unclassified": None
            }
            libtype = 'rna'
        t7_read_cnt = np.nan
        viral_summary = [path for path in lib['diagnosticOutput'] if 'dna.viral.dxsm.out.summary.gz' in path]
        if viral_summary:
            viral_summary_path = os.path.join(results_dir, viral_summary[0])
            with gzip.open(viral_summary_path, 'rt') as viral_summary_file:
                for line in viral_summary_file:
                    viral_obj = json.loads(line)
                    if viral_obj['reporting_id'] == '26706_10760':
                        t7_read_cnt = viral_obj['read_count']
                        break
        else:
            t7_read_cnt = 0
        bacterial_summary_dna = [path for path in lib['diagnosticOutput'] if 'dna.bacterial.dxsm.out.summary.gz' in path]
        bacterial_summary_rna = [path for path in lib['diagnosticOutput'] if 'rna.bacterial.dxsm.out.summary.gz' in path]
        ecoli_read_cnt = 0
        ecoli_coverage = 0
        if bacterial_summary_dna:
            bacterial_summary_path_dna = os.path.join(results_dir, bacterial_summary_dna[0])
            with gzip.open(bacterial_summary_path_dna, 'rt') as bacterial_summary_file_dna:
                for line in bacterial_summary_file_dna:
                    bacterial_obj = json.loads(line)
                    if bacterial_obj['taxid'] == 562:
                        ecoli_read_cnt = bacterial_obj['read_count']
                        break
        if bacterial_summary_rna:
            bacterial_summary_path_rna = os.path.join(results_dir, bacterial_summary_rna[0])
            with gzip.open(bacterial_summary_path_rna, 'rt') as bacterial_summary_file_rna:
                for line in bacterial_summary_file_rna:
                    bacterial_obj = json.loads(line)
                    if bacterial_obj['taxid'] == 562:
                        for gene in bacterial_obj['gene_info']:
                            if gene['geneid'] == 0:
                                ecoli_coverage = gene['coverage']
                                break
        results.update({
            seq_sple: {
                "Accession": accession,
                "lib_type": libtype,
                "sample_composition": sample_comp,
                "total_reads": total_reads,
                "t7_reads": t7_read_cnt,
                "t7_normalized_reads": t7_read_cnt * 10e6 / total_reads,
                "ecoli_reads": ecoli_read_cnt,
                "ecoli_normalized_reads": ecoli_read_cnt * 10e6 / total_reads,
                "ecoli_coverage": ecoli_coverage,
            }
        })

results_df = pd.DataFrame.from_dict(results, orient='index')
sample_comp_df = pd.DataFrame(results_df['sample_composition'].apply(pd.Series))
concat = pd.concat([results_df, sample_comp_df], axis=1)
concat = concat.drop(columns='sample_composition')
concat.index.name = 'SeqSple'
concat.to_csv('200119-2-1_results.csv')
