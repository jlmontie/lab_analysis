import json
from bin_sample_composition import bin_reads
import glob
import os
from ncbi_taxonomy_utils import ncbi_taxonomy
ncbi = ncbi_taxonomy()
import gzip
import pandas as pd

results_dir = '/data/analysis_group1/idbd_rnd/results/191104_NB551543_0156_AHTLL7AFXY'
batch_file_paths = glob.glob(os.path.join(results_dir, 'batch', '*'))

results = {}
for batch_file in batch_file_paths:
    with open(batch_file) as input:
        batch_obj = json.load(input)
    for lib in batch_obj['libraries']:
        seq_sple = lib['seqSple']
        accession = lib['bioSple']
        print(accession)
        total_reads = 0
        for batch_lib in batch_obj['batch']['readsDist']['DNA']:
            if batch_lib['bioSple'] == accession:
                total_reads = batch_lib['postQualityReads']
                print(total_reads)
        if total_reads == 0:
            continue
        summary_paths = lib['diagnosticOutput']
        composition_path = glob.glob(os.path.join(results_dir, 'tax', seq_sple + '*sample_composition.out'))
        if composition_path:
            sample_comp = bin_reads(composition_path[0], ncbi_class=ncbi, ctrl_taxids=[10760], quantification='relative')
            # for key in sample_comp:
            #     sample_comp[key] = sample_comp[key] * 10e6 / total_reads
        else:
            sample_comp = {
            "Human": None,
            "Bacteria": None,
            "Virus": None,
            "Parasite": None,
            "Fungus": None,
            "Unclassified": None
        }
        viral_summary = [path for path in lib['diagnosticOutput'] if 'dna.viral.dxsm.out.summary.gz' in path]
        viral_summary_path = os.path.join(results_dir, viral_summary[0])
        with gzip.open(viral_summary_path, 'rt') as viral_summary_file:
            for line in viral_summary_file:
                viral_obj = json.loads(line)
                if viral_obj['reporting_id'] == '26706_10760':
                    t7_read_cnt = viral_obj['read_count']
                    break
        bacterial_summary = [path for path in lib['diagnosticOutput'] if 'dna.bacterial.dxsm.out.summary.gz' in path]
        bacterial_summary_path = os.path.join(results_dir, bacterial_summary[0])
        ecoli_read_cnt = 0
        ecoli_coverage = 0
        kpneumo_read_cnt = 0
        kpneumo_coverage = 0
        with gzip.open(bacterial_summary_path, 'rt') as bacterial_summary_file:
            for line in bacterial_summary_file:
                bacterial_obj = json.loads(line)
                if bacterial_obj['taxid'] == 562:
                    ecoli_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            ecoli_coverage = gene['coverage']
                if bacterial_obj['taxid'] == 573:
                    kpneumo_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            kpneumo_coverage = gene['coverage']
        results.update({
            accession: {
                "sample_composition": sample_comp,
                "total_reads": total_reads,
                "t7_reads": t7_read_cnt,
                "t7_normalized_reads": t7_read_cnt * 10e6 / total_reads,
                "ecoli_reads": ecoli_read_cnt,
                "ecoli_normalized_reads": ecoli_read_cnt * 10e6 / total_reads,
                "ecoli_coverage": ecoli_coverage,
                "kpneumo_reads": kpneumo_read_cnt,
                "kpneumo_normalized_reads": kpneumo_read_cnt * 10e6 / total_reads,
                "kpneumo_coverage": kpneumo_coverage,
            }
        })

results_df = pd.DataFrame.from_dict(results, orient='index')
sample_comp_df = pd.DataFrame(results_df['sample_composition'].apply(pd.Series))
concat = pd.concat([results_df, sample_comp_df], axis=1)
concat = concat.drop(columns='sample_composition')
concat.index.name = 'Accession'
concat.to_csv('191104-1-1_results_relative.csv')
