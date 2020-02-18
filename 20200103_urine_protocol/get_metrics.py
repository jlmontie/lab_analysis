import json
from bin_sample_composition import bin_reads
import glob
import os
from ncbi_taxonomy_utils import ncbi_taxonomy
ncbi = ncbi_taxonomy()
import gzip
import pandas as pd

results_dir = '/data/analysis_group1/idbd_rnd/results/200104_NB551543_0189_AHY2Y5AFXY/'
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
        composition_path = glob.glob(os.path.join(results_dir, 'tax', seq_sple + '*sample_composition.out'))
        if composition_path:
            sample_comp = bin_reads(composition_path[0], ncbi_class=ncbi, quantification='absolute', ctrl_taxids=[10760])
            for key in sample_comp:
                sample_comp[key] = sample_comp[key] * 10e6 / total_reads
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
        bacterial_summary = [path for path in lib['diagnosticOutput'] if 'dna.bacterial.dxsm.out.summary.gz' in path]
        bacterial_summary_path = os.path.join(results_dir, bacterial_summary[0])
        ecoli_read_cnt = 0
        ecoli_coverage = 0
        kpneumo_read_cnt = 0
        kpneumo_coverage = 0
        calbicans_read_cnt = 0
        calbicans_coverage = 0
        saureus_read_cnt = 0
        saureus_coverage = 0
        sagal_read_cnt = 0
        sagal_coverage = 0
        sepiderm_read_cnt = 0
        sepiderm_coverage = 0
        efaecalis_read_cnt = 0
        efaecalis_coverage = 0
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
                if bacterial_obj['taxid'] == 1280:
                    saureus_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            saureus_coverage = gene['coverage']
                if bacterial_obj['taxid'] == 1311:
                    sagal_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            sagal_coverage = gene['coverage']
                if bacterial_obj['taxid'] == 1282:
                    sepiderm_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            sepiderm_coverage = gene['coverage']
                if bacterial_obj['taxid'] == 1351:
                    efaecalis_read_cnt = bacterial_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            efaecalis_coverage = gene['coverage']
                fungal_summary = [path for path in lib['diagnosticOutput'] if 'dna.fungal_parasite.fungal.dxsm.out.summary.gz' in path]
        fungal_summary_path = os.path.join(results_dir, fungal_summary[0])
        with gzip.open(fungal_summary_path, 'rt') as fungal_summary_file:
            for line in fungal_summary_file:
                fungal_obj = json.loads(line)
                if fungal_obj['taxid'] == 5476:
                    calbicans_read_cnt = fungal_obj['read_count']
                    for gene in bacterial_obj['gene_info']:
                        if gene['geneid'] == 0:
                            calbicans_coverage = gene['coverage']
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
                "calbicans_reads": calbicans_read_cnt,
                "calbicans_normalized_reads": calbicans_read_cnt * 10e6 / total_reads,
                "calbicans_coverage": calbicans_coverage,
                "saureus_reads": saureus_read_cnt,
                "saureus_normalized_reads": saureus_read_cnt * 10e6 / total_reads,
                "saureus_coverage": saureus_coverage,
                "sagal_reads": sagal_read_cnt,
                "sagal_normalized_reads": sagal_read_cnt * 10e6 / total_reads,
                "sagal_coverage": sagal_coverage,
                "sepiderm_reads": sepiderm_read_cnt,
                "sepiderm_normalized_reads": sepiderm_read_cnt * 10e6 / total_reads,
                "sepiderm_coverage": sepiderm_coverage,
                "efaecalis_reads": efaecalis_read_cnt,
                "efaecalis_normalized_reads": efaecalis_read_cnt * 10e6 / total_reads,
                "efaecalis_coverage": efaecalis_coverage
            }
        })

results_df = pd.DataFrame.from_dict(results, orient='index')
sample_comp_df = pd.DataFrame(results_df['sample_composition'].apply(pd.Series))
concat = pd.concat([results_df, sample_comp_df], axis=1)
concat = concat.drop(columns='sample_composition')
concat.index.name = 'Accession'
concat.to_csv('200103-1-1_results.csv')
