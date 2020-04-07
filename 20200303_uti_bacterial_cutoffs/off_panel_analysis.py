import json
from collections import defaultdict

import pandas as pd
import numpy as np

on_panel = pd.read_csv('org_taxids_uti_all.txt', sep='\t')
on_panel_taxids = on_panel['taxid'].unique().tolist()

# arup
arup = pd.read_csv('arup_coverages_with_review.csv')
arup['on_panel'] = arup['taxid'].isin(on_panel_taxids)
arup['source'] = 'ARUP'

arup_fungpar = pd.read_csv('arup_coverages_fungpar.csv')
arup_fungpar['on_panel'] = arup_fungpar['taxid'].isin(on_panel_taxids)
arup_fungpar['source'] = 'ARUP'

arup_vir = pd.read_csv('arup_coverages_vir.csv')
arup_vir['on_panel'] = arup_vir['taxid'].isin(on_panel_taxids)
arup_vir['source'] = 'ARUP'

# synergy
synergy = pd.read_csv('synergy_coverages_with_review.csv')
synergy['on_panel'] = synergy['taxid'].isin(on_panel_taxids)
synergy['source'] = 'Synergy'

synergy_fungpar = pd.read_csv('arup_coverages_fungpar.csv')
synergy_fungpar['on_panel'] = synergy_fungpar['taxid'].isin(on_panel_taxids)
synergy_fungpar['source'] = 'Synergy'

synergy_vir = pd.read_csv('arup_coverages_vir.csv')
synergy_vir['on_panel'] = synergy_vir['taxid'].isin(on_panel_taxids)
synergy_vir['source'] = 'Synergy'

merged = pd.concat([arup, arup_fungpar, arup_vir, synergy, synergy_fungpar,
    synergy_vir])

# add class type from reporting names file
reporting_names = pd.read_csv('explify_uti_reporting_name_info_table.txt',
    sep='\t')
reporting_names = reporting_names.rename(
    columns={'reporting_name': 'organism'}
)

merged_class = merged.merge(
    reporting_names[['organism', 'class_type', 'subclass']],
    on='organism', how='left'
)

# print(merged_class.info())

# print(merged_class['class_type'].value_counts())
# print(merged_class['subclass'].value_counts())
# print(synergy['on_panel'].value_counts())


# print(arup.columns)
# print(synergy.columns)

merged_filtered = merged_class[
    (
        (merged_class['class_type'].isin(['bacterial', 'fungal', 'parasite']))
        & (merged_class['rna_coverage'] >= 0.98)
    )
    | (
        (merged_class['class_type'] == 'viral')
        & (
            (merged_class['rna_coverage'] >= 0.5)
            | (merged_class['dna_coverage'] >= 0.5)
        )
    )
]

print(merged_filtered.columns)
print(merged_filtered['class_type'].value_counts())
merged_filtered = merged_filtered[merged_filtered['on_panel'] == False]
merged_filtered = merged_filtered[
    ~(merged_filtered['organism'] == 'Enterobacteria phage T7')
]

# medical relevance
mr = pd.read_csv(
    'explify_uti_mr_pubmed_disease_category_intersection_counts.txt', sep='\t'
)
mr = mr.rename(columns={'reporting_name': 'organism'})

merged_mr = merged_filtered.merge(
    mr[['organism', 'uti', 'total_mr_pub_count']], on='organism', how='left'
)
print(merged_mr.info())
merged_mr = merged_mr.drop(
    columns=['subclass', 'on_panel', 'seq_id', 'above_cutoff', 'review',
        'reporting_id']
)
merged_mr.to_csv('uti_off_panel_high_coverage.csv', index=False)
