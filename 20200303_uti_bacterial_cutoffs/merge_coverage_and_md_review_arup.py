import argparse

import pandas as pd
import numpy as np


def main(coverage_path, review_path, outpath):
    coverages = pd.read_csv(coverage_path,
        dtype={'accession': str, 'seq_id': str, 'batch_id': str, 'seq_sple': str,
            'organism': str, 'taxid': str},
        index_col=0)
    dump_orgs = pd.read_excel(review_path,
        sheet_name='organisms', dtype={'Accession': str})
    dump_orgs = dump_orgs.rename(columns={'Accession': 'accession',
        'Organism Name': 'organism', 'Review': 'review', 'Batch ID': 'batch_id',
        'Reporting ID': 'reporting_id'})
    merged = coverages.merge(
        dump_orgs[['batch_id', 'accession', 'organism', 'reporting_id', 'review']],
        on=['batch_id', 'accession', 'organism'], how='left'
    )
    merged['above_cutoff'] = True
    merged.loc[merged['reporting_id'].isna(), 'above_cutoff'] = False
    lab_report = pd.read_csv('Summary_Accuracy_200226.csv')
    lab_report = lab_report.rename(columns={'Accession number': 'accession',
        'ARUP results': 'organism', 'Sample ID': 'sample_id'})
    merged_2 = merged.merge(lab_report[['accession', 'organism', 'sample_id']],
        on=['accession', 'organism'], how='left')
    merged_2.loc[~merged_2['sample_id'].isna(), 'review'] = 'positive'
    merged_2 = merged_2.drop(columns='sample_id')
    print(f"coverages.shape: {coverages.shape}")
    print(f"merged.shape: {merged_2.shape}")
    print(f"merged no na reporting_id: {np.sum(~merged_2['reporting_id'].isna())}")
    print(f"merged no na above_cutoff: {np.sum(merged_2['above_cutoff'])}")
    merged_2.to_csv(outpath, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_coverages')
    parser.add_argument('input_review')
    parser.add_argument('merged_output')
    args = parser.parse_args()
    main(args.input_coverages, args.input_review, args.merged_output)