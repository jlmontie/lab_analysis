import os
import glob
import json
import gzip
import argparse

import pandas as pd
import numpy as np

from sample_composition_utils import SampleCompParser
from bin_sample_composition import bin_reads
from idbd_bio_utils import NcbiTaxonomy


def parse_sample_comp(composition_path, ctrl_taxa=None, ncbi_tax=None,
    org_taxa=None, genus_taxid=None):
    comp_parser = SampleCompParser(composition_path, ctrl_taxa=ctrl_taxa,
        ncbi_tax=ncbi_tax)
    total_reads = comp_parser.total_reads
    return_dict = {}
    if org_taxa is not None:
        for taxid in org_taxa:
            return_dict.update(
                {ncbi.get_name(taxid): comp_parser.get_taxid_nr(taxid,
                    normalizer=1e7)}
            )
    if ctrl_taxa is not None:
        for taxid in ctrl_taxa:
            return_dict.update(
                {ncbi.get_name(taxid): comp_parser.get_taxid_nr(taxid,
                    normalizer=1e7)}
            )
    if genus_taxid is not None:
        return_dict.update(
            {ncbi.get_name(genus_taxid): comp_parser.get_genus_nr(
                genus_taxid, normalizer=1e7)}
        )

    org_composition = comp_parser.get_org_comp_nr(normalizer=1e7)
    return_dict = {**return_dict, **org_composition}
    return return_dict


def process_rundir(rundir, libtype=None, org_taxa=None, ctrl_taxa=None,
        genus_taxid=None):
    batch_glob = os.path.join(rundir, 'batch', '*')
    batch_paths = glob.glob(batch_glob)
    nr_dict_ls = []
    for batch_path in batch_paths:
        with open(batch_path) as infile:
            batch = json.load(infile)
        for lib in batch['libraries']:
            if libtype == 'rna' or libtype == 'RNA':
                if lib['libType'] == 'DNA':
                    continue
            if libtype == 'dna' or libtype == 'DNA':
                if lib['libType'] == 'RNA':
                    continue
            batch_id = batch['batch']['libBatchId']
            seq_sple = lib['seqSple']
            accession = lib['bioSple']
            first_dx_path = '.'.join(lib['diagnosticOutput'][0].split('.')[:2])
            composition_path = os.path.join(rundir,
                first_dx_path) + '.sample_composition.out'
            if not os.path.isfile(composition_path):
                print(f"No composition file found for {seq_sple}")
                continue
            if ctrl_taxa is None:
                ctrl_taxa = [int(org['reportingId'].split('_')[-1]) for org in
                    lib['internalControls']['organisms']]
            nr_dict = parse_sample_comp(composition_path,
                ctrl_taxa=ctrl_taxa, org_taxa=org_taxa, ncbi_tax=ncbi,
                genus_taxid=genus_taxid)
            nr_dict.update({'batch id': batch_id})
            nr_dict.update({'accession': accession})
            if libtype is None:
                nr_dict.update({'lib type': lib['libType']})
            nr_dict_ls.append(nr_dict)
            df = pd.DataFrame(nr_dict_ls)
            cols = list(df.columns)
            if libtype is not None:
                cols = list(df.columns)
                cols = [cols[-2:]] + cols[:-2]
                df = df[cols]
                df.to_csv(f"{batch_id}_NR_{libtype.upper()}.csv", index=False)
            else:
                cols = list(df.columns)
                cols = cols[-3:] + cols[:-3]
                df = df[cols]
                df.to_csv(f"output/{batch_id}_NR.csv", index=False)


ncbi = NcbiTaxonomy()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=('Get normalized reads from sample composition file in'
            'a run directory')
    )
    parser.add_argument('run_dir', type=str, help='Run directory')
    parser.add_argument(
        '-c', '--ctrl_taxa', type=str,
        help='(optional) comma-separated list of control taxids'
    )
    parser.add_argument(
        '-o', '--org_taxa', type=str,
        help='(optional) comma-separated list of organism taxids'
    )
    parser.add_argument(
        '-l', '--lib_type', type=str,
        help='(optional) only process this lib type, either "dna" or "rna"'
    )
    parser.add_argument(
        '-g', '--genus_taxid', type=int,
        help='(optional) genus taxid for which to get sum of all children'
    )
    args = parser.parse_args()
    if args.ctrl_taxa is not None:
        ctrl_taxa = [int(taxid) for taxid in args.ctrl_taxa.split(',')]
    else:
        ctrl_taxa = None
    if args.org_taxa is not None:
        org_taxa = [int(taxid) for taxid in args.org_taxa.split(',')]
    else:
        org_taxa = None
    os.chdir(os.path.dirname(__file__))
    process_rundir(args.run_dir, libtype=args.lib_type,
        org_taxa=org_taxa, ctrl_taxa=ctrl_taxa, genus_taxid=args.genus_taxid)
