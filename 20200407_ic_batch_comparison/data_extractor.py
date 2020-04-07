import os
import glob
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from idbd_bio_utils import NcbiTaxonomy

from sample_composition_utils import SampleCompParser

class DataExtractor(object):
    """
    Extract data from summary files and sample composition files for
    batches specified. Batches may be provided in a list in which case
    data from all samples in the batch will be pulled. Alternatively,
    a dictionary may be supplied mapping batch names to accessions for
    which data should be pulled. If both are supplied, data will be
    pulled for all batch samples except for the batches in the
    dictionary.

    Anticipated updates:
    1. Allow specifying specific run directories rather than an entire
    project directory.
    2. Allow for Synergy or other by dynamically changing
    'diagnosticOutput' fieldname.
    """
    def __init__(self, project_dir=None, batch_list=None, accession_dict=None,
            lib_type=None, ctrl_reporting_ids=None):
        """
        args:
            project_dir (str): a directory containing sequencing run
            directories
            batch_list (list): batch ids for the batches to be
            pulled, not case-sensitive
            accession_dict (dict): maps batches to samples that should
            be pulled
            lib_type (str): either 'rna' or 'dna'
            ctrl_reporting_ids (list): reporting ids for the internal
            controls. This overrides the control reporting_id in the
            batch file.
        """
        assert all((project_dir is not None, lib_type is not None)), \
                (f"project_dir and lib_type must be supplied")
        assert any((batch_list is not None, accession_dict is not None)), \
            (f"either batch_list or accession_dict must be supplied")
        self.project_dir = project_dir
        if batch_list is not None:
            self.batch_list = [batch.lower() for batch in batch_list]
        else:
            self.batch_list = [batch.lower() for batch in
                list(accession_dict.keys())]
        if accession_dict is not None:
            accessions = []
            for item in accession_dict.values():
                accessions.extend(item)
            self.accessions = list(set([accession.lower() for accession in
                accessions]))
        else:
            self.accessions = []
        self.accession_dict = accession_dict
        self.lib_type = lib_type
        self.ctrl_reporting_ids = ctrl_reporting_ids
        self.rundirs = glob.glob(os.path.join(self.project_dir, '*'))
        assert self.rundirs, f"No run directories found in {self.project_dir}"
        self.batches_found = set()
        self.ncbi = NcbiTaxonomy()
        self.df = pd.DataFrame()

    def collect_data(self):
        """
        Pull data for the batches
        returns: pandas DataFrame with data for all batches
        """
        rundir_dict_ls = []
        for rundir in tqdm(self.rundirs):
            batch_files = glob.glob(os.path.join(rundir, 'batch', '*'))
            if not batch_files:
                continue
            matching_batches = []
            for batch in batch_files:
                batch_basename = os.path.basename(batch).split('.')[0]
                if batch_basename.lower() in self.batch_list:
                    matching_batches.append(batch)
                    self.batches_found.add(batch_basename.lower())
            if not matching_batches:
                continue
            batch_dict_ls = self._iterate_batches(matching_batches)
            rundir_dict_ls.extend(batch_dict_ls)
        self.df = self._create_dataframe(rundir_dict_ls)
        self._print_unprocessed_batches()

    def save_data(self, outpath):
        self.df.to_csv(outpath, index=False)

    def _iterate_batches(self, matching_batches):
        batch_dict_ls = []
        for batch in matching_batches:
            lib_dict_ls = self._parse_batch_file(batch)
            batch_dict_ls.extend(lib_dict_ls)
        return batch_dict_ls

    def _parse_batch_file(self, batch_path):
        with open(batch_path) as batch_file:
            batch = json.load(batch_file)
            batch_id = batch['batch']['libBatchId']
        lib_dict_ls = []
        for lib in batch['libraries']:
            accession = lib['bioSple']
            # skip samples not in accession_dict if it exists
            if self.accession_dict is not None:
                if batch_id.lower() in [batch.lower() for batch in
                                        self.accession_dict.keys()]:
                    if accession.lower() not in self.accessions:
                        continue
            lib_type = lib['libType']
            # skip samples if not matching lib_type if supplied
            if self.lib_type is not None:
                if not lib_type.lower() == self.lib_type.lower():
                    continue
                # partial paths for viral summary and composition files
                vir_partial = self.lib_type + '.viral.dxsm.out.summary.gz'
                comp_partial = '*' + self.lib_type + '.sample_composition.out'
            else:
                vir_partial = '.viral.dxsm.out.summary.gz'
                comp_partial = '*.sample_composition.out'
            seq_sple = lib['seqSple']
            run_date = batch['analysis']['timeCompleted'][:10]
            sample_name = lib['spleName']
            total_reads = lib['qualityFilterInfo']['readsOut']
            tax_paths = lib['diagnosticOutput']
            run_dir = Path(batch_path).parent
            ### Control data ###
            vir_paths = [os.path.join(run_dir, path) for path in tax_paths if
                vir_partial in path]
            if self.ctrl_reporting_ids:
                ctrl_reporting_ids = self.ctrl_reporting_ids
            else:
                ctrl_reporting_ids = [item['reportingId'] for item in
                                      lib['internalControls']['organisms']]
            norm_ctrl_read_cnt_summary = self._get_ctrl_cnts(vir_paths,
                ctrl_reporting_ids, total_reads)
            ### Composition data ###
            composition_path_ls = glob.glob(os.path.join(rundir, 'tax',
                seq_sple + comp_partial))
            norm_ctrl_read_cnt_comp, org_composition = \
                self._get_composition_data(composition_path_ls,
                    ctrl_reporting_ids)
            lib_dict_ls.append(self._build_lib_dict())
        return lib_dict_ls

    def _get_ctrl_cnts(self, vir_paths, ctrl_reporting_ids, total_reads):
        ctrl_read_cnt = 0
        if vir_paths:
            vir_path = vir_paths[0]
            with gzip.open(vir_path) as vir_file:
                for line in vir_file:
                    obj = json.loads(line.strip())
                    if obj['reporting_id'] in ctrl_reporting_ids:
                        ctrl_read_cnt += int(obj['read_count'])
        norm_ctrl_read_cnt = 10e6 * ctrl_read_cnt / total_reads
        return norm_ctrl_read_cnt

    def _get_composition_data(self, composition_path_ls, ctrl_reporting_ids):
        taxids = [item.split('_')[1] for item in ctrl_reporting_ids]
        if not composition_path_ls:
            norm_ctrl_read_cnt = 0
            org_composition = {
                "Human": 0,
                "Bacteria": 0,
                "Virus": 0,
                "Parasite": 0,
                "Fungus": 0,
                "Unclassified": 0
            }
        else:
            composition_path = composition_path_ls[0]
            composition_file = os.path.basename(composition_path)
            composition_parser = SampleCompParser(composition_path, self.ncbi,
                taxids)
            norm_ctrl_read_cnt = np.sum([
                (taxid, composition_parser.get_taxid_reads(taxid))
                for taxid in taxids
            ])
            org_composition = composition_parser.get_org_comp_nr()
        return norm_ctrl_read_cnt, org_composition

    def _build_lib_dict(self, accession, seq_sple, batch_id, run_date,
            sample_name, lib_type, total_reads, norm_ctrl_read_cnt_summary,
            norm_ctrl_read_cnt_comp, org_composition):
        lib_dict = {
            'Accession': accession,
            'Seq Sple': seq_sple,
            'Batch ID': batch_id,
            'Run Date': run_date,
            'Sample Name': sample_name,
            'Library Type': lib_type,
            'Total Reads': total_reads,
            'Total Ctrl NR': norm_ctrl_read_cnt_summary,
            'Log10 Total Ctrl NR': np.log10(norm_ctrl_read_cnt_summary),
            'Total Ctrl NR (Composition File)': norm_ctrl_read_cnt_comp,
            'Log10 Total Ctrl NR (Composition File)': \
                np.log10(norm_ctrl_read_cnt_comp),
            "Human": org_composition['human'],
            "Bacteria": org_composition['bacteria'],
            "Virus": org_composition['virus'],
            "Parasite": org_composition['parasite'],
            "Fungus": org_composition['fungus'],
            "Unclassified": org_composition['unclassified']
        }
        return lib_dict

    def _create_dataframe(self, rundir_dict_ls):
        df = pd.DataFrame(rundir_dict_ls)
        return df

    def _print_unprocessed_batches(self):
        batch_set = self.batch_list
        batches_not_found = batch_set.difference(self.batches_found)
        print(f"Batches not found: {list(batches_not_found)}")
