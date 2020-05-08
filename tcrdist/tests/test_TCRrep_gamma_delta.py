import pytest
import pandas as pd
import numpy as np
import tcrdist as td
from tcrdist.repertoire import TCRrep
import os

def test_TCRrep_gamma_delta():
    """
    Simple Execution of Workflow on small gamma delta dataset
    """
    # sant.csv was generated from Sant et al. 2019
    # Using
        # 1. sant_et_all_clean_stables.py and 
        # 2. clean.py 
    df = pd.read_csv(os.path.join('tcrdist','test_files_compact','sant.csv'), sep = ",")
    tr = TCRrep(cell_df = df, organism = "human", chains = ['gamma','delta'], db_file='gammadelta_db.tsv')  

    tr.infer_cdrs_from_v_gene(chain = 'gamma', imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain = 'delta',  imgt_aligned=True)

    tr.index_cols = ['clone_id', 'subject', 'v_g_gene', "v_d_gene", 
                    'cdr3_g_aa', 'cdr3_d_aa',
                    'cdr1_g_aa', 'cdr2_g_aa', 'pmhc_g_aa',
                    'cdr1_d_aa', 'cdr2_d_aa', 'pmhc_d_aa']

    tr.deduplicate()
    tr.clone_df

    # AttributeError: 'TCRrep' object has no attribute '_tcrdist_legacy_method_gamma_delta'
    tr._tcrdist_legacy_method_gamma_delta()

