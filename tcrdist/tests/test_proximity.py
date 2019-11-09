import unittest
import os.path as op
import inspect
import pandas as pd
import numpy as np

import tcrdist as td
from tcrdist.repertoire import TCRrep
from tcrdist.proximity import TCRproximity

class test_proximity(unittest.TestCase):
    def test_proximity_with_data(self):
        filename = op.join(td.__path__[0], 'datasets', 'vdjDB_PMID28636592.tsv')
        pd_df = pd.read_csv(filename, sep='\t')
        t_df = td.mappers.vdjdb_to_tcrdist2(pd_df=pd_df)

        t_df = t_df.loc[(t_df.organism == 'HomoSapiens') & ((t_df.epitope == 'M1') | (t_df.epitope == 'BMLF'))]

        tr = TCRrep(cell_df=t_df, organism='human')
        tr.infer_cdrs_from_v_gene(chain='alpha')
        tr.infer_cdrs_from_v_gene(chain='beta')
        tr.index_cols =['subject', 'cdr3_b_aa', 'epitope']
        tr.deduplicate()
        tr.compute_pairwise_all(chain='beta', metric='nw', proceses=1)

        prox_M1 = TCRproximity(tr, 'M1', ['BMLF'], 10, chain='beta', cdrs='cdr3')
        prox_M1.plot_NN_score_distribution()
        prox_M1.plot_ROC()

if __name__ == '__main__':
    unittest.main()