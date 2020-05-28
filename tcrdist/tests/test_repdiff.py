import unittest
import os.path as op
import pandas as pd
import numpy as np
import warnings

import tcrdist as td
from tcrdist.repertoire import TCRrep

class test_stats(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        filename = op.join(td.__path__[0], 'test_files', 'vdjDB_PMID28636592.tsv')
        pd_df = pd.read_csv(filename, sep='\t')
        t_df = td.mappers.vdjdb_to_tcrdist2(pd_df=pd_df)

        t_df = t_df.loc[(t_df.organism == 'HomoSapiens') & (t_df.epitope == 'M1')]

        tr = TCRrep(cell_df=t_df, organism='human')
        tr.infer_cdrs_from_v_gene(chain='alpha')
        tr.infer_cdrs_from_v_gene(chain='beta')
        tr.index_cols =['subject',
                        'cdr3_b_aa']
        tr.deduplicate()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tr.compute_pairwise_all(chain='beta', metric='nw', proceses=1)
        self.pw = tr.cdr3_b_aa_pw

        np.random.seed(110820)
        self.clone_df = tr.clone_df.assign(Visit=np.random.choice(['Pre', 'Post'], size=tr.clone_df.shape[0], p=[0.4, 0.6]),
                                            Stim=np.random.choice(['A', 'B', 'C'], size=tr.clone_df.shape[0], p=[0.4, 0.1, 0.5]))

    def test_chm_NN(self):
        print(self.clone_df.shape)
        print(self.clone_df.head())
        print(self.pw.shape)
        import hierdiff as hd
        print(dir(hd))
        print(hd.__spec__)
        res = hd.neighborhood_tally(df=self.clone_df,
                                  pwmat=self.pw,
                                  x_cols=['Visit', 'Stim'],
                                  count_col='count',
                                  knn_neighbors=50,
                                  knn_radius=None,
                                  subset_ind=None,
                                  cluster_ind=None)

        res = td.stats.neighborhood_diff(self.clone_df, self.pw, x_cols=['Visit', 'Stim'], test_method='chm')
        self.assertTrue(res.shape[0] == self.clone_df.shape[0])

    def test_fishers_NN(self):
        res = td.stats.neighborhood_diff(self.clone_df, self.pw, x_cols=['Visit'], test_method='fishers')
        self.assertTrue(res.shape[0] == self.clone_df.shape[0])

    def test_chi2_NN(self):
        res = td.stats.neighborhood_diff(self.clone_df, self.pw, x_cols=['Visit'], test_method='chi2')
        self.assertTrue(res.shape[0] == self.clone_df.shape[0])

    def test_fishers_HC(self):
        res = td.stats.hcluster_diff(self.clone_df, self.pw, x_cols=['Visit'], test_method='fishers')

    def test_member_summ(self):
        res, Z = td.stats.hcluster_diff(self.clone_df, self.pw, x_cols=['Visit'], test_method='fishers')
        
        summ_df = td.stats.member_summ(res_df=res, clone_df=self.clone_df, count_col='count')
        res = res.join(summ_df)

        summ_df = td.stats.member_summ(res_df=res, clone_df=self.clone_df, count_col='count', addl_cols=['Stim'])
        res = res.join(summ_df)

if __name__ == '__main__':
    unittest.main()
