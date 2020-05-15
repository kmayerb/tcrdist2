import unittest
from tcrdist import pairwise
from tcrdist import repertoire
from tcrdist.repertoire import TCRrep, load_hdf5
import parasail
import pandas as pd
import numpy as np

import os
import tempfile

example_df = pd.DataFrame({'count': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1, 18: 1, 19: 1}, 'j_b_gene': {0: 'TRBJ1-2*01', 1: 'TRBJ1-2*01', 2: 'TRBJ1-2*01', 3: 'TRBJ1-2*01', 4: 'TRBJ1-2*01', 5: 'TRBJ1-2*01', 6: 'TRBJ1-2*01', 7: 'TRBJ1-5*01', 8: 'TRBJ1-2*01', 9: 'TRBJ1-2*01', 10: 'TRBJ1-2*01', 11: 'TRBJ1-2*01', 12: 'TRBJ1-2*01', 13: 'TRBJ1-2*01', 14: 'TRBJ2-3*01', 15: 'TRBJ1-5*01', 16: 'TRBJ2-7*01', 17: 'TRBJ1-1*01', 18: 'TRBJ2-7*01', 19: 'TRBJ2-7*01'}, 'j_a_gene': {0: 'TRAJ42*01', 1: 'TRAJ42*01', 2: 'TRAJ42*01', 3: 'TRAJ50*01', 4: 'TRAJ42*01', 5: 'TRAJ42*01', 6: 'TRAJ42*01', 7: 'TRAJ20*01', 8: 'TRAJ42*01', 9: 'TRAJ42*01', 10: 'TRAJ42*01', 11: 'TRAJ42*01', 12: 'TRAJ42*01', 13: 'TRAJ42*01', 14: 'TRAJ49*01', 15: 'TRAJ33*01', 16: 'TRAJ42*01', 17: 'TRAJ49*01', 18: 'TRAJ31*01', 19: 'TRAJ37*02'}, 'cdr3_a_aa': {0: 'CAGQASQGNLIF', 1: 'CAGQASQGNLIF', 2: 'CAGQASQGNLIF', 3: 'CAGPRETSYDKVIF', 4: 'CAGQASQGNLIF', 5: 'CAGQASQGNLIF', 6: 'CAGQASQGNLIF', 7: 'CAETRSRDYKLSF', 8: 'CAGQASQGNLIF', 9: 'CAGQASQGNLIF', 10: 'CAGQASQGNLIF', 11: 'CAGQASQGNLIF', 12: 'CAGQASQGNLIF', 13: 'CAGQASQGNLIF', 14: 'CAVADTGNQFYF', 15: 'CLVGSMDSNYQLIW', 16: 'CAVPKGSQGNLIF', 17: 'CAVSDSGTGNQFYF', 18: 'CAGPFGRLMF', 19: 'CAGPDGSSNTGKLIF'}, 'epitope': {0: 'pp65', 1: 'pp65', 2: 'pp65', 3: 'pp65', 4: 'pp65', 5: 'pp65', 6: 'pp65', 7: 'pp65', 8: 'pp65', 9: 'pp65', 10: 'pp65', 11: 'pp65', 12: 'pp65', 13: 'pp65', 14: 'pp65', 15: 'M1', 16: 'M1', 17: 'M1', 18: 'M1', 19: 'M1'}, 'cdr3_b_aa': {0: 'CASSIQALLTF', 1: 'CASSIQALLTF', 2: 'CASSIQALLTF', 3: 'CASSSAYYGYTF', 4: 'CASSIQALLTF', 5: 'CASSIQALLTF', 6: 'CASSIQALLTF', 7: 'CASSQEEGPGNQPQHF', 8: 'CASSIQALLTF', 9: 'CASSIQALLTF', 10: 'CASSIQALLTF', 11: 'CASSIQALLTF', 12: 'CASSIQALLTF', 13: 'CASSIQALLTF', 14: 'CATAITSTQYF', 15: 'CASSSQSNQPQHF', 16: 'CASSIRSSYEQYF', 17: 'CASSQMTGLNTEAFF', 18: 'CASSLFPGFGEQYF', 19: 'CASSLIFPSGEQYF'}, 'v_b_gene': {0: 'TRBV12-3*01', 1: 'TRBV12-3*01', 2: 'TRBV12-3*01', 3: 'TRBV12-3*01', 4: 'TRBV12-3*01', 5: 'TRBV12-3*01', 6: 'TRBV12-3*01', 7: 'TRBV4-1*01', 8: 'TRBV12-3*01', 9: 'TRBV12-3*01', 10: 'TRBV12-3*01', 11: 'TRBV12-3*01', 12: 'TRBV12-3*01', 13: 'TRBV12-3*01', 14: 'TRBV12-3*01', 15: 'TRBV25-1*01', 16: 'TRBV19*01', 17: 'TRBV28*01', 18: 'TRBV27*01', 19: 'TRBV27*01'}, 'id': {0: 'human_tcr0001', 1: 'human_tcr0002', 2: 'human_tcr0003', 3: 'human_tcr0004', 4: 'human_tcr0005', 5: 'human_tcr0006', 6: 'human_tcr0007', 7: 'human_tcr0008', 8: 'human_tcr0009', 9: 'human_tcr0010', 10: 'human_tcr0011', 11: 'human_tcr0012', 12: 'human_tcr0013', 13: 'human_tcr0014', 14: 'human_tcr0015', 15: 'human_tcr0016', 16: 'human_tcr0017', 17: 'human_tcr0018', 18: 'human_tcr0019', 19: 'human_tcr0020'}, 'v_a_gene': {0: 'TRAV35*01', 1: 'TRAV35*01', 2: 'TRAV35*01', 3: 'TRAV35*02', 4: 'TRAV35*01', 5: 'TRAV35*01', 6: 'TRAV35*01', 7: 'TRAV5*01', 8: 'TRAV35*01', 9: 'TRAV35*01', 10: 'TRAV35*01', 11: 'TRAV35*01', 12: 'TRAV35*01', 13: 'TRAV35*01', 14: 'TRAV22*01', 15: 'TRAV4*01', 16: 'TRAV8-3*02', 17: 'TRAV8-6*02', 18: 'TRAV27*01', 19: 'TRAV35*02'}, 'subject': {0: 'human_subject0010', 1: 'human_subject0010', 2: 'human_subject0010', 3: 'human_subject0010', 4: 'human_subject0010', 5: 'human_subject0010', 6: 'human_subject0010', 7: 'human_subject0010', 8: 'human_subject0010', 9: 'human_subject0010', 10: 'human_subject0010', 11: 'human_subject0010', 12: 'human_subject0010', 13: 'human_subject0010', 14: 'human_subject0010', 15: 'human_subject0007', 16: 'human_subject0015', 17: 'human_subject0007', 18: 'human_subject0007', 19: 'human_subject0007'}})

class test_repertoire(unittest.TestCase):

    def test_TCRrep___init(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        self.assertTrue(testrep)

    def test_TCRrep___init_creates_cell_df_as_pdDF(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        self.assertTrue(isinstance(testrep.cell_df, pd.DataFrame))

    def test_TCRrep___init_creates_cdr3_a_aa_smat(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        """Changed this test now that matrix for nw_metric is required
        to be the name of a parasail matrix (e.g. "blosum62")"""
        self.assertTrue(isinstance(testrep.cdr3_a_aa_smat, str))

    def test_TCRrep___init_creates_cdr3_b_aa_smat(self):
        """Changed this test now that matrix for nw_metric is required
        to be the name of a parasail matrix (e.g. "blosum62")"""
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        self.assertTrue(isinstance(testrep.cdr3_b_aa_smat, str))

    def test_TCRrep___deduplicate_creates_clone_df_as_pdDF(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        testrep.deduplicate()
        self.assertTrue(isinstance(testrep.clone_df, pd.DataFrame))


    def test_repertoire_infer_cdrs_from_v_gene(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        testrep.infer_cdrs_from_v_gene(chain = "alpha")
        testrep.infer_cdrs_from_v_gene(chain = "beta")
        testrep.index_cols =['epitope', 'subject', 'cdr3_a_aa', 'cdr1_a_aa',
                             'cdr2_a_aa', 'pmhc_a_aa', 'cdr3_b_aa', 'cdr1_b_aa',
                             'cdr2_b_aa', 'pmhc_b_aa']
        testrep.deduplicate()
        self.assertTrue(np.all([isinstance(testrep.clone_df['cdr1_a_aa'], pd.Series),
                                isinstance(testrep.clone_df['cdr1_b_aa'], pd.Series),
                                isinstance(testrep.clone_df['cdr2_a_aa'], pd.Series),
                                isinstance(testrep.clone_df['cdr2_b_aa'], pd.Series)]))

    def test_repertoire_infer_olga_pgen_from_cdr3_and_vj_genes(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        testrep.infer_cdrs_from_v_gene(chain = "alpha")
        testrep.infer_cdrs_from_v_gene(chain = "beta")
        testrep.index_cols =['epitope', 'subject', 'v_a_gene', 'j_a_gene',
                             'v_b_gene', 'j_b_gene', 'cdr3_a_aa', 'cdr1_a_aa',
                             'cdr2_a_aa', 'pmhc_a_aa', 'cdr3_b_aa', 'cdr1_b_aa',
                             'cdr2_b_aa', 'pmhc_b_aa']
        testrep.deduplicate()
        testrep.infer_olga_aa_cdr3_pgens(chain = "alpha")
        testrep.infer_olga_aa_cdr3_pgens(chain = "beta")
        #print(testrep.clone_df.cdr3_a_aa_pgen)
        self.assertTrue(isinstance(testrep.clone_df.cdr3_a_aa_pgen, pd.Series))
        self.assertTrue(np.all([ isinstance(testrep.clone_df.cdr3_a_aa_pgen, pd.Series),
                                isinstance(testrep.clone_df.cdr3_b_aa_pgen, pd.Series)]))

    def test_repertoire_infer_olga_pgen_from_cdr3_only_genes(self):
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"])
        testrep.infer_cdrs_from_v_gene(chain = "alpha")
        testrep.infer_cdrs_from_v_gene(chain = "beta")
        testrep.index_cols =['epitope', 'subject', 'v_a_gene', 'j_a_gene',
                             'v_b_gene', 'j_b_gene', 'cdr3_a_aa', 'cdr1_a_aa',
                             'cdr2_a_aa', 'pmhc_a_aa', 'cdr3_b_aa', 'cdr1_b_aa',
                             'cdr2_b_aa', 'pmhc_b_aa']
        testrep.deduplicate()
        testrep.infer_olga_aa_cdr3_pgens(chain = "alpha",
                                         cdr3_only = True)
        testrep.infer_olga_aa_cdr3_pgens(chain = "beta",
                                         cdr3_only = True)
        #print(testrep.clone_df.cdr3_a_aa_pgen)
        self.assertTrue(np.all([ isinstance(testrep.clone_df.cdr3_a_aa_pgen, pd.Series),
                                 isinstance(testrep.clone_df.cdr3_b_aa_pgen, pd.Series)]))

    def test_repertoire_full_use_case(self):
        """
        This is not a unit test persay! This is a test of a use_case of an instance of
        the TCRrep() class used for pairwise sequence comparison

        """

        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"]) # (1)
        testrep.index_cols.append("epitope")                         # (2)
        testrep.index_cols.append("subject")
        testrep.deduplicate()                                    # (3)
        testrep.cdr3_a_aa_smat = 'blosum62'               # (4)
        testrep.cdr3_b_aa_smat = 'blosum62'
        testrep.compute_pairwise(chain = "alpha")                # (5)
        testrep.compute_pairwise(chain = "beta")                 # (6)
        tcrdist = testrep.cdr3_a_aa_pw + testrep.cdr3_b_aa_pw    # (7)

        expected_tcrdist = np.array([[   0.,  222.,  210.,  223.,  231.,  239.,  219.,  231.,  175.],
               [ 222.,    0.,  116.,  175.,  173.,  185.,  131.,  209.,  205.],
               [ 210.,  116.,    0.,  175.,  169.,  183.,  145.,  201.,  221.],
               [ 223.,  175.,  175.,    0.,  154.,  200.,  162.,  234.,  202.],
               [ 231.,  173.,  169.,  154.,    0.,  152.,  120.,  182.,  192.],
               [ 239.,  185.,  183.,  200.,  152.,    0.,  146.,  112.,  192.],
               [ 219.,  131.,  145.,  162.,  120.,  146.,    0.,  178.,  172.],
               [ 231.,  209.,  201.,  234.,  182.,  112.,  178.,    0.,  220.],
               [ 175.,  205.,  221.,  202.,  192.,  192.,  172.,  220.,    0.]])

        self.assertTrue((tcrdist == expected_tcrdist).all())


    def test_repertoire_full_use_case_hamming(self):
        """
        This is not a unit test persay! This is a test of a use_case of an instance of
        the TCRrep() class used for pairwise sequence comparison
        """
        testrep = TCRrep(cell_df = example_df, chains = ["alpha", "beta"]) # (1)
        testrep.index_cols.append("epitope")                         # (2)
        testrep.index_cols.append("subject")
        testrep.deduplicate()                                    # (3)
        testrep.cdr3_a_aa_smat = 'blosum62'               # (4)
        testrep.cdr3_b_aa_smat = 'blosum62'
        testrep.compute_pairwise(chain = "alpha", metric = "hamming")                # (5)
        testrep.compute_pairwise(chain = "beta", metric = "hamming")                 # (6)
        tcrdist = testrep.cdr3_a_aa_pw + testrep.cdr3_b_aa_pw    # (7)

        expected_tcrdist = np.array([[  0.,  18.,  17.,  18.,  19.,  22.,  19.,  18.,  16.],
           [ 18.,   0.,  11.,  15.,  15.,  17.,  10.,  18.,  18.],
           [ 17.,  11.,   0.,  18.,  15.,  17.,  13.,  18.,  20.],
           [ 18.,  15.,  18.,   0.,  14.,  19.,  14.,  20.,  18.],
           [ 19.,  15.,  15.,  14.,   0.,  14.,  11.,  17.,  16.],
           [ 22.,  17.,  17.,  19.,  14.,   0.,  14.,  13.,  18.],
           [ 19.,  10.,  13.,  14.,  11.,  14.,   0.,  17.,  15.],
           [ 18.,  18.,  18.,  20.,  17.,  13.,  17.,   0.,  19.],
           [ 16.,  18.,  20.,  18.,  16.,  18.,  15.,  19.,   0.]])

        self.assertTrue((tcrdist == expected_tcrdist).all())

    def test_repertoire____use_case_hamming_paired_tcrdist(self):
        tr = TCRrep(cell_df = example_df.copy(), organism = "human", chains= ["alpha", "beta"])
        tr.infer_cdrs_from_v_gene(chain = "alpha")
        tr.infer_cdrs_from_v_gene(chain = "beta")
        tr.index_cols =['cdr3_a_aa', 'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa', 'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']
        tr.deduplicate()
        #tr.clone_df
        tr.compute_pairwise_all(chain = "alpha", metric = "hamming")
        tr.compute_pairwise_all(chain = "beta", metric = "hamming")
        r = tr.compute_paired_tcrdist(chains = ['alpha', 'beta'])

        expected = {'paired_tcrdist': np.array([[  0.,  50.,  49.,  48.,  49.,  49.,  51.,  44.,  46.],
                [ 50.,   0.,  21.,  29.,  29.,  47.,  40.,  45.,  44.],
                [ 49.,  21.,   0.,  42.,  39.,  47.,  42.,  44.,  50.],
                [ 48.,  29.,  42.,   0.,  14.,  35.,  48.,  52.,  46.],
                [ 49.,  29.,  39.,  14.,   0.,  30.,  45.,  49.,  44.],
                [ 49.,  47.,  47.,  35.,  30.,   0.,  46.,  43.,  48.],
                [ 51.,  40.,  42.,  48.,  45.,  46.,   0.,  36.,  41.],
                [ 44.,  45.,  44.,  52.,  49.,  43.,  36.,   0.,  47.],
                [ 46.,  44.,  50.,  46.,  44.,  48.,  41.,  47.,   0.]]),
         'paired_tcrdist_weights': {'cdr1_a_aa_pw': 1,
          'cdr1_b_aa_pw': 1,
          'cdr2_a_aa_pw': 1,
          'cdr2_b_aa_pw': 1,
          'cdr3_a_aa_pw': 1,
          'cdr3_b_aa_pw': 1,
          'pmhc_a_aa_pw': 1,
          'pmhc_b_aa_pw': 1}}
        #print(r['paired_tcrdist'][1, 2])
        #print(expected['paired_tcrdist'][1, 2])
        #print(tr.clone_df.iloc[1])
        #print(tr.clone_df.iloc[2])
        #print(r['paired_tcrdist'] == expected['paired_tcrdist'])

        self.assertTrue((r['paired_tcrdist'] == expected['paired_tcrdist']).all())

    def test_save_as_hdf5(self):
        testrep = TCRrep(cell_df=example_df, chains=["alpha", "beta"]) # (1)
        testrep.index_cols.append("epitope")                         # (2)
        testrep.index_cols.append("subject")
        testrep.deduplicate()                                    # (3)
        testrep.cdr3_a_aa_smat = 'blosum62'               # (4)
        testrep.cdr3_b_aa_smat = 'blosum62'
        testrep.compute_pairwise(chain="alpha")                # (5)
        testrep.compute_pairwise(chain="beta")                 # (6)
        tcrdist = testrep.cdr3_a_aa_pw + testrep.cdr3_b_aa_pw    # (7)

        tmpf = tempfile.NamedTemporaryFile(mode='w+b', suffix='.h5', prefix='tcrrep_test', delete=False)
        tmpf.close()
        testrep.save_as_hdf5(tmpf.name)
        self.assertTrue(os.path.exists(tmpf.name))

        incoming = load_hdf5(tmpf.name)
        os.unlink(tmpf.name)
        incoming_tcrdist = incoming.cdr3_a_aa_pw + incoming.cdr3_b_aa_pw

        self.assertTrue((tcrdist == incoming_tcrdist).all())
        self.assertTrue(testrep.chains == incoming.chains)


if __name__ == '__main__':
    unittest.main()
