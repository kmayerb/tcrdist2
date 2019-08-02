"""
This is an important test because it confirms that tcrdist2 reproduces
the metric results used in the paper described as the 'tcr-dist'.

There is one point that requires additional clarification.

The original tcr-dist pipeline acknowledged uncertainty in v-gene alpha
v-gene beta inferences. In some cases multipe genes could have likely produced
the sequenced read. When multiple alleles were plausible, the distance metric
between two TCRs was based on a comparison of all possible pairwise distances.

For example suppose both TRAV1*01 and TRAV1*02 were plausible v-gene candidates for a
TCR1, and TRAV2*02 was a the only plausible v-gene for TCR2. Then, the distance
D(TCR1, TCR2) was the minimum { D(TRAV1*01,TRAV2*02 ), D(TRAV1*02,TRAV2*02)}

This was made possible by the fact that there are only so many possible
distances between v-gene combinations given the fixed alignment assumption and a
known v-gene database.

By contrast, tcrdist2 works differently. One benefit is that
it can actually re-align the cdr1, cdr2, pmhc regions inferred
from the v-gene calls. This is done when using the
method 'hamming' or 'nw'.  However, the approach for multiprocessing
these comparisons is not well suited to making
a variable number of pairwise sequences comparisons for each clone.

While tcrdist2 incorporates legacy metrics 'tcrdist_cdr1, tcrdist_cdr3' and
also allows the use of the IMGT alignments for the input comparisons,

    infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned = True) # <---

it does not currently calculate the distance metric for multiple possible
v-gene alleles. In practices these alleles are often different by only 1 or 2
AA substitutions and the impact of choosing one of two similar alleles
on the final pairwise distance metric is small (i.e < 10 tcrdist units)

The test here compares 267 sequences that have unambiguous v-gene calls
where the tcrdist2 and tcr-dist legacy outputs should match.

truth matrix - is a distance matrix output by legacy tcrdist
python compute_distances.py --clones_file  mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv --organism mouse

mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_AB.dist

"""

from tcrdist.tests.legacy_truth_values import truth_matrix
from tcrdist.tests.legacy_truth_values import test_df

import unittest
import pandas as pd
import numpy as np
import tcrdist as td
from tcrdist import tcr_distances
from tcrdist import mappers
from tcrdist import objects
from tcrdist.repertoire import TCRrep

class test_comparability(unittest.TestCase):

    def test_reproduction_of_tcrdist_metric(self):

        pb1 = TCRrep(cell_df = test_df, organism = "mouse")             # 6
        pb1.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned = True) # <------- # 7
        pb1.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned = True) # <------- # 8
        pb1.index_cols =['epitope',                                      # 9
                        'clone_id',
                        'subject',
                        'cdr3_a_aa',
                        'cdr1_a_aa',
                        'cdr2_a_aa',
                        'pmhc_a_aa',
                        'cdr3_b_aa',
                        'cdr1_b_aa',
                        'cdr2_b_aa',
                        'pmhc_b_aa',
                        'v_a_gene',
                        'v_b_gene']
        pb1.deduplicate()                                                # 10
        pb1.compute_pairwise_all(chain = "alpha",                        # <------- 11
                                 metric = 'tcrdist_cdr3',
                                 compute_specific_region = 'cdr3_a_aa',
                                 #user_function = tcrdist_metric_align_cdr3s_false,
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_a_aa',
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_a_aa',
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_a_aa',
                                 processes = 6)

        pb1.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = 'tcrdist_cdr3',
                                 #user_function = tcrdist_metric_align_cdr3s_false,
                                 compute_specific_region = 'cdr3_b_aa',
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_b_aa',
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_b_aa',
                                 processes = 6)
        pb1.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_b_aa',
                                 processes = 6)

        pb1.compute_paired_tcrdist()

        assert(np.all(pb1.paired_tcrdist.shape == truth_matrix.shape))

        test_order = pb1.clone_df.copy()
        test_order['sort_index'] = test_order.index
        test_order.index = test_order.clone_id
        test_order_sorted = test_order.loc[pb1.cell_df.clone_id,]
        sort_index = test_order_sorted['sort_index']

        # Sort so that  clone_df matches input cell_df order
        # truth matric is based on order in original_df

        test_matrix = pd.DataFrame(pb1.paired_tcrdist).loc[sort_index, sort_index].values
        self.assertTrue(np.all(test_matrix == truth_matrix))



if __name__ == '__main__':
    unittest.main()
