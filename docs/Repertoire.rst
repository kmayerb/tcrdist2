Repertoire
==========

This page documents the class :py:class:`tcrdist.repertoire.TCRrep` within the
:py:mod:`tcrdist.repertoire` module
(`repertoire.py <https://github.com/kmayerb/tcrdist2/blob/API2/tcrdist/repertoire.py>`_ source code.)

The application programming interface for tcrdist2 involves a set of step-wise
commands centered around the :py:class:`tcrdist.repertoire.TCRrep` class, which
stores input data, methods, user-parameters, and results. It calls methods for
parallelized alignment and distance computation between amino acid strings in the
:py:mod:`tcrdist.pairwise` module
(`pairwise.py <https://github.com/kmayerb/tcrdist2/blob/API2/tcrdist/pairwise.py>`_ source code.)

.. code-block:: python

  import pandas as pd
  import tcrdist as td
  from tcrdist import mappers
  from tcrdist.repertoire import TCRrep


  pd_df = pd.read_csv("DMJVdb_PMID28636592.tsv", sep = "\t")      # 1
  t_df = td.mappers.vdjdb_to_tcrdist2(pd_df = pd_df)              # 2
  t_df.organism.value_counts                                      # 3
  index_mus = t_df.organism == "MusMusculus"                      # 4
  t_df_mus = t_df.loc[index_mus,:].copy()                         # 5


  tr = TCRrep(cell_df = t_df_mus, organism = "mouse")             # 6

  tr.infer_cdrs_from_v_gene(chain = 'alpha')                      # 7

  tr.infer_cdrs_from_v_gene(chain = 'beta')                       # 8

  tr.index_cols =['epitope',                                      # 9
                  'subject',
                  'cdr3_a_aa',
                  'cdr1_a_aa',
                  'cdr2_a_aa',
                  'pmhc_a_aa',
                  'cdr3_b_aa',
                  'cdr1_b_aa',
                  'cdr2_b_aa',
                  'pmhc_b_aa']

  tr.deduplicate()                                                # 10


  tr.compute_pairwise_all(chain = "alpha",                        # 11
                          metric = "hamming",
                          processes = 6)

  tr.compute_pairwise_all(chain = "beta",                         # 12
                          metric = "hamming",
                          processes = 6)

  tr.compute_paired_tcrdist(chains = ['alpha','beta'],            # 13
                            replacement_weights= {'cdr3_a_aa_pw': 3,
                                                  'cdr3_b_aa_pw': 3})




.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: tcrdist.repertoire

.. autoclass:: TCRrep
    :members: infer_cdrs_from_v_gene, deduplicate, compute_pairwise_all, compute_paired_tcrdist
    :show-inheritance:
