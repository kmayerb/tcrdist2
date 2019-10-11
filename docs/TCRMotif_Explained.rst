
TCRMotif_Explained
==================

October 10, 2019

`TCRMotif` is a tcrdist2 class for running the find_cdr3_motif_ from Dash et al. 
in tcrdist2.

.. _find_cdr3_motif: https://github.com/phbradley/tcr-dist/blob/master/find_cdr3_motifs.py



This the most time-consuming process currently incorporated into tcrdist2.
It can take 10-15 minutes per epitope with 1000-5000 clones.

In tcrdist2, candidate CDR3 motifs are returned as a DataFrame rather than a
log-file.

Using TCRMotif Directly on a Clones File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import pandas as pd
    from tcrdist.cdr3_motif import TCRMotif
    fn = "mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv"
    clones_df_test = pd.read_csv(fn, sep="\t")
    tm = TCRMotif(clones_df = clones_df_test,
                  organism = "mouse",
                  chains = ["A","B"],
                  epitopes =["PA"])
    discovered_motifs = tm.find_cdr3_motifs()

Changing TCRMotif Defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are number of parameters in the original tcrdist
``find_cdr3_script.py``.

The defaults can be viewed by

.. code:: python

   tm = TCRMotif().params
   tm.params

.. code:: bash

   {'min_count': 10,
    'max_ng_lines': 5000000,
    'max_motif_len': 100,
    'nsamples': 25,
    'min_expected': 0.25,
    'max_overlap': 0.8,
    'max_useful_expected': 1.0,
    'chi_squared_threshold': 50.0,
    'verbose' : false
    ...
    }

These can all be modified by directly changing the class attributes. e.g.,

.. code:: python

   tm.verbose = True

The updated parameters will be used at runtime ``TCRMotif.find_cdr3_motifs()``. e.g.,


.. code:: ipython3

    import pandas as pd
    from tcrdist.cdr3_motif import TCRMotif

    clones_df_test = pd.read_csv("mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv", sep="\t")
    tm = TCRMotif(clones_df = clones_df_test, organism = "mouse", chains = ["A","B"], epitopes =["PA"])
    tm.very_verbose = True
    discovered_motifs = tm.find_cdr3_motifs()

Combining TCRrep and TCRMotif
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TCRMotif can be combined with the tcrdist2 interface. Here is an example
of a motif scan on data take from vdjDB

.. code:: ipython3

    import pandas as pd

    import tcrdist as td
    from tcrdist import mappers
    from tcrdist.repertoire import TCRrep
    from tcrdist.cdr3_motif import TCRMotif

    pd_df = pd.read_csv("vdjDB_PMID28636592.tsv", sep = "\t")        # 1
    t_df = td.mappers.vdjdb_to_tcrdist2(pd_df = pd_df)               # 2
    t_df.organism.value_counts                                       # 3
    index_mus = t_df.organism == "MusMusculus"                       # 4
    t_df_mus = t_df.loc[index_mus,:].copy()                          # 5

    tr = TCRrep(cell_df = t_df_mus, organism = "mouse")              # 6

    tr.infer_cdrs_from_v_gene(chain = 'alpha')                       # 7
    tr.infer_cdrs_from_v_gene(chain = 'beta')                        # 8

    tr.index_cols = ['subject', 'epitope',                           # subject and epitope
                     'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',# gene usage, not this essential info for TCRMotif
                     'cdr3_a_aa', 'cdr3_b_aa',                       # CDR 3
                     'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',          # alpha CDR 1, 2, and 2.5
                     'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']          # beta CDR 1, 2, and 2.5

    tr.deduplicate()                                                 # 10

    tm = TCRMotif(clones_df = tr.tcr_motif_clones_df(), organism = "mouse", chains = ["A","B"], epitopes = ["PA"]) # 11
    discovered_motifs = tm.find_cdr3_motifs()
    tm.motifs_df.head()

From Raw Paired Nucleotide Sequences to Motifs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tcrdist2 can produce the 'clones_df' DataFrame from raw paired end sequences.
This is not recommended for large dataset > 1000 paired sequeneces
because its very time consuming, but here is an example for reference.

Also note that identifyClones function may be more stringent compared
with the deduplicate method in tcrdist2.

.. code:: ipython3

    import pandas as pd

    import tcrdist as td
    from tcrdist.cdr3_motif import TCRMotif
    ps_df = td.processing.readPairedSequences(paired_seqs_file = "tcrdist/datasets/test_mouse_pairseqs.tsv",
                                              organism = "mouse", use_parasail= False)
    prob_df = td.processing.computeProbs(ps_df)
    assert prob_df.shape[0] == ps_df.shape[0]
    ps_prob_df = pd.concat([ps_df, prob_df], axis=1)
    clones_df = td.processing.identifyClones(ps_prob_df, min_quality_for_singletons=0)
    clones_df.head()
    tm_ex = TCRMotif(clones_df = clones_df, organism = "mouse", chains = ["A","B"], epitopes = ["PA"])
    tm_ex.find_cdr3_motifs()
