Complete Example
================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
############

This page illustrates the integration of the major features of tcrdist2.

There are three major classes in tcrdist2:

* **TCRrep** :py:class:`tcrdist.repertoire.TCRrep` - specifies metric and computes distances between members of a TCR repertoire

* **TCRsubset** :py:class:`tcrdist.subset.TCRsubset` - analyzes epitope specificity for a specific subset of a TCR repertoire

* **TCRmotif** :py:class:`tcrdist.cdr3_motif.TCRMotif` - identifies candidate motifs characteristic of a subset of receptors

In this example, we emphasize how several tcrdist2 modules are integrated.
We reserve a more detailed expositio for separate sections,
each with its own detailed page in this documentation.

Before proceeding, we also offer a note on the
syntax used in the following coded examples. We instantiate
instances of **TCRrep**, **TCRsubset**, and **TCRMotif** classes as
`tr`, `ts`, and `tm`, respectively.

Load tcrdist2
#############

.. code-block:: python

  import pandas as pd
  import numpy as np
  import tcrdist as td

  from tcrdist import mappers
  from tcrdist.repertoire import TCRrep
  from tcrdist.cdr3_motif import TCRMotif
  from tcrdist.subset import TCRsubset
  from tcrdist.storage import StoreIOMotif, StoreIOEntropy


Preprocessing
#############

Example Data
************
This example uses the mouse clones files from the Dash et al. 2017 study.
Clones file can be downloaded (`clones <https://www.dropbox.com/s/l0z12f8lc752wfx/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv?dl=1>`_ 2MB).
The search for motifs in tcrdist2 is a slow step. Users who wish to use
precomputed candidate motifs can download a motifs file
(`PA motifs <https://www.dropbox.com/s/z7wwmwb1n6dpq74/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log?dl=1>`_ 32KB).



.. code-block:: python

  #1
  tcrdist_clone_fn = 'mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
  tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")

  # 2
  ind = (tcrdist_clone_df.epitope == "PA") | (tcrdist_clone_df.epitope == "F2")
  tcrdist_clone_df = tcrdist_clone_df[ind].copy()

  #3
  mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping
  #4
  tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,
                                              mapping = mapping)

The steps shown above are:

1. Load the clones file as a Pandas DataFrame.
2. Select only those clones specific to epitope F2 or PA. We choose two here, so we can later illustrate subsetting by epitope specificity.
3. load an appropriate mapping (i.e. a dictionary of old and new column names).
4. Apply the mapping to the the original DataFrame, to select and rename columns.


Pairwise Distance
#################

TCRrep defines metric for comparing T Cell receptors.

.. code-block:: python

  #1
  tr = TCRrep(cell_df = tcrdist2_df, organism = "mouse")

  #2
  tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
  tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)

  #3
  tr.index_cols = ['clone_id', 'subject', 'epitope',
                   'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                   'cdr3_a_aa', 'cdr3_b_aa',
                   'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                   'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                   'cdr3_b_nucseq', 'cdr3_a_nucseq',
                   'va_countreps', 'ja_countreps',
                   'vb_countreps', 'jb_countreps',
                   'va_gene', 'vb_gene',
                   'ja_gene', 'jb_gene']

  #4
  tr.deduplicate()

  #5
  tr._tcrdist_legacy_method_alpha_beta()

  #6
  distA = tr.dist_a
  distB = tr.dist_b
  assert np.all(((distA + distB) - tr.paired_tcrdist) == 0)

The steps shown above are:

1. Instantiate an instance of TCRrep class: `tr`.
2. Infer CDR2, CDR2, and CDR2.5. Setting imgt_aligned = True is crucial if using the legacy TCRdist metric.
3. Set index columns.
4. Remove duplicates (duplicates are clones identical across all index columns).
5. `_tcrdist_legacy_method_alpha_beta()` is a convenience function for calculating a legacy TCRdist in the method described by Dash et al. 2017.
6. Save alpha- and beta-specific distance matrices.

Note: CDR2.5 alpha and CDR2.5 beta are referred to in tcrdist2 as pmhc_a and phmc_b, respectively.

tcrdist2 offers a lot flexibility in the calculation of "tcrdistances,"
and the reader is encouraged to consult the X,Y, and Z documentation page.
The specificiation of `_tcrdist_legacy_method_alpha_beta()` using tcrdist2
can be illuminated by typing `??TCRrep._tcrdist_legacy_method_alpha_beta`
in ipython.

Subset
######

TCRsubset defines a subset of interest.

.. code-block:: python

  # 1
  criteria = tr.clone_df.epitope == "PA"
  clone_df_subset = tr.clone_df[criteria]

  # 2
  distA_subset = distA.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()
  distB_subset = distB.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()

  # 3
  ts = TCRsubset(clone_df_subset,
               organism = "mouse",
               epitopes = ["PA"],
               epitope = "PA",
               chains = ["A","B"],
               dist_a = distA_subset,
               dist_b = distB_subset)

The steps shown above are:

1. Set Some Criteria for the sequences we want to look at.
2. Select subsets of the original distance matrices
3. Instantiate an instance of TCRsubset: `ts`


Motifs
######

Searchs for Motifs
******************

Warning this step can take a long time! Save or serialize results to disk.

.. code-block:: python

  tm = TCRMotif(ts.tcr_motif_clones_df(),    #1
               organism = "mouse",
               chains = ["A"],
               epitopes = ["PA"])
  tm.find_cdr3_motifs()                      #2

The steps shown above are:

1. Instantiate an instance of TCRMotif,  `tm`
2. Search for motifs.


Analyze Motifs
##############

Load Motifs From File
*********************

If using a previously generated motifs file:

.. code-block:: python

  cnames = ["file_type","count", "expect_random","expect_nextgen", "chi_squared", "nfixed",
  "showmotif", "num", "othernum", "overlap", "ep", "ab",
  "nseqs", "v_rep_counts", "j_rep_counts"]

  motif_fn = 'mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log'
  x = open(motif_fn, "r").readlines()
  motif_df = pd.DataFrame([l.split() for l in x], columns = cnames)
  i = 0
  row = motif_df.iloc[i,:].to_dict()


Load Motifs From a TCRMotif Instance
************************************

If using a motifs generated by TCRMotif:

.. code-block:: python


  i = 0
  row = tm.motif_df.iloc[i,:].to_dict()


Analyzing Motifs
****************

.. code-block:: python

  from tcrdist.plotting import plot_pwm
  from tcrdist.storage import StoreIOMotif, StoreIOEntropy
  import pandas as pd
  import IPython

  # 1
  StoreIOMotif_instance = StoreIOMotif(**row)
  # 2
  StoreIOMotif_instance = ts.analyze_motif(s = StoreIOMotif_instance)
  # 3
  StoreIOMotif_instance = ts.analyze_matches(s = StoreIOMotif_instance)

Notice that analysis of motifs methods which handle the comparisons of
epitope specific cdr3 to non-epitope specific reference datasets are
containted within TCRsubset.

The steps shown above are:

1. Initialize instance of the information carrier class StoreIOMotif
2. Analyze_motif (`ts.analyze_motif`) to determine matches and neighbors of matches,  with new attributes are appended to the StoreIOMotifinstance.
3. Analyze matches (`ts.analyze_matches`) to identify the the relative entropy between position wise matrices


Visualizing Motifs
******************

.. code-block:: python

  svg = plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600)
  IPython.display.SVG(svg)

.. image:: PA_alpha_MotifPlot0.svg


We plotted just on of many candidate motifs.

All suggested motifs, or a subset of motifs, can be plotted from a the motifs
DataFrame, as is illustrated below:

.. code-block:: python

  svg_list = []
    for i,row in motif_df.iterrows():
      StoreIOMotif_instance = StoreIOMotif(**row)
      StoreIOMotif_instance = ts.analyze_motif(  s = StoreIOMotif_instance)
      StoreIOMotif_instance = ts.analyze_matches(s = StoreIOMotif_instance)
      svg = plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600)
      svg_list.append(svg)

  IPython.display.SVG(svg_list[1])
  IPython.display.SVG(svg_list[2])
