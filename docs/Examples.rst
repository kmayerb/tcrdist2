Examples
========


These vignette use tcrdist2 to explore real TCR data to demonstrate some of the
core features of tcrdist2.


Example 1: Compute Distances for paired alpha/beta chain sequences from Dash et al. 2017
----------------------------------------------------------------------------------------

The data for this example comes from the paper "Quantifiable predictive features
define epitope-specific T cell receptor repertoires", which can be downloaded
from the `VDJdb <https://vdjdb.cdr3.net/search>`_. Data were selected using the
search filter (references == PMID:28646592). This query returns CDR3 sequences
and predicted V and J gene usage for 2442 paired alpha/beta T-cell receptors
with known epitope specificities. These were downloaded as a tab separated flat
file (*DMJVdb_PMID28636592.tsv*) which can be downloaded here [URL]

The application programming interface for tcrdist2 involves a set of step-wise
commands centered around the TCRrep class, which stores input data,
user-parameters, and results.

Before explaining each step in detail, it is useful to present these steps as a
block of stepwise code. Below is the entire set of commands needed to calculate
a tcrdist - a weighted pairwise distance measure based on comparison across
multiple T-cell receptor complementarity-determining regions (CDRs),
In addition to annotations describing each step, we follow this
introductory example with more information on alternative analysis configurations.
This page illustrates the use of tcrdist2.


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



Stepwise Explanation
^^^^^^^^^^^^^^^^^^^^

#. Load data to a Pandas DataFrame.

#. Use `td.mappers.vdjdb_to_tcrdist2` to convert data from vdjDB to a Pandas
   DataFrame with correct tcrdist2 headers.

#. Count instances of human and mouse sequences in the data.

#. Index the sequences that come from MusMusculus (mouse).

#. Create a copy of the subset DataFrame `t_df_mus`, including only mouse TCRs.

#. Create an instance of the TCRrep class initialized with the `t_df_mus` DataFrame.
   Upon initialization organism is set to "mouse" and formatted input data is passed to `cell_df` argument

#. Use `infer_cdrs_from_v_gene` to populate CDR1, CDR2 and pMHC loop fields fo the 'alpha' chain.
    - `chain` argument is set to either 'alpha', 'beta', 'delta', 'gamma'

#. Repeat step 7, setting chain to 'beta'.
    - Because of hypermutation that occurs in the CDR3 region, the CDR3 sequence
      must be directly supplied. However, for the other complementarity-determining
      regions the sequence is not provided in the vdjDB data. However, tcrdist2
      uses the predicted v-gene variant call to infer the likely amino acid
      sequence at the remaining complementarity-determining regions: CDR1, CDR2,
      and the pMHC loop positions (the pMHC loop is between the CDR2 and CDR3).

#. Specify index columns. Any sequence identical across all the index columns
   will be grouped at the following step.

#. Call `tr.deduplicate` to remove duplicates and create the `tr.clones_df` attribute.
   *Even if there are no duplicates this step is necessary to produce the `clones_df` attribute of the TCRrep instance.* Any row of the DataFrame missing any of the CDRs specified in the `index_col` list will not be included in the `clones_df` DataFrame.

#. Call `tr.compute_pairwise_all` specifying the chain, metric, and number of parallel processes to use
    - `chain` argument is set to either 'alpha', 'beta', 'delta', 'gamma'
    - `metric` argument is set to either 'hamming', 'nw' or 'custom' In this
      first example we are using the Hamming Distance, which is the number of
      mismatching positions between two aligned strings. We will later show
      how tcrdist2 can incorporate the use of substitution matrices.
    - `processes` argument is optional for specifying the number of
      available CPUs. tcrdist2 uses python's multiprocessing package to
      parallelize pairwise distance computation.

#. Repeat the previous step setting `chain` argument to 'beta'. We will show how
   individual CDR computations can be specified in a later example.

#. Once the site-specific individual pairwise distances are computed across all
   of the complementarity determining regions, `tr.compute_paired_tcrdist` computes
   the 'tcrdist'- a weighted sum of the distances at each of the CDRs.
  - The default to weight all CDRs equally but the argument `replacement_weights`
    takes a dictionary and can be used to place greater weight on
    mismatches occurring in certain CDRs

That's it! If you've followed along you've computed a tcrdist from real data.


Examining the Results
^^^^^^^^^^^^^^^^^^^^^

Without considering the extensive plotting tools developed in the
original version of TCRdist. Let's take quick look at the results.
You might be wondering how 'tcrdistances' are distributed among receptors
with shared or distinct epitope specificities.


Clustered Heatmap
^^^^^^^^^^^^^^^^^

.. image:: f2.png

The code for producing this figure directly from tcrdist output is shown below.


Example 1A: Accessing Individual CDR Results
--------------------------------------------

Here we extend the first example to show the flexibility of tcrdist2. In the
introductory example, we counted the number of mismatches
between 8 total CDRs and combined the results into a single distance metric.

This individual information is readily available within the
instance of the TCRrep class that was created.

A common naming convention is used to store a number of objects within the TCRrep class.

TCRrep.[cdr1|cdr2|cdr3|pmhc]_[a|b|d|g]_aa_pw

- the first position references the CDR.

- the second position references a: alpha, b: beta, d: delta, g: gamma chains

- the third position references the molecular type aa: amino acid or nuc: nucleotide

- the final position reference the object pw: pairwise, sm: substitution matrix, etc.


To examine only the the cdr3 regions, the pairwise results can be directly accessed:

.. code-block:: python

  tr.cdr3_a_aa_pw
  tr.cdr3_b_aa_pw

These can be analyzed directly. Moreover, the `tcrdistances` can be
recalculated with different CDR weights. By default the last generated
tcrdist is stored as TRCrep.paired_tcrdist. However, multiple results can be
stored in the store_tcrdist attribute. The following example illustrates the point.

Example 1B: Custom Weights and Stored Results
---------------------------------------------

.. code-block:: python

  tcrdist0 = tr.compute_paired_tcrdist(chains = ['alpha','beta'])

  replacement_weights = {'cdr1_a_aa_pw':1,
                         'cdr2_a_aa_pw':1,
                         'cdr3_a_aa_pw':2,
                         'pmhc_a_aa_pw':1,
                         'cdr1_b_aa_pw':2,
                         'cdr2_b_aa_pw':2,
                         'cdr3_b_aa_pw':4,
                         'pmhc_b_aa_pw':0}
  # 1
  tcrdist1 = tr.compute_paired_tcrdist(chains = ['alpha','beta'],
                          replacement_weights= replacement_weights)
  # 2
  tcrdist0 = tr.compute_paired_tcrdist(chains = ['alpha','beta'])


  # 3
  tr.stored_tcrdist.append(tcrdist0)
  tr.stored_tcrdist.append(tcrdist1)

  # 4
  tr.stored_tcrdist[1]



#. Repeat step 13 from the previous example using the default weights of 1

#. Repeat step 13 using new weights.

#. Store both results by appending them to `sotred_tcrdist` list attribute

#. Access either result. The weights are stored along with the pairwise distances.


    {'paired_tcrdist': array([[  0.,  76.,  80., ...,  89.,  89.,  87.],
          [ 76.,   0.,  60., ...,  81.,  75.,  43.],
          [ 80.,  60.,   0., ...,  59.,  81.,  77.],
          ...,
          [ 89.,  81.,  59., ...,   0.,  60.,  58.],
          [ 89.,  75.,  81., ...,  60.,   0.,  40.],
          [ 87.,  43.,  77., ...,  58.,  40.,   0.]]),
    'paired_tcrdist_weights': {'cdr1_a_aa_pw': 1,
    'cdr1_b_aa_pw': 2,
    'cdr2_a_aa_pw': 1,
    'cdr2_b_aa_pw': 2,
    'cdr3_a_aa_pw': 2,
    'cdr3_b_aa_pw': 4,
    'pmhc_a_aa_pw': 1,
    'pmhc_b_aa_pw': 2}}


Example 1C: Using a Substitution Matrix
----------------------------------------

The introductory example used the Hamming Distance (number of aligned positions
with mismatching information) to calculate pairwise distance between each CDR3.

In the original investigation “Quantifiable predictive features define
epitope-specific T cell receptor repertoires”, took a different approach.

    The TCRdist distance between two TCRs is defined to be the similarity-weighted mismatch distance between the potential pMHC-contacting loops of the two receptors (Extended Data Fig. 3). The loop definitions used are based on the IMGT CDR definitions (http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html) with the following modifications: (1) we include the pMHC-facing loop between CDR2 and CDR3 (IMGT alignment columns 81–86) since residues in this loop have been observed making pMHC contacts in solved structures; (2) we use the ‘trimmed CDR3’ defined above rather than the full IMGT CDR3. The mismatch distance is defined based on the BLOSUM62 (ref. 37) substitution matrix as follows: distance (a, a) = 0; distance (a, b) = min (4, 4-BLOSUM62 (a, b)), where 4 is 1 unit greater than the most favourable BLOSUM62 score for a mismatch, and a and b are amino acids. This has the effect of reducing the mismatch distance penalty for amino acids with positive (that is, favourable) BLOSUM62 scores (for example,: dist(I, V) = 1; dist(D, E) = 2; dist(Q, K) = 3), where I, V, D, E, Q and K are the single letter amino acid codes for isoleucine, valine, aspartate, glutamate, glutamine and lysine, respectively. A gap penalty of 4 (8 for the CDR3) is used as the distance between a gap position and an amino acid. To account for the greater role of the CDR3 regions in peptide recognition and offset the larger number (3) of non-CDR3 loops, a weight of 3 is applied to mismatches in the CDR3s.

    For each epitope-specific repertoire, we computed a TCRdist distance matrix between all receptors. This distance matrix was used for clustering and dimensionality reduction as described below as well as in the TCRdiv diversity calculation. The sampling density nearby each receptor was estimated by taking the weighted average distance to the nearest-neighbour receptors in the repertoire: a small nearest-neighbours distance (NN-distance) indicates that there are many other nearby receptors and hence greater local sampling density. For analyses reported here we used the nearest 10 per cent of the repertoire with a weight that linearly decreases from nearest to farthest neighbours. Values smaller than 10 focus on the very nearest neighbours, enhancing detection of rare clusters, while increasing the sensitivity to noise or... *




















Additional Code for Plots
-------------------------


Code For Clustered Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

  # Convert the ndarray containing tcrdist to a DataFrame
  px = pd.DataFrame(tr.paired_tcrdist)


  g = sns.clustermap(data= px,
                     row_colors = row_colors,
                     col_colors = row_colors,
                     row_cluster=True,
                     col_cluster=True,
                     yticklabels=False,
                     xticklabels=False,
                    )

  # bostock3 is Mike Bostock's 3rd Categorical Set
  bostock3 = ["#8dd3c7","#ffffb3","#bebada","#fb8072",
              "#80b1d3","#fdb462","#b3de69","#fccde5",
              "#d9d9d9","#bc80bd","#ccebc5","#ffed6f"]
  lut = dict(zip(tr.clone_df.epitope.unique(), bostock3))
  row_colors = tr.clone_df.epitope.map(lut)

  # set legend
  for label in tr.clone_df.epitope.unique():
      g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                              label=label, linewidth=0)

  g.ax_row_dendrogram.set_visible(False)
  g.ax_col_dendrogram.legend(loc="center", ncol = 4)
  g.cax.set_position([.97, .2, .03, .45])

Distribution of Distances
^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: f1.png

Code For Distribution of Distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Python

  import matplotlib
  import matplotlib.pyplot as plt
  import seaborn as sns
  %matplotlib inline


  def epitope_to_epitope(e1,
                         e2,
                         clone_df = tr.clone_df,
                         paired_tcrdist = tr.paired_tcrdist,
                         var = "epitope"):
    """
    A function for subsetting distances to TCRs with shared or distinct or
    shared epitope specificity.
    """
    e1_ind = clone_df[var] == e1
    e2_ind = clone_df[var] == e2
    tr_df = pd.DataFrame(paired_tcrdist)
    e1_to_e2 = tr_df.loc[e1_ind , e2_ind].values.flatten()
    return(e1_to_e2)

  sns.kdeplot(epitope_to_epitope(e1 = "M45", e2 = "M45"), bw = 4, label = "tcrdist(M45,M45)")
  sns.kdeplot(epitope_to_epitope(e1 = "PB1", e2 = "PB1"), bw = 4, label = "tcrdist(PB1,PB1)")
  sns.kdeplot(epitope_to_epitope(e1 = "M45", e2 = "PB1"), bw = 4, label = "tcrdist(M45,PB1)")
  plt.legend(loc = 2);
