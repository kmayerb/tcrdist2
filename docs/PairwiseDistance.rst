.. _PairwiseDistance:

Pairwise Distance
=================

Introduction
############



Pairwise Distance

.. tip ::
  Amino acid substitution matrices (i.e., BLOSUM62) are used by tcrdist2
  to find the highest scoring alignment between two CDRs. Once the alignment
  is determined different tcrdistance metrics exist for estimating a distance
  between the two sequences.

.. tip ::
  'Legacy tcrdist' refers to the algorithms used in the Dash et al. 2016 paper:
  "The mismatch distance is defined based on the BLOSUM62 substitution matrix
  as follows: distance (a, a) = 0; distance (a, b) = min (4, 4-BLOSUM62 (a, b)),
  where 4 is 1 unit greater than the most favourable BLOSUM62 score for a mismatch,
  and a and b are amino acids.
  This has the effect of reducing the mismatch distance penalty for amino acids
  with positive (that is, favourable) BLOSUM62 scores (for example,:
  dist(I, V) = 1; dist(D, E) = 2; dist(Q, K) = 3), where I, V, D, E, Q and K
  are the single letter amino acid codes for isoleucine, valine, aspartate,
  glutamate, glutamine and lysine, respectively.
  A gap penalty of 4 (8 for the CDR3) is used as the distance between a gap
  position and an amino acid. To account for the greater role of the CDR3
  regions in peptide recognition and offset the larger number (3)
  of non-CDR3 loops, a weight of 3 is applied to mismatches in the CDR3s."



TCRdist for Prediction
######################

.. code-block:: python

  import pandas as pd
  import numpy as np
  import tcrdist as td

  from tcrdist import mappers
  from tcrdist.repertoire import TCRrep

  tcrdist_clone_fn = 'mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
  tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")

  mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping
  tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,
                                              mapping = mapping)


  tr = TCRrep(cell_df = tcrdist2_df, organism = "mouse")


  tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
  tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)


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


  tr.deduplicate()

  tr._tcrdist_legacy_method_alpha_beta()

  distA = tr.dist_a
  distB = tr.dist_b
  assert np.all(((distA + distB) - tr.paired_tcrdist) == 0)



  # K NEAREST NEIGHBORS
  pr = namedtuple("perf", ["observed", "predicted", "dist"])
  obsereved = tr.clone_df.epitope.to_list()
  performance = list()

  k = 5
  for i,row in tr.clone_df.iterrows():
      ind = (tr.clone_df.subject != row.subject)                          # Index hold out all data from that subject
      distances = tr.paired_tcrdist[i,ind]                                # Get Distances from the ith row, holding out subject
      sorted_indices = np.argsort(distances)                              # Get index of storted distances small to large
      sorted_epitopes = tr.clone_df.epitope.iloc[sorted_indices].to_list()# Get epitopes associated wtih those indices
      sorted_distances =  distances[sorted_indices]                       # Get distances associated with those neighbors
      predicted = sorted_epitopes[0:k]                                    # Get Predicted epitopes for K nearest neighbors
      predicted_distance = sorted_distances[0:k]                          # Get distances for K nearest neighbots
      performance.append(pr(obsereved[i], predicted, predicted_distance)) # Save Performance Information

  performance[1:10]
  #  [perf(observed='PA', predicted=['PA', 'm139', 'NP', 'PB1', 'M45'], dist=array([132., 150., 153., 165., 171.])),
  #  perf(observed='PA', predicted=['m139', 'PA', 'PA', 'PA', 'PA'], dist=array([36., 42., 48., 60., 66.])),
  #  perf(observed='PA', predicted=['PA', 'PA', 'PA', 'PA', 'PA'], dist=array([33., 38., 42., 42., 54.])),
  #  perf(observed='PA', predicted=['NP', 'PA', 'PA', 'PA', 'PA'], dist=array([84., 84., 90., 90., 96.])),
  #  perf(observed='PA', predicted=['PA', 'M38', 'PA', 'PA', 'PA'], dist=array([84., 84., 84., 84., 84.])),
  #  perf(observed='PA', predicted=['PA', 'PA', 'PA', 'PA', 'PA'], dist=array([24., 30., 32., 42., 45.])),
  #  perf(observed='PA', predicted=['PA', 'M45', 'PA', 'PA', 'PA'], dist=array([ 0., 53., 54., 56., 56.])),
  #  perf(observed='PA', predicted=['NP', 'PA', 'PA', 'PA', 'PA'], dist=array([116., 116., 116., 126., 126.])),
  #  perf(observed='PA', predicted=['PA', 'PA', 'PA', 'PA', 'PA'], dist=array([15., 32., 33., 51., 51.]))]
