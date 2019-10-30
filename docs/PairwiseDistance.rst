.. _PairwiseDistance:

Pairwise Distance
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

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
