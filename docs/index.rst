tcrdist2
========

`tcrdist2 <https://github.com/kmayerb/tcrdist2>`_ is a python
API-enabled toolkit expanding on the T cell receptor analysis pipeline
developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. Nature (2017). The manuscript
(`doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_)
outlines a basis for T cell receptor comparison by
incorporating dissimilarity across multiple complementarity determining regions.
Extended Data `Figure 3 <https://www.nature.com/articles/nature22383/figures/7>`_
is particularly illustrative.

.. image:: E3.jpg

**Extended Data Figure 3 : Schematic overview of the TCRdist calculation.**


"Each TCR is mapped to the amino acid sequences of the loops within the
receptor that are known to provide contacts to the pMHC (commonly referred
to as CDR1, CDR2, and CDR3, as well as an additional variable loop between
CDR2 and CDR3). The distance between two TCRs is computed by comparing these
concatenated CDR sequences using a similarity-weighted Hamming distance,
with a gap penalty introduced to capture variation in length and a higher
weight given to the CDR3 sequence in recognition of its disproportionate
role in epitope specificity (see Methods and Extended Data Fig. 3)."


.. _getting-started:


.. toctree::
   :caption: Getting started
   :maxdepth: 2

   Installation
   HotStart
   Inputs
   CompleteExample
   Saving
   Plotting
   

.. toctree::
   :caption: Pre-Processing
   :maxdepth: 2

   Mixcr
   Preprocess_10X
   Preprocess_Adaptive
   Preprocessing
   tcrdist1_to_tcrdist2

.. _more-details:

.. toctree::
   :caption: More details
   :maxdepth: 2
   
   Motifs
   Mappers
   Neighborhood
   OlgaModel
   Pairwise
   PairwiseDistance
   pGen
   Processing
   Proximity 
   Repertoire
   Storage
   Subset
   TCRMotif_Explained
   WhatIsNew


.. _in-depth-topics:

.. toctree::
   :caption: In-Depth
   :maxdepth: 2

   Example_2
   Example_1






Lookups
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
