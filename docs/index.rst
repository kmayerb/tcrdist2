tcrdist2
========

2020-06-27

**tcrdist2 provides flexible distance measures for comparing T cell receptors across complementarity determining regions (CDRs)**

`tcrdist2 <https://github.com/kmayerb/tcrdist2>`_ is a python API-enabled toolkit for analyzing T-cell receptor repertoires. Some of the functionality and code is adapted from the original tcr-dist package which was released with the publication of Dash et al. Nature (2017). The manuscript
(`doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_). This package contains a new API for accessing the features of tcr-dist, as well as many new features that expand the T cell receptor analysis pipeline.


.. image:: E3.jpg

**Extended Data Figure 3 : Schematic overview of the TCRdist calculation.**

Extended Data `Figure 3 <https://www.nature.com/articles/nature22383/figures/7>`_
is particularly illustrative.

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
   Inputs
   tcrdist2

.. _quick-examples:
  
.. toctree::
   :caption: Quick Examples
   :maxdepth: 2
   
   ex1


.. _longk-examples:

.. toctree::
   :caption: Long Examples
   :maxdepth: 2
   
   HotStart
   CompleteExample
   Saving
   Plotting

.. _Pre-Processing:

.. toctree::
   :caption: Pre-Processing
   :maxdepth: 2

   Mixcr
   Preprocess_10X
   Preprocess_Adaptive
   Preprocessing
   tcrdist1_to_tcrdist2

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
