.. tcrdist2 documentation master file, created by
   sphinx-quickstart on Mon Jul  8 11:24:44 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

tcrdist2
========

tcrdist2 is a API-enabled toolkit for T-cell receptor analysis inspired by the
TCRdist pipeline developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. Nature (2017).
`doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_

The original investigation “Quantifiable predictive features define
epitope-specific T cell receptor repertoires”, outlines an approach based
on incorporating multiple complementarity determining regions as a
basis for receptor comparison.

    "Each TCR is mapped to the amino acid sequences of the loops within the
    receptor that are known to provide contacts to the pMHC (commonly referred
    to as CDR1, CDR2, and CDR3, as well as an additional variable loop between
    CDR2 and CDR3). The distance between two TCRs is computed by comparing these
    concatenated CDR sequences using a similarity-weighted Hamming distance,
    with a gap penalty introduced to capture variation in length and a higher
    weight given to the CDR3 sequence in recognition of its disproportionate
    role in epitope specificity (see Methods and Extended Data Fig. 3)."


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Examples
   Repertoire
   Pairwise
   Processing





Lookups
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
