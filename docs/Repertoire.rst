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


.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: tcrdist.repertoire

.. autoclass:: TCRrep
    :members: infer_cdrs_from_v_gene, deduplicate, compute_pairwise_all, compute_paired_tcrdist, infer_olga_aa_cdr3_pgens, reduce_file_size, save_as_hdf5, rebuild_from_hdf5
    :show-inheritance:
