Pairwise
========


The application programming interface for tcrdist2 involves a set of step-wise
commands centered around the :py:class:`tcrdist.repertoire.TCRrep` class, which
stores input data, methods, user-parameters, and results. It calls methods for
parallelized alignment and distance computation between amino acid strings in the
:py:mod:`tcrdist.pairwise` module.


This page documents the key methods in the :py:mod:`tcrdist.pairwise` module.

Source code: (`pairwise.py <https://github.com/kmayerb/tcrdist2/blob/API2/tcrdist/pairwise.py>`_)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: tcrdist.pairwise

.. autofunction:: apply_pw_distance_metric_w_multiprocessing

.. autofunction:: _f_pwdist_parallel_using_distance_wrapper

.. autofunction:: nw_metric

.. autofunction:: hm_metric

.. autofunction:: function_factory

.. autofunction:: get_pwdist_indices

.. autofunction:: get_chunked_pwdist_indices

.. autofunction:: unpack_pooled_dd_to_kkv
