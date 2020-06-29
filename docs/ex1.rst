.. toctree::
   :maxdepth: 3
   :caption: Contents:

**tcrdist2 provides flexible distance measures for comparing T cell receptors across complementarity determining regions (CDRs)**

The following short examples use the :py:class:`TCRrep` class and :py:meth:`tcrdist2` method, which is appropriate for small and medium-sized datasets of 10-10,000 TCR clones. 
The data used in these examples can be downloaded here: `dash.csv <https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/tcrdist/test_files_compact/dash.csv>`_ (375KB). The :py:meth:`tcrdist2`
method automates some steps that can be seen in the longer examples.


TCR Distances
=============

Hamming Distance
----------------

Hamming distance between receptors, using the CDR1, CDR2, CDR2.5 and CDR3.
By default, the CDRs are aligned with the Needleman–Wunsch
algorithm using BLOSUM62 scoring matrix. 

.. literalinclude:: ../tcrdist/tests/ex1.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50



.. code-block:: python

	In [2]: tr.pw_tcrdist
	Out[2]:
	array([[ 0, 45, 41, ..., 42, 37, 26],
	       [45,  0, 49, ..., 46, 47, 45],
	       [41, 49,  0, ..., 29, 28, 45],
	       ...,
	       [42, 46, 29, ...,  0, 28, 27],
	       [37, 47, 28, ..., 28,  0, 39],
	       [26, 45, 45, ..., 27, 39,  0]], dtype=int16)


Reciprocal Alignment Distance
-----------------------------

Reciprocal allignment distance between receptors, 
using the CDR1, CDR2, CDR2.5 and CDR3. 
By default, the CDRs are aligned with the Needleman–Wunsch
algorithm using BLOSUM62 scoring matrix. 



.. tip::

	The reciprocal distance is calucated using the Needlman-Wunsch algorith and a scoring matrix, such that
	the distance is score(x,x) + score(y,y) - 2 * score(x,y)   

.. literalinclude:: ../tcrdist/tests/ex2.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50

.. code-block:: python

	In [2]: tr.pw_tcrdist
	Out[2]:
	array([[  0, 504, 449, ..., 440, 407, 280],
	       [504,   0, 573, ..., 506, 533, 518],
	       [449, 573,   0, ..., 313, 336, 479],
	       ...,
	       [440, 506, 313, ...,   0, 319, 268],
	       [407, 533, 336, ..., 319,   0, 455],
	       [280, 518, 479, ..., 268, 455,   0]], dtype=int16)


Chain-Specific Distance
-----------------------

Chain specific distance attributes are provided as :py:attr:`TCRrep.pw_alpha` and 
:py:attr:`TCRrep.pw_beta`.


.. literalinclude:: ../tcrdist/tests/ex3.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


.. code-block:: python

	In [2]: tr.pw_tcrdist
	Out[2]:
	array([[ 0, 45, 41, ..., 42, 37, 26],
	       [45,  0, 49, ..., 46, 47, 45],
	       [41, 49,  0, ..., 29, 28, 45],
	       ...,
	       [42, 46, 29, ...,  0, 28, 27],
	       [37, 47, 28, ..., 28,  0, 39],
	       [26, 45, 45, ..., 27, 39,  0]], dtype=int16)

	In [3]: tr.pw_alpha
	Out[3]:
	array([[ 0, 25, 23, ..., 27, 22, 23],
	       [25,  0, 28, ..., 27, 27, 27],
	       [23, 28,  0, ..., 24, 23, 28],
	       ...,
	       [27, 27, 24, ...,  0, 26, 14],
	       [22, 27, 23, ..., 26,  0, 25],
	       [23, 27, 28, ..., 14, 25,  0]], dtype=int16)

	In [4]: tr.pw_beta
	Out[4]:
	array([[ 0, 20, 18, ..., 15, 15,  3],
	       [20,  0, 21, ..., 19, 20, 18],
	       [18, 21,  0, ...,  5,  5, 17],
	       ...,
	       [15, 19,  5, ...,  0,  2, 13],
	       [15, 20,  5, ...,  2,  0, 14],
	       [ 3, 18, 17, ..., 13, 14,  0]], dtype=int16)


CDR-Specific Distance
---------------------

CDR-specific tcr-distances.


.. literalinclude:: ../tcrdist/tests/ex4.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


.. tip::

  CDR2.5 alpha and CDR2.5 beta, the pMHC-facing loop between CDR2
  and CDR3, are referred to in tcrdist2 as pmhc_a and phmc_b, respectively.


.. code-block:: python

	In [2]: tr.cdr3_b_aa_pw
	Out[2]:
	array([[0, 7, 5, ..., 2, 2, 2],
	       ...,
	       [2, 6, 5, ..., 1, 2, 0]], dtype=int16)

	In [3]: tr.pmhc_b_aa_pw
	Out[3]:
	array([[0, 4, 4, ..., 4, 4, 0],
	       ...,
	       [0, 4, 4, ..., 4, 4, 0]], dtype=int16)

	In [4]: tr.cdr2_b_aa_pw
	Out[4]:
	array([[0, 6, 5, ..., 5, 5, 0],
	       ...,
	       [0, 6, 5, ..., 5, 5, 0]], dtype=int16)

	In [5]: tr.cdr1_b_aa_pw
	Out[5]:
	array([[0, 3, 4, ..., 4, 4, 1],
	       ...,
	       [1, 2, 3, ..., 3, 3, 0]], dtype=int16)



Custom CDR Weights
------------------

Customize weight of CDR-specific contributions to the overall :py:attr:`.pw_tcrdist` matrix. 
For instance, here a weight of 3 is applied to the CDR3 compared to other CDRs.

.. literalinclude:: ../tcrdist/tests/ex4.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50

.. tip::

  Check weightings used attribute :py:attr:`TCRrep.tr.paired_tcrdist_weights`. New weightings 
  can be applied wiht replacment weights or by matrix addition. 


Custom Distances
================


pwseqdist
---------
tcrdist2 uses the pip installable package `pwseqdist <https://pypi.org/project/pwseqdist/>`_, to compute distances using parrallel processing. 

In addition to the methods above, almost any metric comparing strings can be used. 
For instance, the package python-Levenshtein calculates the edit distance between 
strings which has become a popular way of comparing large sets of CDR3s.

.. literalinclude:: ../tcrdist/tests/ex6.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


.. code-block:: python

	In [2]: dmat
	Out[2]:
	array([[0, 1, 2],
	       [1, 0, 1],
	       [2, 1, 0]])


Levenshtein Distance
--------------------

pwseqdist functionality can be applied to any string column in the TCRrep clone DataFrame with :py:meth:`.custom_dmat`.

.. literalinclude:: ../tcrdist/tests/ex7.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50

.. code-block:: python

	In [2]: tr.custom_dmat(cdr = 'cdr3_b_aa', metric =Levenshtein.distance, processes = 1)
	Out[2]:
	array([[0, 7, 5, ..., 2, 2, 2],
	       [7, 0, 8, ..., 6, 7, 6],
	       [5, 8, 0, ..., 5, 5, 5],
	       ...,
	       [2, 6, 5, ..., 0, 2, 1],
	       [2, 7, 5, ..., 2, 0, 2],
	       [2, 6, 5, ..., 1, 2, 0]])


.. tip::

	Using :py:meth:`.add_custom_dmat` will automatically create a new attribute with the name `custom_[selected_cdr]_pw`. For instance,


.. code-block:: python

	In [2]: tr.add_custom_dmat(cdr = 'cdr3_b_aa', metric =Levenshtein.distance, processes = 1)
	
	In [3]: tr.custom_cdr3_b_aa_pw
	Out[3]:
	array([[0, 7, 5, ..., 2, 2, 2],
	       [7, 0, 8, ..., 6, 7, 6],
	       [5, 8, 0, ..., 5, 5, 5],
	       ...,
	       [2, 6, 5, ..., 0, 2, 1],
	       [2, 7, 5, ..., 2, 0, 2],
	       [2, 6, 5, ..., 1, 2, 0]])



Levenshtein Distance on CDRs
----------------------------

With, :py:meth:`.add_custom_dmat` custom metrics can be combined to compute distances across multipe CDRs. 
For instance,

.. literalinclude:: ../tcrdist/tests/ex8.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


.. code-block:: python

	In [2]: tr.pw_levenshtein_tcrdist_beta
	Out[2]:
	array([[ 0, 20, 18, ..., 15, 15,  3],
	       [20,  0, 21, ..., 19, 20, 18],
	       [18, 21,  0, ...,  5,  5, 17],
	       ...,
	       [15, 19,  5, ...,  0,  2, 13],
	       [15, 20,  5, ...,  2,  0, 14],
	       [ 3, 18, 17, ..., 13, 14,  0]])


Saving
======

tcrdist2 uses the pip installable package `zipdist <https://pypi.org/project/zipdist/>`_, to compute to archive and reload Numpy and Pandas TCRrep attributes. 

An archive can be created at the time :py:meth:`tcrdist2` is called by setting :py:attr:`tcrdist2.save` to True.

.. literalinclude:: ../tcrdist/tests/ex9.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


An archive can be created at any stage using  :py:meth:`TCRrep.archive`. All JSON serializable, Numpy array, and Pandas DataFrame attributes are saved.

.. literalinclude:: ../tcrdist/tests/ex10.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


Reloading
=========

A :py:class:`TCRrep` instance can be rebuilt from the .tar.gz archive file. For instance:

.. literalinclude:: ../tcrdist/tests/ex11.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50

All attributes can be accessed as before.

.. code-block:: python

	In [2]: tr.pw_tcrdist
	Out[2]:
	array([[ 0, 45, 41, ..., 42, 37, 26],
	       [45,  0, 49, ..., 46, 47, 45],
	       [41, 49,  0, ..., 29, 28, 45],
	       ...,
	       [42, 46, 29, ...,  0, 28, 27],
	       [37, 47, 28, ..., 28,  0, 39],
	       [26, 45, 45, ..., 27, 39,  0]], dtype=int16)

