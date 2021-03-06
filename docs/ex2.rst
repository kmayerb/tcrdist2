.. toctree::
   :maxdepth: 2
   :caption: Contents:

Clustering
==========

tcrdist2 recognizes two general cluster attributes:

1. :py:attr:`TCRrep.clone_index` 

2. :py:attr:`TCRrep.clone_df`

Builtin hierarchical clustering methods :py:meth:`TCRrep.simple_clone_index` and 
:py:meth:`TCRrep.clone_index_to_df` methods are provided. 
However, clone cluster attribuets can be generated by user 
selected methods. Cluster membership need not unique. 
That is, a clone may be in multiple clusters. 

.. literalinclude:: ../tcrdist/tests/ex13.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50


.. code-block:: python

	In [2]: tr.cluster_index
	Out[2]:
	array([104, 135,  76,  64,  72,  57,  64,  71,  62,  81,  84, 104, 122,
		.......
	       141,  56, 113,  92,  94, 120,  77,  13,  56,  52,  58,   4],
	      dtype=int32)


.. code-block:: python

	In [3]: tr.cluster_df
	Out[3]:
	     cluster_id                                          neighbors  K_neighbors
	3             4  [16, 25, 26, 29, 32, 50, 61, 68, 69, 94, 103, ...           24
	91           92  [35, 38, 41, 105, 131, 146, 181, 186, 189, 206...           18
	80           81  [9, 13, 70, 74, 81, 85, 104, 106, 133, 148, 21...           17
	55           56  [18, 22, 42, 91, 98, 187, 191, 195, 217, 231, ...           15
	93           94  [15, 77, 123, 144, 173, 175, 205, 212, 259, 29...           11
	57           58  [83, 159, 203, 220, 237, 243, 249, 261, 300, 3...           11
	103         104        [0, 11, 27, 46, 65, 67, 110, 162, 170, 275]           10
	76           77                [124, 160, 185, 227, 230, 306, 318]            7
	78           79                  [78, 90, 194, 199, 251, 270, 280]            7



Sampling
========

tcrdist2 uses the pip installable package `tcrsampler <https://pypi.org/project/tcrsampler/>`_, to sample CDR3s from user-specified background. 

An example background dataset can be downloaded here: `britanova_chord_blood.csv <https://www.dropbox.com/s/rkbce72njcei4y8/?dl=1>`_.

The tcrsampler can be used to get V-gene, J-gene or join V-J-gene frequency estimates.  

.. literalinclude:: ../tcrdist/tests/ex12.py
    :linenos:
    :dedent: 1
    :language: python
    :lines: 4-50

.. code-block:: python

	In [2]: t.v_freq
	Out[2]:
	{'TRBV10-1*01': 0.001986268228829416,
		.....
	 'TRBV7-8*01': 0.01020904355302685,
	 'TRBV7-9*01': 0.03903704172001506,
	 'TRBV9*01': 0.020669023478171192}

.. code-block:: python

	In [3]: t.j_freq
	Out[3]:
	{'TRBJ1-1*01': 0.056913580882579605,
		.....
	 'TRBJ2-6*01': 0.01458019406497144,
	 'TRBJ2-7*01': 0.22612204954887244}

.. code-block:: python

	In [4]: t.vj_freq
	Out[4]: {.... 
	 ('TRBV28*01', 'TRBJ1-5*01'): 0.004184454053794057,
	 ('TRBV28*01', 'TRBJ1-6*01'): 0.002038884189985508,
	 ('TRBV28*01', 'TRBJ2-1*01'): 0.011147572827334612,
	 ('TRBV28*01', 'TRBJ2-2*01'): 0.001618765348999049}

Moreover, tcrsampler can return CDR3s based on specified V-J gene usage. 


.. code-block:: python

	In [7]: t.sample_background(v ='TRBV10-1*01', j ='TRBJ1-1*01',n=3, depth = 1, seed =1, use_frequency= True )
	   ...:
	Out[7]: ['CASSPRGDTEAFF', 'CASSEGATEAFF', 'CASSPRGDTEAFF']

