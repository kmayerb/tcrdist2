
Clustering
==========

Clones may be assigned to clusters using user-defined methods. 



Hierachical Clustering
----------------------



Neigborhood Clustering
----------------------



Network Clustering
-------------------




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



Enrichment
==========