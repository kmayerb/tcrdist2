Saving Your Work
================

Save
####

Save with :py:meth:`tcrdist.repertoire.TCRrep.save_as_hdf5`:

.. code-block:: python

  tr.save_as_hdf5("TCRrep_file.h5")


Reload
######

Reload from hdf5 with :py:meth:`tcrdist.repertoire.TCRrep.rebuild_from_hdf5`

.. code-block:: python

  tr2 = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
  tr2.rebuild_from_hdf5("TCRrep_file.h5")


Reduce File Size
################

Most tcrdistance can be expressed as integers without loss of information.
Th int16 data types native to numpy are a great choice and
can reduce dramatically the file size of a saved
TCRrep instance. We've made that easy. If you want to convert numpy arrays
from float64 type to int16 type (or any other type) storage use
:py:meth:`tcrdist.repertoire.TCRrep.reduce_file_size`:


.. code-block:: python

  tr.reduce_file_size(data_type = "int16")
  tr.save_as_hdf5("TCRrep_file.h5")



Details
#######

Suppose you have done some work with a tcrdist2 objects and you want to save
them for further use. Computation of pairwise distances can take
a few minutes as the number of TCR clones increases.

Individual elements of the TCRrep class can be saved to
text files. See
`pandas.DataFrame.to_csv <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html#pandas-dataframe-to-csv>`_
for dataframes and, for
numpy arrays, see `numpy.savetxt <https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html>`_

However, a :py:class:`tcrdist.repertoire.TCRrep` class instance can be
serialized in it entirety. The preferred method for saving a TCRrep instance is to serialize the crucial
attributes to an HDF5 file

.. tip ::

  If you are unfamiliar with HDF5 files,
  consider downloading the free `HDF5 viewer <https://www.hdfgroup.org/downloads/hdfview/>`_
  It enables inspection of the HDF5 file's contents.

Here is how it all works:

.. code-block:: python

  import pandas as pd
  from tcrdist.repertoire import TCRrep
  from tcrdist import mappers

  tcrdist_clone_fn = 'tcrdist/test_files/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
  tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")
  ind = (tcrdist_clone_df.epitope == "PA") | (tcrdist_clone_df.epitope == "F2")
  tcrdist_clone_df = tcrdist_clone_df[ind].copy()
  mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping
  tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,
                                              mapping = mapping)
  tr = TCRrep(tcrdist2_df , organism = "mouse")
  tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
  tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)
  tr.index_cols = ['clone_id', 'subject', 'epitope',
                 'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                 'cdr3_a_aa', 'cdr3_b_aa']
  tr.deduplicate()
  tr._tcrdist_legacy_method_alpha_beta()
  tr.stored_tcrdist = None
  tr.reduce_file_size()
  tr.save_as_hdf5("TCRrep_file.h5")


When returning, initialize a new TCRrep instance (i.e., tr2 in the code below).
Provide a empty pd.DataFrame for the cell_df and specify the correct
model organism. Then call :py:meth:`tcrdist.repertoire.TCRrep.rebuild_from_hdf5`

.. code-block:: python

  tr2 = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
  tr2.rebuild_from_hdf5("TCRrep_file.h5")


You can inspect the new object. Here we confirm that the crucial attributes
are identical in tr and tr2.

.. code-block:: python

  >>> import numpy as np
  >>> {x : np.all(getattr(tr, x) == getattr(tr2, x)) for x in tr.__dict__.keys()}
  {'cell_df': True,
 'chains': True,
 'organism': True,
 'pwdist_df': True,
 'clone_df': True,
 'index_cols': True,
 'stored_tcrdist': False,
 'paired_tcrdist': True,
 'paired_tcrdist_weights': True,
 'meta_cols': True,
 'project_id': True,
 'all_genes': False,
 'imgt_aligned_status': True,
 'cdr3_a_aa_smat': True,
 'cdr2_a_aa_smat': True,
 'cdr1_a_aa_smat': True,
 'pmhc_a_aa_smat': True,
 'cdr3_b_aa_smat': True,
 'cdr2_b_aa_smat': True,
 'cdr1_b_aa_smat': True,
 'pmhc_b_aa_smat': True,
 'cdr3_a_aa_pw': True,
 'cdr3_b_aa_pw': True,
 'dist_a': True,
 'dist_b': True}


Pickle
######

.. tip::

 READ: TCRrep instances and all there contents can be pickled. That's good.
 **But pickling is cursed!** That's bad. In fact, we advise against
 using pickle for long-term storage as upgraded versions of tcrdist2 may not
 recognize pickled file made from a prior version!!! But the pickle comes
 with your choice of toppings. That's good. But the toppings are also
 `cursed <https://www.youtube.com/watch?v=KaWdKSR2rtA>`_

If you wish to pickle a TCRrep instance:

.. code-block:: python

  tr._pickle("TCRrep_file.p")

To get it back:

.. code-block:: python

  import pickle
  tr3 = pickle.load(open("TCRrep_file.p", "rb"))
  tr3._initialize_chain_specific_attributes()
  {x : np.all(getattr(tr, x) == getattr(tr3, x)) for x in tr.__dict__.keys()}

.. tip::

  If you intend to
  calcuate more tcrdistances, it is also necessary to call
  :py:meth:`tcrdist.repertoire.TCRrep._initialize_chain_specific_attributes` which restores
  parasail distance matrices which can not be serialized.
