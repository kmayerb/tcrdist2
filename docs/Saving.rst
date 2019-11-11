Saving Your Work
================

In contrast to the original tcrdist code, which wrote out a series of intermediate
flat files, tcrdist2 uses Pandas DataFrames and numpy objects,
which are associated with tcrdist2 classes and held in memory.

Certain steps in the pipeline, such as the calculation of
pairwise distances and search for candidate cdr3 motifs
can take time to compute. This page documents methods for saving
computed tcrdist2 objects for later use.

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

Most tcrdistances can be expressed as integers without loss of information.
The `int16` data type, native to numpy, is a choice that
can reduce the file size of a saved TCRrep instance.

To convert numpy arrays with float64 data to arrays with
int16 data storage, use
:py:meth:`tcrdist.repertoire.TCRrep.reduce_file_size`:


.. code-block:: python

  tr.reduce_file_size(data_type = "int16")
  tr.save_as_hdf5("TCRrep_file.h5")



Details
#######

This section provides more details on methods for saving tcrdist2 objects for later
use.  Rather than write each element to a separate file,
an instance of the :py:class:`tcrdist.repertoire.TCRrep` class can be
serialized in it entirety. The preferred method for saving a TCRrep instance
is to save it's attributes to an HDF5 file.

.. tip ::

  If you are unfamiliar with HDF5 files,
  consider downloading the free `HDF5 viewer <https://www.hdfgroup.org/downloads/hdfview/>`_.
  It enables inspection of the HDF5 file's contents.

An example shows how it all works:

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


To start where you left off:

* 1. initialize a new TCRrep instance  (i.e., tr2 in the code below), providing an empty pd.DataFrame for the cell_df argument. Specify the correct
* 2. Call :py:meth:`tcrdist.repertoire.TCRrep.rebuild_from_hdf5`

.. code-block:: python

  tr2 = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
  tr2.rebuild_from_hdf5("TCRrep_file.h5")


You can inspect and confirm that the crucial attributes
are identical in tr and tr2, with the following code:

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


 Flat Files
 ##########

 Individual elements stored in the
 TCRrep class can be saved to text files directly. (For DataFrame, see
 `pandas.DataFrame.to_csv() <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html#pandas-dataframe-to-csv>`_
 and for numpy arrays, see: `numpy.savetxt() <https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html>`_)


Pickle
######

.. tip::

 READ: TCRrep instances and their contents can be pickled. That's good.
 **But pickling is cursed!** That's bad. In fact, we advise against
 using pickle for long-term storage of complex objects.
 This is because future versions of tcrdist2 may not
 recognize pickled files made from a prior version!!! But the pickle comes
 with your choice of toppings. That's good. But the toppings are also
 `cursed <https://youtu.be/Krbl911ZPBA>`_

If you wish to pickle a TCRrep instance (*caveat emptor*):

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
  calculate more tcrdistances with this TCRrep instance,
  it is also necessary to call
  :py:meth:`tcrdist.repertoire.TCRrep._initialize_chain_specific_attributes` which restores
  parasail distance matrices which can not be pickled.
