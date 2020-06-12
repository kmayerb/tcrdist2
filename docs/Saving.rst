Saving Your Work
================

In contrast to the original tcrdist code, which wrote out a series of intermediate
flat files, tcrdist2 uses Pandas DataFrames and NumPy objects,
which are associated with tcrdist2 classes and held in memory.

Certain steps in the pipeline, such as the calculation of
pairwise distances and search for candidate cdr3 motifs
can take time to compute. This page documents methods for saving
computed tcrdist2 objects for later use.


Archive
#######

Archive a TCRrep instance with :py:meth:`tcrdist.repertoire.TCRrep.archive`.

.. code-block:: python
  
  tr = TCRrep(tcrdist2_df , organism = "mouse")
  tr.archive(dest = "default_archive", dest_tar_name = "default_archive.tar.gz")



Reload
######

You can recreate a TCRrep instance from an archive.tar.gz file using
:py:meth:`tcrdist.repertoire.TCRrep.rebuild`

.. code-block:: python

  tr_new = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
  tr_new.rebuild(dest_tar_name = "default_archive.tar.gz")
  

.. tip::

  :py:meth:`tcrdist.repertoire.TCRrep.archive` creates a archive which provides pairwise-matrices
  in .csv and .npy and .feather formats. For large datasets, consider manually archiving only the 
  objects of the TCRrep instance that you care about. :py:meth:`tcrdist.repertoire.TCRrep.archive` 
  and :py:meth:`tcrdist.repertoire.TCRrep.rebuild` depend on the package zipdist. 
  More information on `zipdist <https://pypi.org/project/zipdist/>`_  can be found on PyPI.



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
