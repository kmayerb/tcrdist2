Preprocessing
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Inputs
######

The most common T cell receptor data type we anticipate supporting in
in tcrdist2:

* tcrdist - current users of TCRdist may wish to use a clones file.
* VDJDB database - T cell receptors with predicted epitope specificity
* 10X Genomics - paired-chain T cell receptor data
* Adaptive Biotechnologies - single chain bulk receptor sequencing data
* iReceptor - compilations of large T cell receptor sequencing efforts

If you are generating custom T cell receptor sequencing data and need help
using tcrdist2, please get in touch (kmayerbl@fredhutch.org).



Bulk Data Processed with MiXCR
##############################

If you have bulk sequencing data, MiXCR is an excellent option for rapidly
annotating the CDRs and identifying gene V and J gene usage. Great documentation
for MiXCR can be found at `mixcr.readthedocs.io/ <https://mixcr.readthedocs.io/en/master/>`_. 

.. code-block::

  mixcr align test_R1.fastq test_R2.fastq test2.vdjca --species hsa
  mixcr assemble test2.vdjca output.clns
  mixcr exportClones output.clns output.clns.txt
  mixcr exportAlignments test2.vdjca test2.vdjca.txt
  more test2.vdjca.txt | cut -d'\t' -f3,5,21,30 > ${res}.v_j_cdr3_n_cdr3_aa.txt


10X Genomics Data Files
#######################




Adaptive Data File
##################




VDJDB Data Files
################





TCRparse
########

Describe how to parse raw fastq data using legacy tcrdist tools.




Mappers
#######

The typical workflow for getting data into tcrdist2 is to load a data file
with Pandas and then apply a mapper to select and rename columns to match
those required by tcrdist2. see :ref:`Mappers`


We have begun writing wrappers for various data sources. Mappers are
just python dictionaries and can be easily emulated if you data headers
are unique.

.. code-block:: python

  import pandas as pd
  from tcrdist import mappers
  mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping
  tcrdist1_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df_PA,
                                              mapping = mapping)


Required Fields
###############

The follow are required fields for tcrdist2 modules:
