Storage
=======

.. toctree::
   :maxdepth: 2
   :caption: Contents:

This page documents the class :py:class:`tcrdist.storage.StoreIO` within the
:py:mod:`tcrdist.storage` module
`storage.py <https://github.com/kmayerb/tcrdist2/blob/API2/tcrdist/storage.py>`_ (source code).
We describe it here in detail as it is used extensively to pass inputs and outputs
between functions in :py:class:`tcrdist.cdr3_motif.TCRMotif` and
:py:class:`tcrdist.subset.TCRsubset`. 

In the analysis of CDR3 motifs, tcrdist2 utilizes the algorithms
developed in tcdist1. These algorithms have been modularized. Input and Output
is facilitated :py:class:`tcrdist.storage.StoreIOMotif` and
:py:class:`tcrdist.storage.StoreIOEntropy`, which inherit from :py:class:`tcrdist.storage.StoreIO`.


Example usage of StoreIO
########################

Initialize
**********

.. code-block:: python

  >>> from tcrdist.storage import StoreIO
  >>> from collections import namedtuple
  >>> S = StoreIO()
  >>> print(S)
  StoreIO Attributes:
	name : <class 'str'> : 0x1c1b52cb20 : StoreIO
	valid_attrs : <class 'list'> : 0x1c62198148 : ['a', 'b', 'c', 'd']
	valid_attrs_type : <class 'list'> : 0x1c3f066b08 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
	a : <class 'NoneType'> : 0x104d6df78 : None
	b : <class 'NoneType'> : 0x104d6df78 : None
	c : <class 'NoneType'> : 0x104d6df78 : None
	d : <class 'NoneType'> : 0x104d6df78 : None

Notice that there are four valid attributes, and the valid types are specified
as integers.


Initialize with Keyword Args
****************************

.. code-block:: python

  >>> from tcrdist.storage import StoreIO
  >>> from collections import namedtuple
  >>> abcd = namedtuple('abcd', ['a','b','c','d'])
  >>> x = abcd('1','2','3','4')
  >>> S = StoreIO(**x._asdict())
  >>> print(S)
  StoreIO Attributes:
  	name : <class 'str'> : 0x1046a01f0 : StoreIO
  	valid_attrs : <class 'list'> : 0x1c19a93288 : ['a', 'b', 'c', 'd']
  	valid_attrs_type : <class 'list'> : 0x1c197423c8 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
  	a : <class 'str'> : 0x102556c00 : 1
  	b : <class 'str'> : 0x102556c70 : 2
  	c : <class 'str'> : 0x102556fb8 : 3
  	d : <class 'str'> : 0x10255e030 : 4

Notice that input types are strings (i.e., '1', '2', '3', '4'),
but the valid types are specified as integers. StoreIO permits a
check, which is useful before passing the StoreIO instance
to a function.

Type Validation and Coercion
****************************

If `_validate_attrs()` is called, a TypeError results.

.. code-block:: python

  >>> S._validate_attrs()
  TypeError: StoreIO.a must be of type: <class 'int'>

If `_coerce_attrs()` is called, type coercion is attempted. These strings can
be coerced to integer, so the type validation passes.

.. code-block:: python

  >>> S._coerce_attrs()
  >>> print(S)
  StoreIO Attributes:
  	name : <class 'str'> : 0x1046a01f0 : StoreIO
  	valid_attrs : <class 'list'> : 0x1c19a93288 : ['a', 'b', 'c', 'd']
  	valid_attrs_type : <class 'list'> : 0x1c197423c8 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
  	a : <class 'int'> : 0x10217d5b0 : 1
  	b : <class 'int'> : 0x10217d5d0 : 2
  	c : <class 'int'> : 0x10217d5f0 : 3
  	d : <class 'int'> : 0x10217d610 : 4
  >>> S._validate_attrs()
  True

Use with Functions
******************

.. code-block:: python

  >>> def ten_x_ab(StoreIOinstance):
  ...   """
  ...   10X the values of attribute a and attribute b
  ...   """
  ...   StoreIOinstance.a = StoreIOinstance.a * 10
  ...   StoreIOinstance.b = StoreIOinstance.b * 10
  ...   StoreIOinstance._validate_attrs()
  ...   return(StoreIOinstance)
  >>> from tcrdist.storage import StoreIO
  >>> from collections import namedtuple
  >>> S = StoreIO()
  >>> abcd = namedtuple('abcd', ['a','b','c','d'])
  >>> x = abcd('1','2','3','4')
  >>> S = StoreIO(**x._asdict())
  >>> S._coerce_attrs()
  >>> S._validate_attrs()
  >>> ten_x_ab(S)
  >>> print(S)
  StoreIO Attributes:
	name : <class 'str'> : 0x1046a01f0 : StoreIO
	valid_attrs : <class 'list'> : 0x1c19a90648 : ['a', 'b', 'c', 'd']
	valid_attrs_type : <class 'list'> : 0x1c19a36588 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
	a : <class 'int'> : 0x10217d6d0 : 10
	b : <class 'int'> : 0x10217d810 : 20
	c : <class 'int'> : 0x10217d5f0 : 3
	d : <class 'int'> : 0x10217d610 : 4

Set Attributes with Keyword Args
********************************

StoreIO attributes can be set after initialization with
:py:func:`tcrdist.storage.StoreIO.set_attrs_with_kwargs`

Notice that this method only assigns valid attributes even if more keyword
args are provided (e.g., e = 5 in the example below).
Additional attributes can be assigned to the
StoreIO instance if `validate_kwargs = False`, as shown
below with the keyword 'e' which was not pre-specified in
`StoreIO.valid_attrs`

.. code-block:: python

  >>> from tcrdist.storage import StoreIO
  >>> from collections import namedtuple
  >>> S = StoreIO()
  >>> print(S)
    StoreIO Attributes:
  	name : <class 'str'> : 0x103c94340 : StoreIO
  	valid_attrs : <class 'list'> : 0x1c19085b88 : ['a', 'b', 'c', 'd']
  	valid_attrs_type : <class 'list'> : 0x1c19085bc8 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
  	a : <class 'NoneType'> : 0x10165ef78 : None
  	b : <class 'NoneType'> : 0x10165ef78 : None
  	c : <class 'NoneType'> : 0x10165ef78 : None
  	d : <class 'NoneType'> : 0x10165ef78 : None
  >>> S.set_attrs_with_kwargs(a = 5, b = 8, c = 2, d = 4, e = 5)
  >>> print(S)
  StoreIO Attributes:
  	name : <class 'str'> : 0x103c94340 : StoreIO
  	valid_attrs : <class 'list'> : 0x1c19085b88 : ['a', 'b', 'c', 'd']
  	valid_attrs_type : <class 'list'> : 0x1c19085bc8 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
  	a : <class 'int'> : 0x1016a3630 : 5
  	b : <class 'int'> : 0x1016a3690 : 8
  	c : <class 'int'> : 0x1016a35d0 : 2
  	d : <class 'int'> : 0x1016a3610 : 4
  >>> S.set_attrs_with_kwargs(validate_kwargs = False, a = 5, b = 8, c = 2, d = 4, e = 5)
  >>> print(S)
  StoreIO Attributes:
  	name : <class 'str'> : 0x103c94340 : StoreIO
  	valid_attrs : <class 'list'> : 0x1c19085b88 : ['a', 'b', 'c', 'd']
  	valid_attrs_type : <class 'list'> : 0x1c19085bc8 : [<class 'int'>, <class 'int'>, <class 'int'>, <class 'int'>]
  	a : <class 'int'> : 0x1016a3630 : 5
  	b : <class 'int'> : 0x1016a3690 : 8
  	c : <class 'int'> : 0x1016a35d0 : 2
  	d : <class 'int'> : 0x1016a3610 : 4
  	e : <class 'int'> : 0x1016a3630 : 5



StoreIO
#######

.. automodule:: tcrdist.storage

.. autoclass:: StoreIO
    :members: __init__, set_attrs_with_kwargs, _validate_attrs, _coerce_attrs, _type_check, _type_coerce

StoreIOMotif
############

.. autoclass:: StoreIOMotif
    :show-inheritance:

StoreIOEntropy
##############

.. autoclass:: StoreIOEntropy
    :show-inheritance:
