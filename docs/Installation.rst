Installation
============

`tcrdist2 <https://github.com/kmayerb/tcrdist2>`_ is a
python API-enabled toolkit expanding on the T cell receptor analysis pipeline
developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. Nature (2017). `doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_

tcrdist2 is designed to be run using Python 3. The latest version can
be installed from GitHub using pip:

.. code-block:: none

  pip install git+https://github.com/kmayerb/tcrdist2.git@API2

Before installing tcrdist2, check that you have the appropriate dependencies.

Dependencies
++++++++++++

Tcrdist2 requires wheel, numpy, scipy, matplotlib, scipy, scikit-learn, tables,
as well as the C-based sequence aligner parasail.

Recommended python 3.6 dependencies are specified in the
`requirements36.txt <https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/requirements36.txt>`_ .
It should be placed the directory where you launch the installation.
The following will install tcrdist on a linux or macOSX within a virtual
environment (Recommended).

.. code-block:: none
  
  pip install wheel
  pip install -r requirements36.txt
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2


Virtual Env
+++++++++++

You may want to use a virtual env to isolate tcrdist within its own environment.
If you are creating a virtual environment on a Linux macine for the first time, the
following preparatory steps may be required.

.. code-block:: none

  apt-get python3-dev
  pip install virtualenv

With virtualenv installed, continue with the installation:

.. code-block:: none

  python3 -m venv ./venv
  pip install wheel
  source venv/bin/activate
  pip install -r requirements36.txt
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2


Test the Installation
+++++++++++++++++++++

.. code-block:: none

  python -c 'import tcrdist as td; td.say_hello()'

You should see, the following:

.. code-block:: none

  Hello: 'By recombination, random insertion, deletion and substitution,
  the small set of genes that encode the T-cell receptor has the potential
  to create between 10^15 and 10^20 TCR clonotypes ...
  However, the actual diversity of a persons TCR repertoire cannot possibly
  lie in this range. There are only an estimated 10^13 cells in the
  human body [3]' -- Laydon et al. 2015. PMC4528489


Optional Installation Files
+++++++++++++++++++++++++++

To install important next-generation reference files:

.. code-block:: none

  python -c "import tcrdist as td; td.setup_db.install_all_next_gen()"

To install blast within your virtual environment
(see Optional Blast Tools for installation on other platforms):

.. code-block:: none

  python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'ncbi_linux')"


Optional Blast Tools
++++++++++++++++++++

tcrdist2 uses `parasail <https://github.com/jeffdaily/parasail-python>`_
for sequence alignments; however, some features have the option to use BLAST instead.

The BLAST version 2.2.16 used in Dash et al. 2017, can be optionally installed with
the followings commands.

After installing tcrdist2, if working in a macOSX environment:

.. code-block:: none

  python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'ncbi_osx');"


After installing tcrdist2, if working in a Linux environment:

.. code-block:: none

  python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'ncbi_linux');"


If the NCBI links change, a backup download link can be accessed by changing the *download_from* argument:

.. code-block:: none

  python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_osx');"


.. code-block:: none

  python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_linux);"
