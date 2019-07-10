Installation
============

tcrdist2 is a API-enabled toolkit for T Cell receptor analysis inspired by the
TCRdist pipeline developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. Nature (2017).
`doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_


Install tcrdist2 and its Dependencies
+++++++++++++++++++++++++++++++++++++
tcrdist2 can be installed using pip:

.. code-block:: none

  pip install git+https://github.com/kmayerb/tcrdist2.git@API2


However, it is highly recommended that **tcrdist2**
is installed within a `python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_.
Using a virtual environment (venv) isolates the program's dependencies so that
installing legacy packages for python (2.7.11) -- numpy (1.10.1), scipy (0.16.0),
scikit-learn (0.17.1), and matplotlib (1.4.3) --
does not interfere with any of your other ongoing python projects.

Setting up a virtual env takes less than 5 minutes using the commands below.

To configure your machine to run **tcrdist2** using the correct dependencies,
use the `*requirements.txt* <https://github.com/kmayerb/tcrdist2/blob/API2/requirements.txt>`_
file provided in the tcrdist2 github repository.

The instructions below assume that you have a working version of condas
installed or can install python 2.7.11 by other means.

.. code-block:: none

  conda create -n py27v python=2.7.11 pip virtualenv
  conda activate py27v
  virtualenv venv
  conda deactivate
  conda deactivate
  source ./venv/bin/activate
  pip install -r requirements.txt
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2


- Using condas, install a base python interpretor (Python version 2.7.11) with pip and virtualenv.
- Activate it: **conda activate py27v**
- Make a virtual env that will contain all of tcrdists dependencies: **virtualenv venv**
- Deactivate condas env (twice to deactivate py27v and base) : **conda deactivate**
- Source venv : **source ./venv/bin/activate.**
- pip install all tcrdists dependencies **pip install -r requirements.txt**
- pip install tcrdist2
- OPTIONAL: Install Blast Tools (see section below)


Or with python 2.7.11, pip, vritualenv already installed:

.. code-block:: none

  virtualenv venv
  source ./venv/bin/activate
  pip install -r requirements.txt
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2


Test the Installation
+++++++++++++++++++++

.. code-block:: none

  python -c 'import tcrdist as td; td.say_hello()'

You should see, the following:

.. code-block:: none

  > Hello: 'By recombination, random insertion, deletion and substitution,
  > the small set of genes that encode the T-cell receptor has the potential
  > to create between 10^15 and 10^20 TCR clonotypes ...
  > However, the actual diversity of a persons TCR repertoire cannot possibly
  > lie in this range. There are only an estimated 10^13 cells in the
  > human body [3]' -- Laydon et al. 2015. PMC4528489

Optional Blast Tools
++++++++++++++++++++

tcrdist2 uses `parasail <https://github.com/jeffdaily/parasail-python>'_
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
