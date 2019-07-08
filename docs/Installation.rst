Installation
============

tcrdist2 is a API-enabled toolkit for T-cell receptor analysis inspired by the
TCRdist pipeline developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. Nature (2017).
`doi:10.1038/nature22383 <https://www.nature.com/articles/nature22383>`_


tcrdist2 Dependencies
+++++++++++++++++++++

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


Install Development Version
+++++++++++++++++++++++++++

With python 2.7.11, pip, vritualenv already installed:

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
