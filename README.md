[![Build Status](https://travis-ci.com/kmayerb/tcrdist2.svg?branch=API2)](https://travis-ci.com/kmayerb/tcrdist2)

# TCRdist pipeline, version 2.0.0

2019-06-17

**tcrdist2** is a API-enabled version of the TCRdist toolkit developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and colleagues as part of a T-cell receptor epitope specificity analysis in Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383). The original code for replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist).

## Recommended Installation Method

tcrdist2 can be installed using pip:

```bash
pip install git+https://github.com/kmayerb/tcrdist2.git@API2
```

However, it is highly recommended that **tcrdist2**
is installed within a [python virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/). Using a virtual environment (venv) isolates the program's dependencies so that installing legacy packages for
python (2.7.11) -- numpy (1.10.1), scipy (0.16.0), scikit-learn (0.17.1), and matplotlib (1.4.3) -- does not interfere with any of your other ongoing python projects. Setting up a virtual env takes less than 5 minutes using the commands below.

To configure your machine to run **tcrdist2** using the correct dependencies, use the [*requirements.txt*](https://github.com/kmayerb/tcrdist2/blob/API2/requirements.txt) file provided in the tcrdist2 github repository. The instructions below assume that you have a working version of condas installed or can install python 2.7.11 by other means.

```bash
conda create -n py27v python=2.7.11 pip virtualenv
conda activate py27v
virtualenv venv
conda deactivate
conda deactivate
source ./venv/bin/activate
pip install -r requirements.txt
pip install git+https://github.com/kmayerb/tcrdist2.git@API2
```

- Using condas, install a base python interpretor (Python version 2.7.11) with pip and virtualenv.
- Activate it: **conda activate py27v**
- Make a virtual env that will contain all of tcrdists dependencies: **virtualenv venv**
- Deactivate condas env (twice to deactivate py27v and base) : **conda deactivate**
- Source venv : **source ./venv/bin/activate.**
- pip install all tcrdists dependencies **pip install -r requirements.txt**
- pip install tcrdist2
- OPTIONAL: Install Blast Tools (see section below)


With python 2.7.11, pip, vritualenv already installed:

```bash
virtualenv venv
source ./venv/bin/activate
pip install -r requirements.txt
pip install git+https://github.com/kmayerb/tcrdist2.git@API2
```

### Test the installation
```bash
python -c 'import tcrdist as td; td.say_hello()'
```
You should see, the following:
```
> Hello: 'By recombination, random insertion, deletion and substitution,
> the small set of genes that encode the T-cell receptor has the potential
> to create between 10^15 and 10^20 TCR clonotypes ...
> However, the actual diversity of a persons TCR repertoire cannot possibly
> lie in this range. There are only an estimated 10^13 cells in the
> human body [3]' -- Laydon et al. 2015. PMC4528489
```

### Optional: Install Blast Tools

tcrdist2 uses [parasail](https://github.com/jeffdaily/parasail-python)
for sequence alignments; however, some features have the option to use BLAST instead.

The BLAST version 2.2.16 used in Dash et al. 2017, can be optionally installed with
the followings commands.

After installing tcrdist2, if working in a macOSX environment:
```bash
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'ncbi_osx');"
```

After installing tcrdist2, if working in a Linux environment:
```bash
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'ncbi_linux');"
```

If the NCBI links change, a backup download link can be accessed by changing the *download_from* argument:

```bash
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_osx');"
```

```bash
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_osx');"
```

### Configure the full dev-env using python 2.7.11
If you want to test or extend the functionality of **tcrdist2** using the same
development environment that we are currently using,
configure your environment with the [requirements-dev.txt](https://github.com/kmayerb/tcrdist2/blob/API2/requirements-dev.txt) file.

```bash
virtualenv venv-dev
source ./venv-dev/bin/activate
pip install -r requirements-dev.txt
git clone https://github.com/kmayerb/tcrdist2.git
```

## Work with tcrdist2 interactively!

**tcrdist2** was works with [Pandas](https://vimeo.com/59324550) DataFrames.
Therefore, you may find it useful to work interactively with ipython.

We are working on providing an ipython notebook with example instructions
[instructions_api.ipyn](https://github.com/kmayerb/tcrdist2/blob/API2/instructions_api.ipynb)

## Example 1: tcrdist2 on a single receptor sequence

```ipython
betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
chain = td.processing.processNT(organism = 'human', chain = 'B', nuc = betaNT, quals = betaQuals, use_parasail = True)
pd.DataFrame(chain)
```

```ipython
In [1]: import tcrdist as td

In [2]: import pandas as pd

In [3]: betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGG
   ...: TCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCT
   ...: TCAGAGTGGAGTCANN'
   ...: betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20
   ...: .25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47
   ...: .53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68
   ...: .68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68
   ...: .68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59
   ...: .59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15
   ...: .21.20.25.3.0.0'

In [4]: chain = td.processing.processNT(organism = 'human',
   ...:                                 chain = 'B',
   ...:                                 nuc = betaNT,
   ...:                                 quals = betaQuals,
   ...:                                 use_parasail = True)

In [5]: pd.DataFrame(chain)
```

## Example 2: tcrdist2 on a batch of sequences

### readPairedSequences

```python
psDf = td.processing.readPairedSequences(paired_seqs_file = "tcrdist/datasets/test_human_pairseqs.tsv",
                                         organism = "human",
                                         use_parasail = True);
```
### computeProbs

```python
probDf = td.processing.computeProbs(psDf)
psDf = psDf.join(probDf)
```

###  identifyClones
```python
clonesDf = td.processing.identifyClones(psDf)                                         
```



## More Information on Dependencies

Following the instructions above and setting up a virtual environment should take care of ensuring the proper
dependencies are available to tcrdist2.  In addition to original **TCRdist** dependencies, **tcrdist2**
requires the following new dependencies:
- futures=3.2.0
- pandas=0.20.3
- parasail-python=1.1.16


If you are familiar with dependency management in condas,
a **tcrdist2** development environment can be created rapidly using conda:

```bash
conda create --name tcrpy27osX python=2.7.11 scipy=0.16.0  matplotlib=1.4.3 numpy=1.10.1 futures=3.2.0 pandas=0.20.3 parasail-python=1.1.16 scikit-learn=0.17.1 jupyterlab jupyter
```
The current dependencies can be found here [tcrpy27osX.yml](https://github.com/kmayerb/tcrdist2/blob/API2/tcrpy27osX.yml), so
recreating the conda env is a one-liner:

```bash
conda env create -n tcrpy27osX -f tcrpy27osX.yml
```
## Blast Functionality and Data Files.
.

BLAST 2.2.16 can easily be downloaded here for [macosx](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-universal-macosx.tar.gz)
or here for [linux](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-x64-linux.tar.gz)
and put in the externals/ folder if the installation commands in the above section are not successful.

folder to place BLAST execuatable: tcrdist/external/blast-2.2.16/bin/
- blastall and formatdb are the only required executables.


The data set and files from the original Dash et al. 2017 paper can be found
[here](https://www.dropbox.com/s/kivfp27gbz2m2st/tcrdist_extras_v2.tgz),
but are not included in the tcrdist2 installation.



---
# Citing

Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas

Nature (2017) doi:10.1038/nature22383

---
# Thanks

Part of this analysis uses parameters derived from nextgen data in publicly released studies. We are grateful to the authors of those studies for making their data available. See the `README_db.txt` file in `./db/` (after running `setup.py`)

The code uses the command line parsing toolkit `blargs`. See the license and info in `external/blargs/`

The tables in the .html results can be sorted thanks to `tablesorter`. See the license and info in `external/tablesorter/`

Sequence parsing relies on the BLAST suite, see info in `external/blast-2.2.16/`

---
# UPDATES

tcrdist2 can be installed, but it is worth noting that tcrdist2 is undergoing heavy development in June 2019.

We anticipate a new release with expanded functionality and vizualization tools in July 2019.
