# TCRdist pipeline, version 2.0.0

**tcrdist2** is new API version of TCRdist original developed by 
Phillip Harlan Bradley, Jeremy Chase Crawford, and colleagues as part of analysis of T-cell receptor specificity in
Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383). 
The original code replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist). 

## Installation Methods

It is highly recommended that you run or (develop) tcrdist2 
within a [python virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/). 
Using a virtual env isolates the program's dependencies so that installing legacy versions of numpy won't 
interfere with any of your other ongoing python projects.  

### Install tcrdist using python 2.7.11:
To run tcrdist using the correct dependencies,
use the *requirements.txt* file provided in the
tcrdist2 github repository.

```bash
virtualenv venv
source ./venv/bin/activate
pip install -r requirements-dev.txt
pip install pip install git+https://github.com/kmayerb/tcrdist2.git@API2
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

### Install the full dev-env using python 2.7.11
If you want to extend the functionality of tcrdist2 using the same 
development environment that we are currently using, 
configure your environment with the 
requirements-dev.txt file.

```bash
virtualenv venv-dev
source ./venv-dev/bin/activate
pip install -r requirements-dev.txt
git clone https://github.com/kmayerb/tcrdist2.git
```

## tcrdist2 is interactive!

Since tcrdist2 was designed to work with pandas DataFrames, you may find it useful to work 
interactively with ipython or jupyterlab following the provided notebook [instructions_api.ipyn]()

## Example 1: tcrdist2 on a single receptor sequence

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

pd.DataFrame(chain)
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
dependencies are available to tcrdist. Here is additional information on the tcrdist dependencies. 


tcrdist2 retains the original core dependencies
 - python (2.7.11), 
 - numpy(1.10.1),
 - scipy(0.16.0), 
 - scikit-learn(0.17.1), 
 - matplotlib(1.4.3), 

New API dependencies include:
- futures=3.2.0 
- pandas=0.20.3 
- parasail-python=1.1.16

If you prefer to use condas over pip, the development environment can be created rapidly in using the conda package manager. 
 
On a macOS or linux machine:

```bash
conda create --name tcrpy27osX python=2.7.11 scipy=0.16.0  matplotlib=1.4.3 numpy=1.10.1 futures=3.2.0 pandas=0.20.3 parasail-python=1.1.16 scikit-learn=0.17.1 jupyterlab jupyter
```
The current dependencies can be found here (tcrpy27osX.yml) and (tcrpy27linux.yml). 

## Blast Functionality and Data Files. 

BLAST 2.2.16 does not come with tcrdist2, but some of its features make use of it.

BLAST 2.2.16 can easily be downloaded here for [macosx](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-universal-macosx.tar.gz) 
or here for [linux](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-x64-linux.tar.gz) 
and put in the externals/ folder: 

tcrdist/external/blast-2.2.16/bin/
blastall and formatdb are the required executables.

TODO: setup_blast.py will take care of this.

Full data set Files from the original paper can be found 
[here](https://www.dropbox.com/s/kivfp27gbz2m2st/tcrdist_extras_v2.tgz), 
but are not included in the tcrdist2 installation.

##  Installation

tcrdist2 is not yet ready for installation.

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

