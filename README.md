# TCRdist pipeline, version 2.0.0 (DEVELOPMENT ONLY)

**tcrdist2** is new API version of TCRdist original developed by 
Phil Bradley, Jeremy Crawford, and colleagues as part of analysis of T-cell receptor specificity in
Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383). 
The original code replicating analysis in the manuscript can be found [here](https://github.com/phbradley/tcr-dist). 

## Future Installation Methods

It is highly recommended that you develop and run tcrdist2 
within a [python virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/). Doing so isolates 
the programs dependencies so installing legacy dependencies won't 
interfere with any of your other python projects. 

### tcrdist using python 2.7.11:
```bash
virtualenv venv
source ./venv/bin/activate
pip install -r requirements.txt
pip install git+https://github.com/kmayerb/--------.git@master
```

### To install the development env for tcrdist2 using python 2.7.11
```bash
virtualenv venv-dev
source ./venv-dev/bin/activate
pip install -r requirements-dev.txt
pip list
```

## Dependencies

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

The development environment can be created rapidly in using the conda package manager. 
 
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

TTODO: setup_blast.py will take care of this.

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

