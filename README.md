[![Build Status](https://travis-ci.com/kmayerb/tcrdist2.svg?branch=API2)](https://travis-ci.com/kmayerb/tcrdist2)
[![Build Docs Status](https://readthedocs.org/projects/tcrdist2/badge/?version=latest)](https://tcrdist2.readthedocs.io/en/latest/)
[![Coverage Status](https://coveralls.io/repos/github/kmayerb/tcrdist2/badge.svg?branch=API2)](https://coveralls.io/github/kmayerb/tcrdist2?branch=API2)

# tcrdist2

#### tcrdist2 provides flexible distance measures for comparing T cell receptors

tcrdist2 is a python API-enabled toolkit for analyzing T-cell receptor repertoires.
Some of the functionality and code is adapted from the original tcr-dist
package which was released with the publication of  Dash et al.
Nature (2017) [doi:10.1038/nature22383](https://doi.org/10.1038/nature22383). This package contains a new API
for accessing the features of tcr-dist, as well as many new features that expand 
the T cell receptor analysis pipeline.

The original code for replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist).

---

## Documentation

Documentation, installation instructions, information about dependencies, and examples
can be found at  [tcrdist2.readthedocs.io](https://tcrdist2.readthedocs.io/en/latest/)

## Installation

The development version of tcrdist2 compatible with Python 3.6 or later. It
can be cloned or installed directly.

```bash
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2
```

To test out code used in the documented examples. Create a venv. Install tcrdist2 with all dependencies. 

```bash
python3 -m venv ./tenv
source tenv/bin/activate
pip install git+https://github.com/kmayerb/tcrdist2.git@0.0.4
pip install git+https://github.com/agartland/fishersapi.git#egg=fishersapi
```

Quickly test installation
```
python -c "import tcrdist"
```

## Development Files

If you wish, use the following set of commands to install testing files and legacy blast (replace the `download_from` argument to `'dropbox_osx'` or `'dropbox_linux'` based on your operating system).

```
python -c "import tcrdist as td; td.install_test_files.install_test_files()"
python -c "import tcrdist as td; td.setup_db.install_all_next_gen()"
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_osx')"
```

## Citing

##### Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas

[Nature (2017)](https://doi.org/10.1038/nature22383).



##### OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs

Zachary Sethna, Yuval Elhanati, Curtis G Callan, Aleksandra M Walczak, Thierry Mora

[Bioinformatics (2019)](https://doi.org/10.1093/bioinformatics/btz035)

(tcrdist2 incorporates OLGA for generation probability estimates)



##### Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments

Jeff Daily  

[BMC Bioinformatics (2016)](http://dx.doi.org/10.1186/s12859-016-0930-z)

(tcrdist2 depends on Parasail for fast sequence alignment)
