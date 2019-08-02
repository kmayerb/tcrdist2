[![Build Status](https://travis-ci.com/kmayerb/tcrdist2.svg?branch=API2)](https://travis-ci.com/kmayerb/tcrdist2)
[![Build Docs Status](https://readthedocs.org/projects/tcrdist2/badge/?version=latest)](https://tcrdist2.readthedocs.io/en/latest/)

# tcrdist2

2019-07-12

#### tcrdist2 produces flexible distance measures for comparing T cell receptors.

tcrdist2 is a python API-enabled toolkit expanding on the T cell receptor analysis pipeline
developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T cell receptor epitope specificity analysis
in Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383).

The original code for replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist).

---

## Documentation

Documentation, installation instructions, information about dependencies, and examples
can be found at  [tcrdist2.readthedocs.io](https://tcrdist2.readthedocs.io/en/latest/)

## Installation

The development version of tcrdist2 compatible with Python 2.7 or Python 3.6
can be cloned or installed directly.

```bash
  pip install git+https://github.com/kmayerb/tcrdist2.git@API2
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
