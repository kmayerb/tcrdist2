[![Build Status](https://travis-ci.com/kmayerb/tcrdist2.svg?branch=API2)](https://travis-ci.com/kmayerb/tcrdist2)
[![Build Docs Status](https://readthedocs.org/projects/tcrdist2/badge/?version=latest)](https://tcrdist2.readthedocs.io/en/latest/)

# tcrdist2

2019-07-10

tcrdist2 is a API-enabled toolkit expanding on the T-cell receptor analysis pipeline
developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and
colleagues as part of a T-cell receptor epitope specificity analysis
in Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383).

The original code for replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist).

---

## Documentation

Documentation, installation instructions, and examples
can be found here:

[tcrdist2.readthedocs.io](https://tcrdist2.readthedocs.io/en/latest/)

---

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



---

## tcrdist2 produces a distance measure based on comparison over multiple T-Cell Receptor complementarity-determining regions (CDRs)

Here is an example of what **tcrdist2** can do. For more information checkout [tcrdist2.readthedocs.io](https://tcrdist2.readthedocs.io/en/latest/)

#### Quick-Start Explanation

0. Load an example_df containing CDR3 sequences parsed from bulk or single-cell DNA sequencing reads
1. Initialize an instance of the `TCRrep` class, specifying `chains` and passing example_df to `cells_df.`
2. Optional: Add epitope and subject to the `TCRrep.index_cols` list
3. Use V-gene allele name to infer most probable amino acid sequences at the CDR1, CDR2, and pMHC loop position between the CDR2 and CDR3
4. Specify `index_cols` to indicate whcih CDRs should be used in the comparison
5. `Deduplicate` groups by index_cols and counts replicate cells (which might occur during clonal expansion).
The result is the `TCRrep.clones_df` object.
6. Optional: Specify an approriate substitution matrix for aligning each region
(if none are specified parasail.blosum62 will be used for all)
7. Run pairwise comparisons on CDR regions of the alpha chain.
(In this case, we use the Hamming Distance, a metric summing the number of mismatched aligned positions).
As a result, numpy pairwise distance matrices are stored as `TCRrep.cdr3_a_aa_pw` , `TCRrep.cdr2_a_aa_pw` , `TCRrep.cdr1_a_aa_pw`, and `TCRrep.pmhc_a_aa_pw1`
8. As in step 7, run pairwise comparisons on the CDR regions of the beta chain.
Numpy pairwise distance matrices are stored as `TCRrep.cdr3_b_aa_pw` , `TCRrep.cdr2_b_aa_pw` , `TCRrep.cdr1_b_aa_pw`, and `TCRrep.pmhc_b_aa_pw1`
9. Use `compute_paired_tcrdist` to sum over pairwise matrices at each CDR to calculate a multi-region tcrdist, supplying optional weights.


```python
import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep

# (0)
example_df = pd.DataFrame({'count': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1, 18: 1, 19: 1}, 'j_b_gene': {0: 'TRBJ1-2*01', 1: 'TRBJ1-2*01', 2: 'TRBJ1-2*01', 3: 'TRBJ1-2*01', 4: 'TRBJ1-2*01', 5: 'TRBJ1-2*01', 6: 'TRBJ1-2*01', 7: 'TRBJ1-5*01', 8: 'TRBJ1-2*01', 9: 'TRBJ1-2*01', 10: 'TRBJ1-2*01', 11: 'TRBJ1-2*01', 12: 'TRBJ1-2*01', 13: 'TRBJ1-2*01', 14: 'TRBJ2-3*01', 15: 'TRBJ1-5*01', 16: 'TRBJ2-7*01', 17: 'TRBJ1-1*01', 18: 'TRBJ2-7*01', 19: 'TRBJ2-7*01'}, 'j_a_gene': {0: 'TRAJ42*01', 1: 'TRAJ42*01', 2: 'TRAJ42*01', 3: 'TRAJ50*01', 4: 'TRAJ42*01', 5: 'TRAJ42*01', 6: 'TRAJ42*01', 7: 'TRAJ20*01', 8: 'TRAJ42*01', 9: 'TRAJ42*01', 10: 'TRAJ42*01', 11: 'TRAJ42*01', 12: 'TRAJ42*01', 13: 'TRAJ42*01', 14: 'TRAJ49*01', 15: 'TRAJ33*01', 16: 'TRAJ42*01', 17: 'TRAJ49*01', 18: 'TRAJ31*01', 19: 'TRAJ37*02'}, 'cdr3_a_aa': {0: 'CAGQASQGNLIF', 1: 'CAGQASQGNLIF', 2: 'CAGQASQGNLIF', 3: 'CAGPRETSYDKVIF', 4: 'CAGQASQGNLIF', 5: 'CAGQASQGNLIF', 6: 'CAGQASQGNLIF', 7: 'CAETRSRDYKLSF', 8: 'CAGQASQGNLIF', 9: 'CAGQASQGNLIF', 10: 'CAGQASQGNLIF', 11: 'CAGQASQGNLIF', 12: 'CAGQASQGNLIF', 13: 'CAGQASQGNLIF', 14: 'CAVADTGNQFYF', 15: 'CLVGSMDSNYQLIW', 16: 'CAVPKGSQGNLIF', 17: 'CAVSDSGTGNQFYF', 18: 'CAGPFGRLMF', 19: 'CAGPDGSSNTGKLIF'}, 'epitope': {0: 'pp65', 1: 'pp65', 2: 'pp65', 3: 'pp65', 4: 'pp65', 5: 'pp65', 6: 'pp65', 7: 'pp65', 8: 'pp65', 9: 'pp65', 10: 'pp65', 11: 'pp65', 12: 'pp65', 13: 'pp65', 14: 'pp65', 15: 'M1', 16: 'M1', 17: 'M1', 18: 'M1', 19: 'M1'}, 'cdr3_b_aa': {0: 'CASSIQALLTF', 1: 'CASSIQALLTF', 2: 'CASSIQALLTF', 3: 'CASSSAYYGYTF', 4: 'CASSIQALLTF', 5: 'CASSIQALLTF', 6: 'CASSIQALLTF', 7: 'CASSQEEGPGNQPQHF', 8: 'CASSIQALLTF', 9: 'CASSIQALLTF', 10: 'CASSIQALLTF', 11: 'CASSIQALLTF', 12: 'CASSIQALLTF', 13: 'CASSIQALLTF', 14: 'CATAITSTQYF', 15: 'CASSSQSNQPQHF', 16: 'CASSIRSSYEQYF', 17: 'CASSQMTGLNTEAFF', 18: 'CASSLFPGFGEQYF', 19: 'CASSLIFPSGEQYF'}, 'v_b_gene': {0: 'TRBV12-3*01', 1: 'TRBV12-3*01', 2: 'TRBV12-3*01', 3: 'TRBV12-3*01', 4: 'TRBV12-3*01', 5: 'TRBV12-3*01', 6: 'TRBV12-3*01', 7: 'TRBV4-1*01', 8: 'TRBV12-3*01', 9: 'TRBV12-3*01', 10: 'TRBV12-3*01', 11: 'TRBV12-3*01', 12: 'TRBV12-3*01', 13: 'TRBV12-3*01', 14: 'TRBV12-3*01', 15: 'TRBV25-1*01', 16: 'TRBV19*01', 17: 'TRBV28*01', 18: 'TRBV27*01', 19: 'TRBV27*01'}, 'id': {0: 'human_tcr0001', 1: 'human_tcr0002', 2: 'human_tcr0003', 3: 'human_tcr0004', 4: 'human_tcr0005', 5: 'human_tcr0006', 6: 'human_tcr0007', 7: 'human_tcr0008', 8: 'human_tcr0009', 9: 'human_tcr0010', 10: 'human_tcr0011', 11: 'human_tcr0012', 12: 'human_tcr0013', 13: 'human_tcr0014', 14: 'human_tcr0015', 15: 'human_tcr0016', 16: 'human_tcr0017', 17: 'human_tcr0018', 18: 'human_tcr0019', 19: 'human_tcr0020'}, 'v_a_gene': {0: 'TRAV35*01', 1: 'TRAV35*01', 2: 'TRAV35*01', 3: 'TRAV35*02', 4: 'TRAV35*01', 5: 'TRAV35*01', 6: 'TRAV35*01', 7: 'TRAV5*01', 8: 'TRAV35*01', 9: 'TRAV35*01', 10: 'TRAV35*01', 11: 'TRAV35*01', 12: 'TRAV35*01', 13: 'TRAV35*01', 14: 'TRAV22*01', 15: 'TRAV4*01', 16: 'TRAV8-3*02', 17: 'TRAV8-6*02', 18: 'TRAV27*01', 19: 'TRAV35*02'}, 'subject': {0: 'human_subject0010', 1: 'human_subject0010', 2: 'human_subject0010', 3: 'human_subject0010', 4: 'human_subject0010', 5: 'human_subject0010', 6: 'human_subject0010', 7: 'human_subject0010', 8: 'human_subject0010', 9: 'human_subject0010', 10: 'human_subject0010', 11: 'human_subject0010', 12: 'human_subject0010', 13: 'human_subject0010', 14: 'human_subject0010', 15: 'human_subject0007', 16: 'human_subject0015', 17: 'human_subject0007', 18: 'human_subject0007', 19: 'human_subject0007'}})
example_df.head()

# (1)
tr = TCRrep(cell_df = example_df.copy(),  
            organism = "human",
            chains= ["alpha","beta"])
# (2)
tr.index_cols.append('epitope')
tr.index_cols.append('subject')
# (3)
tr.infer_cdrs_from_v_gene(chain = "alpha")
tr.infer_cdrs_from_v_gene(chain = "beta")
# (4)
tr.index_cols =['cdr3_a_aa',
                'cdr1_a_aa',
                'cdr2_a_aa',
                'pmhc_a_aa',
                'cdr3_b_aa',
                'cdr1_b_aa',
                'cdr2_b_aa',
                'pmhc_b_aa']
# (5)
tr.deduplicate()
tr.clone_df.head()

# (7)
tr.compute_pairwise_all(chain = "alpha",
                        metric = "hamming")
# (8)
tr.compute_pairwise_all(chain = "beta",
                        metric = "hamming")
# (9)
tr.compute_paired_tcrdist(chains = ['alpha','beta'],
                          replacement_weights= {'cdr3_a_aa_pw': 3,
                                                'cdr3_b_aa_pw': 3})
```




    {'paired_tcrdist': array([[  0.,  86.,  83.,  84.,  87.,  93.,  89.,  80.,  78.],
            [ 86.,   0.,  43.,  59.,  59.,  81.,  60.,  81.,  80.],
            [ 83.,  43.,   0.,  78.,  69.,  81.,  68.,  80.,  90.],
            [ 84.,  59.,  78.,   0.,  42.,  73.,  76.,  92.,  82.],
            [ 87.,  59.,  69.,  42.,   0.,  58.,  67.,  83.,  76.],
            [ 93.,  81.,  81.,  73.,  58.,   0.,  74.,  69.,  84.],
            [ 89.,  60.,  68.,  76.,  67.,  74.,   0.,  70.,  71.],
            [ 80.,  81.,  80.,  92.,  83.,  69.,  70.,   0.,  85.],
            [ 78.,  80.,  90.,  82.,  76.,  84.,  71.,  85.,   0.]]),
     'paired_tcrdist_weights': {'cdr1_a_aa_pw': 1,
      'cdr1_b_aa_pw': 1,
      'cdr2_a_aa_pw': 1,
      'cdr2_b_aa_pw': 1,
      'cdr3_a_aa_pw': 3,
      'cdr3_b_aa_pw': 3,
      'pmhc_a_aa_pw': 1,
      'pmhc_b_aa_pw': 1}}

---


# Citing

Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas

Nature (2017) doi:10.1038/nature22383


Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z

http://dx.doi.org/10.1186/s12859-016-0930-z

---
