[![Build Status](https://travis-ci.com/kmayerb/tcrdist2.svg?branch=API2)](https://travis-ci.com/kmayerb/tcrdist2)

# tcrdist2

2019-07-04

**tcrdist2** is a API-enabled toolkit of T-Cell Receptor analysis inspired and built on the TCRdist pipeline developed by Phillip Harlan Bradley, Jeremy Chase Crawford, and colleagues as part of a T-cell receptor epitope specificity analysis in Dash et al. [Nature (2017) doi:10.1038/nature22383](https://doi.org/10.1038/nature22383). The original code for replicating analysis performed in the manuscript can be found [here](https://github.com/phbradley/tcr-dist).

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




## tcrdist2 produces a distance measure based on comparison over multiple T-Cell Receptor complementarity-determining regions (CDRs)

Here is an example of what **tcrdist2** can do. A more detailed explanation of
the tools and their customization follows below.

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




# tcrdist2 Vignette: Tools for Computing Pairwise Distance

In tcrdist2, `TCRrep` (T-Cell Receptor Repertoire) is the main object class for repertoire analysis. This vignette illustrates its use.

```python
from collections import OrderedDict
import tcrdist as td
from tcrdist.repertoire import TCRrep
import pandas as pd
import numpy as np
import parasail
```

## Preliminary: Load Example Data


```python
df = td.processing.readPairedSequences(paired_seqs_file = "tcrdist/datasets/test_human_pairseqs.tsv",
                                         organism = "human", use_parasail = True);

tcrdist_to_tcrdist2_mapping = OrderedDict([('id' , 'id'),
                                           ('epitope' , 'epitope'),
                                           ('subject' , 'subject'),
                                           ('cdr3a' , 'cdr3_a_aa'),
                                           ('cdr3b' , 'cdr3_b_aa'),
                                           ('ja_gene', 'j_a_gene'),
                                           ('va_gene','v_a_gene'),
                                           ('jb_gene', 'j_b_gene'),
                                           ('vb_gene', 'v_b_gene')])
example_df = df[tcrdist_to_tcrdist2_mapping.keys()].rename(columns = tcrdist_to_tcrdist2_mapping)
example_df.head()
example_df['count'] = 1
```
```python
example_df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>epitope</th>
      <th>subject</th>
      <th>cdr3_a_aa</th>
      <th>cdr3_b_aa</th>
      <th>j_a_gene</th>
      <th>v_a_gene</th>
      <th>j_b_gene</th>
      <th>v_b_gene</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>human_tcr0001</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>human_tcr0002</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>human_tcr0003</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>human_tcr0004</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGPRETSYDKVIF</td>
      <td>CASSSAYYGYTF</td>
      <td>TRAJ50*01</td>
      <td>TRAV35*02</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>human_tcr0005</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




### `example_df` Data Columns

Data for paired alpha/beta chained single cell TCR data will contain the following headers

#### Columns:
#### identifying data
- id - (e.g. human_tcr0001) unique identifier
- subject - (e.g. human_subject0001) subject from whom the cell was isolated
- epitope - (e.g. pp65) epitope used to bind TCR

#### alpha chain data
- cdr3_a_aa - (e.g. CAGQASQGNLIF) cdr3 alpha amino acid sequence
- v_a_gene - (e.g. TRAV35*01) V gene usage for alpha chain
- j_a_gene - (e.g. TRAJ42*01) J gene usage for alpha chain

#### beta chain data
- cdr3_b_aa - cdr3 beta amino acid sequence
- v_b_gene -  V gene usage for beta chain
- j_b_gene - J gene usage for beta chain

# Using the `TCRrep()` Class:

### 1. Initialize an Instance of the `TCRrep()` Class

The user must provide the cell_df (pandas.DataFrame) and specify with a list which chains are to be analyzed. The only currently supported choices for the `chains` argument are:
- ["alpha"],
- ["beta"],
- ["gamma"],
- ["delta"],
- ["alpha", "beta"]
- ["gamma", "delta"]


```python
tcr = TCRrep( cell_df = example_df, chains = ["alpha", "beta"])
tcr
```




    tcrdist.repertoire.TCRrep for <Your TCR Repertoire Project>
     with index_cols: ['cdr3_a_aa', 'cdr3_b_aa']



### 2. Run `TCRrep.deduplicate()`

Running `.deduplicate` finds redundant copies grouped by index columns and creates the non-redundant `clones_df`. Notice how the count column in `clones_df` shows how frequently the alpha/beta pairs occurred.


```python
tcr.deduplicate()
tcr
```




    tcrdist.repertoire.TCRrep for <Your TCR Repertoire Project>
     with index_cols: ['cdr3_a_aa', 'cdr3_b_aa']




```python
tcr.clone_df
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cdr3_a_aa</th>
      <th>cdr3_b_aa</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CAETRSRDYKLSF</td>
      <td>CASSQEEGPGNQPQHF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CAGPDGSSNTGKLIF</td>
      <td>CASSLIFPSGEQYF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CAGPFGRLMF</td>
      <td>CASSLFPGFGEQYF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CAGPRETSYDKVIF</td>
      <td>CASSSAYYGYTF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>12</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CAVADTGNQFYF</td>
      <td>CATAITSTQYF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CAVPKGSQGNLIF</td>
      <td>CASSIRSSYEQYF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CAVSDSGTGNQFYF</td>
      <td>CASSQMTGLNTEAFF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>CLVGSMDSNYQLIW</td>
      <td>CASSSQSNQPQHF</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



#### 2b. Notice that epitope or subject fields can be added prior to running .`deduplicate()`.


```python
tcr = TCRrep( cell_df = example_df, chains = ["alpha", "beta"])
tcr.index_cols.append("epitope")
tcr.index_cols.append("subject")
tcr
tcr.deduplicate()
tcr.clones_df
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cdr3_a_aa</th>
      <th>cdr3_b_aa</th>
      <th>epitope</th>
      <th>subject</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CAETRSRDYKLSF</td>
      <td>CASSQEEGPGNQPQHF</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CAGPDGSSNTGKLIF</td>
      <td>CASSLIFPSGEQYF</td>
      <td>M1</td>
      <td>human_subject0007</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CAGPFGRLMF</td>
      <td>CASSLFPGFGEQYF</td>
      <td>M1</td>
      <td>human_subject0007</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CAGPRETSYDKVIF</td>
      <td>CASSSAYYGYTF</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>12</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CAVADTGNQFYF</td>
      <td>CATAITSTQYF</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CAVPKGSQGNLIF</td>
      <td>CASSIRSSYEQYF</td>
      <td>M1</td>
      <td>human_subject0015</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CAVSDSGTGNQFYF</td>
      <td>CASSQMTGLNTEAFF</td>
      <td>M1</td>
      <td>human_subject0007</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>CLVGSMDSNYQLIW</td>
      <td>CASSSQSNQPQHF</td>
      <td>M1</td>
      <td>human_subject0007</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



### 3. Run `.compute_pairwise()` on alpha chains

The current default option for calculating pairwise distances is Needleman-Wunsch alignment using the blosum62 substitution matrix. Penalties for a gap in the alignment (open: 3) and a further gap extension (extend: 3).  When `metric = "nw"` Distance are calculated from alignment score using the following formula:

```python
xx = parasail.nw_stats(s1, s1, open=open, extend=extend, matrix=matrix).score
yy = parasail.nw_stats(s2, s2, open=open, extend=extend, matrix=matrix).score
xy = parasail.nw_stats(s1, s2, open=open, extend=extend, matrix=matrix).score

D = xx + yy - 2 * xy
```


Compute_pairwise is run independently on each chain to allow user to customize the distance metric to different regions of the T-cell receptor (We elaborate on these options in the next section).

#### 3a. Run `.compute_pairwise()` on alpha chains


```python
tcr.compute_pairwise(chain = "alpha", metric="nw")
tcr.cdr3_a_aa_pw
```




    array([[   0.,  113.,  103.,   90.,   96.,  107.,  103.,  123.,  114.],
           [ 113.,    0.,   84.,   85.,   75.,  112.,   72.,  106.,  115.],
           [ 103.,   84.,    0.,   83.,   67.,   94.,   82.,  112.,  119.],
           [  90.,   85.,   83.,    0.,   84.,  111.,   93.,  127.,  106.],
           [  96.,   75.,   67.,   84.,    0.,   81.,   45.,   93.,  106.],
           [ 107.,  112.,   94.,  111.,   81.,    0.,   90.,   28.,  109.],
           [ 103.,   72.,   82.,   93.,   45.,   90.,    0.,   86.,  105.],
           [ 123.,  106.,  112.,  127.,   93.,   28.,   86.,    0.,  119.],
           [ 114.,  115.,  119.,  106.,  106.,  109.,  105.,  119.,    0.]])



#### 3b. Run `.compute_pairwise()` on beta chains




```python
tcr.compute_pairwise(chain = "beta", metric="nw")
tcr.cdr3_b_aa_pw
```




    array([[   0.,  109.,  107.,  133.,  135.,  132.,  116.,  108.,   61.],
           [ 109.,    0.,   32.,   90.,   98.,   73.,   59.,  103.,   90.],
           [ 107.,   32.,    0.,   92.,  102.,   89.,   63.,   89.,  102.],
           [ 133.,   90.,   92.,    0.,   70.,   89.,   69.,  107.,   96.],
           [ 135.,   98.,  102.,   70.,    0.,   71.,   75.,   89.,   86.],
           [ 132.,   73.,   89.,   89.,   71.,    0.,   56.,   84.,   83.],
           [ 116.,   59.,   63.,   69.,   75.,   56.,    0.,   92.,   67.],
           [ 108.,  103.,   89.,  107.,   89.,   84.,   92.,    0.,  101.],
           [  61.,   90.,  102.,   96.,   86.,   83.,   67.,  101.,    0.]])



### 4. Combine Pairwise Distances
Because each region is represented as a pairwise matrix weights can be applied to given specific regions greater emphasis.


```python
tcrdist = 1*tcr.cdr3_a_aa_pw + 1*tcr.cdr3_b_aa_pw
tcrdist
```




    array([[   0.,  222.,  210.,  223.,  231.,  239.,  219.,  231.,  175.],
           [ 222.,    0.,  116.,  175.,  173.,  185.,  131.,  209.,  205.],
           [ 210.,  116.,    0.,  175.,  169.,  183.,  145.,  201.,  221.],
           [ 223.,  175.,  175.,    0.,  154.,  200.,  162.,  234.,  202.],
           [ 231.,  173.,  169.,  154.,    0.,  152.,  120.,  182.,  192.],
           [ 239.,  185.,  183.,  200.,  152.,    0.,  146.,  112.,  192.],
           [ 219.,  131.,  145.,  162.,  120.,  146.,    0.,  178.,  172.],
           [ 231.,  209.,  201.,  234.,  182.,  112.,  178.,    0.,  220.],
           [ 175.,  205.,  221.,  202.,  192.,  192.,  172.,  220.,    0.]])



### Default Steps Combined


```python
tcr = TCRrep( cell_df = example_df, chains = ["alpha", "beta"])
tcr.index_cols.append("epitope")
tcr.index_cols.append("subject")
tcr.deduplicate()
tcr.compute_pairwise(chain = "alpha", metric="nw")
tcr.compute_pairwise(chain = "beta", metric="nw")
tcrdist = tcr.cdr3_a_aa_pw + tcr.cdr3_b_aa_pw
tcrdist
```




    array([[   0.,  222.,  210.,  223.,  231.,  239.,  219.,  231.,  175.],
           [ 222.,    0.,  116.,  175.,  173.,  185.,  131.,  209.,  205.],
           [ 210.,  116.,    0.,  175.,  169.,  183.,  145.,  201.,  221.],
           [ 223.,  175.,  175.,    0.,  154.,  200.,  162.,  234.,  202.],
           [ 231.,  173.,  169.,  154.,    0.,  152.,  120.,  182.,  192.],
           [ 239.,  185.,  183.,  200.,  152.,    0.,  146.,  112.,  192.],
           [ 219.,  131.,  145.,  162.,  120.,  146.,    0.,  178.,  172.],
           [ 231.,  209.,  201.,  234.,  182.,  112.,  178.,    0.,  220.],
           [ 175.,  205.,  221.,  202.,  192.,  192.,  172.,  220.,    0.]])



## Customizing the Pairwise Default Methods

There are a number of ways to customize the computation of pairwise distances.

1. Distances can be calculated via reciprocal alignment scores.
2. Distance can also be calculated on aligned strings (i.e. the Hamming Distance)
3.  Any custom distance metric taking two strings can be supplied by the user.

### 1. Customizing Alignment Parameters

When pairwise distance is calculated using alignment scores, the amino acid substitution matrix and gap open/extension penalties influence the optimal alignment and resulting scores. These parameters are available in  `compute_pairwise()` method using `matrix`, `open` and `extend` arguments. When `metric = "nw"` Distance are calculated from alignment score using the following formula:

```
xx = parasail.nw_stats(s1, s1, open=open, extend=extend, matrix=matrix).score
yy = parasail.nw_stats(s2, s2, open=open, extend=extend, matrix=matrix).score
xy = parasail.nw_stats(s1, s2, open=open, extend=extend, matrix=matrix).score

D = xx + yy - 2 * xy
```

The default `compute_pairwise()`:
```
tcr.compute_pairwise(chain = "alpha", metric="nw", matrix = parasail.blosum62, open = 3, extend = 3)
```

`compute_pairwise()` customized with gap penalty = 8:
```
tcr.compute_pairwise(chain = "alpha", metric="nw", matrix = parasail.blosum62, open = 8, extend = 8)
```

`compute_pairwise()` customized with a different substitution matrix:
```
tcr.compute_pairwise(chain = "alpha", metric="nw", matrix = parasail.pam10, open = 8, extend = 8)
```



```python
tcr = TCRrep( cell_df = example_df, chains = ["alpha", "beta"])
tcr.deduplicate()
```




    tcrdist.repertoire.TCRrep for <Your TCR Repertoire Project>
     with index_cols: ['cdr3_a_aa', 'cdr3_b_aa']



#### `tcr.compute_pairwise(chain = 'alpha', metric='nw'`


```python
tcr.compute_pairwise(chain = "alpha", metric="nw")
default_pw = tcr.cdr3_a_aa_pw.copy()
print("tcr.compute_pairwise(chain = 'alpha', metric='nw'")
print(default_pw)
```

    tcr.compute_pairwise(chain = 'alpha', metric='nw'
    [[   0.  113.  103.   90.   96.  107.  103.  123.  114.]
     [ 113.    0.   84.   85.   75.  112.   72.  106.  115.]
     [ 103.   84.    0.   83.   67.   94.   82.  112.  119.]
     [  90.   85.   83.    0.   84.  111.   93.  127.  106.]
     [  96.   75.   67.   84.    0.   81.   45.   93.  106.]
     [ 107.  112.   94.  111.   81.    0.   90.   28.  109.]
     [ 103.   72.   82.   93.   45.   90.    0.   86.  105.]
     [ 123.  106.  112.  127.   93.   28.   86.    0.  119.]
     [ 114.  115.  119.  106.  106.  109.  105.  119.    0.]]


#### `tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.blosum62, open = 8, extend = 8)`


```python
tcr.compute_pairwise(chain = "alpha", metric="nw", matrix = parasail.blosum62, open = 8, extend = 8)
custom_pw_blosum62 = tcr.cdr3_a_aa_pw.copy()
print("tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.blosum62, open = 8, extend = 8)")
print(custom_pw_blosum62)
```

    tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.blosum62, open = 8, extend = 8)
    [[   0.  133.  133.  100.  106.  125.  107.  133.  128.]
     [ 133.    0.  134.   95.  105.  142.   92.  130.  141.]
     [ 133.  134.    0.  123.   87.  114.  112.  156.  159.]
     [ 100.   95.  123.    0.  104.  131.  103.  137.  126.]
     [ 106.  105.   87.  104.    0.  101.   55.  113.  126.]
     [ 125.  142.  114.  131.  101.    0.  116.   48.  133.]
     [ 107.   92.  112.  103.   55.  116.    0.  102.  127.]
     [ 133.  130.  156.  137.  113.   48.  102.    0.  129.]
     [ 128.  141.  159.  126.  126.  133.  127.  129.    0.]]


#### `tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.pam10, open = 8, extend = 8)`


```python
tcr.compute_pairwise(chain = "alpha", metric="nw", matrix = parasail.pam10, open = 8, extend = 8)
custom_pw_pam10 = tcr.cdr3_a_aa_pw.copy()
print("tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.pam10, open = 8, extend = 8)")
print(custom_pw_pam10)
```

    tcr.compute_pairwise(chain = 'alpha', metric='nw', matrix = parasail.pam10, open = 8, extend = 8)
    [[   0.  265.  231.  231.  227.  269.  248.  315.  299.]
     [ 265.    0.  200.  192.  176.  258.  167.  236.  268.]
     [ 231.  200.    0.  216.  170.  250.  201.  278.  292.]
     [ 231.  192.  216.    0.  216.  270.  233.  310.  274.]
     [ 227.  176.  170.  216.    0.  208.  107.  236.  262.]
     [ 269.  258.  250.  270.  208.    0.  203.   66.  270.]
     [ 248.  167.  201.  233.  107.  203.    0.  205.  245.]
     [ 315.  236.  278.  310.  236.   66.  205.    0.  278.]
     [ 299.  268.  292.  274.  262.  270.  245.  278.    0.]]


## 2. Hamming Distance
An alternative to distance computed by alignment scores is to compute a metric on the aligned strings. The Hamming Distance between two strings or vectors of equal length is the number of positions with mismatching information. The method `.compute_pairwise()` can return Hamming Distance using `metric = 'hamming'`.


```python
tcr.compute_pairwise(chain = "alpha", metric="hamming")
tcr.cdr3_a_aa_pw
```




    array([[  0.,   9.,   8.,   8.,   8.,  10.,   9.,  10.,  10.],
           [  9.,   0.,   8.,   7.,   7.,  10.,   6.,   9.,  10.],
           [  8.,   8.,   0.,   9.,   6.,   8.,   7.,  10.,  11.],
           [  8.,   7.,   9.,   0.,   8.,  10.,   8.,  11.,  11.],
           [  8.,   7.,   6.,   8.,   0.,   7.,   4.,   9.,   9.],
           [ 10.,  10.,   8.,  10.,   7.,   0.,   8.,   3.,  10.],
           [  9.,   6.,   7.,   8.,   4.,   8.,   0.,   8.,   9.],
           [ 10.,   9.,  10.,  11.,   9.,   3.,   8.,   0.,  10.],
           [ 10.,  10.,  11.,  11.,   9.,  10.,   9.,  10.,   0.]])



#### An Visual Example of Hamming Distance


```python
algn = parasail.nw_trace("CAGQASQGNLIF", "CATTAQASQGNLIF", open = 3, extend = 3, matrix = parasail.blosum62)
print(algn.traceback.query)
print(algn.traceback.comp)
print(algn.traceback.ref)
print("HAMMING DISTANCE" , td.SequencePair("CAGQASQGNLIF", "CATTAQASQGNLIF").hamming_distance)
```

    CA--GQASQGNLIF
    ||  .|||||||||
    CATTAQASQGNLIF
    ('HAMMING DISTANCE', 3.0)


## 3. Custom Distance

The method `.compute_pairwise()` permits the user to define a custom function that computes a distance based on any input strings. Below is a trivial example. Once the function is defined, it can be used by setting `metric = 'custom'` and supplying the function to the `user_function` argument:

```
tcr.compute_pairwise(chain = "alpha", metric= "custom", user_function= my_user_function)
```


```python
# A trivial example of custom metric
def my_user_function(s1,s2):
    if s1.startswith("CAV") and s2.startswith("CAV"):
        return(0)
    else:
        return(1)
```


```python
tcr = TCRrep( cell_df = example_df, chains = ["alpha", "beta"])
tcr.deduplicate()
tcr.compute_pairwise(chain = "alpha", metric= "custom", user_function= my_user_function)
tcr.cdr3_a_aa_pw
```




    array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]])



# 4. Multiprocessing
The method `.compute_pairwise()` has been designed to take advantage of multiple CPUs to reduced computation time by parallelization. The default number of processes is 2. However for large jobs. The user can take advantage of more CPUs as follows:

```
import multiprocessing
tcr.compute_pairwise(chain = "alpha", metric= "nw", processes = multiprocessing.cpu_count()
```

# 5. Generating Example Data

### Example Data

Data for paired alpha/beta chained single cell TCR data will contain the following headers

#### Columns:
#### identifying data
- id - (e.g. human_tcr0001) unique identifier
- subject - (e.g. human_subject0001) subject from whom the cell was isolated
- epitope - (e.g. pp65) epitope used to bind TCR

#### alpha chain data
- cdr3_a_aa - (e.g. CAGQASQGNLIF) cdr3 alpha amino acid sequence
- v_a_gene - (e.g. TRAV35*01) V gene usage for alpha chain
- j_a_gene - (e.g. TRAJ42*01) J gene usage for alpha chain

#### beta chain data
- cdr3_b_aa - cdr3 beta amino acid sequence
- v_b_gene -  V gene usage for beta chain
- j_b_gene - J gene usage for beta chain

## 5a. data generated from paired-reads (legacy method)
tcrdist2 has method `td.processing.readPairedSequences()` for parsing full length paired-chain sequences. These tools are discussed elsewhere. Because of complexity of parsing task, these operations take some time and are not recommended on more than 5000 sequences. We will parallelize these method in future versions.


```python
import tcrdist as td
from collections import OrderedDict

df = td.processing.readPairedSequences(paired_seqs_file = "tcrdist/datasets/test_human_pairseqs.tsv",
                                         organism = "human", use_parasail = True);

tcrdist_to_tcrdist2_mapping = OrderedDict([('id' , 'id'),
                                           ('epitope' , 'epitope'),
                                           ('subject' , 'subject'),
                                           ('cdr3a' , 'cdr3_a_aa'),
                                           ('cdr3b' , 'cdr3_b_aa'),
                                           ('ja_gene', 'j_a_gene'),
                                           ('va_gene','v_a_gene'),
                                           ('jb_gene', 'j_b_gene'),
                                           ('vb_gene', 'v_b_gene')])
example_df = df[tcrdist_to_tcrdist2_mapping.keys()].rename(columns = tcrdist_to_tcrdist2_mapping)
example_df['count'] = 1
```

    RESULTS BASED ON PARASAIL




```python
#pd.to_pickle(path="./quick_load_example_df", obj=example_df)
```


```python
example_df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>epitope</th>
      <th>subject</th>
      <th>cdr3_a_aa</th>
      <th>cdr3_b_aa</th>
      <th>j_a_gene</th>
      <th>v_a_gene</th>
      <th>j_b_gene</th>
      <th>v_b_gene</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>human_tcr0001</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>human_tcr0002</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>human_tcr0003</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>human_tcr0004</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGPRETSYDKVIF</td>
      <td>CASSSAYYGYTF</td>
      <td>TRAJ50*01</td>
      <td>TRAV35*02</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>human_tcr0005</td>
      <td>pp65</td>
      <td>human_subject0010</td>
      <td>CAGQASQGNLIF</td>
      <td>CASSIQALLTF</td>
      <td>TRAJ42*01</td>
      <td>TRAV35*01</td>
      <td>TRBJ1-2*01</td>
      <td>TRBV12-3*01</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



## 5b. data generated from Adaptive Immunoseq Data


```python
from collections import OrderedDict
adaptive_delta_to_tcrdist2_mapping = OrderedDict( [('count (templates/reads)', 'count'),
                                                   ('aminoAcid' , 'cdr3_d_aa'),
                                                   ('jMaxResolved', 'j_d_gene'),
                                                   ('vMaxResolved','v_d_gene')])
```


```python
adaptive_df = pd.read_csv("tcrdist/datasets/adaptive_delta_template.tsv", sep = "\t")
example_adaptive_df = adaptive_df[adaptive_delta_to_tcrdist2_mapping.keys()].\
    rename(columns = adaptive_delta_to_tcrdist2_mapping)
example_adaptive_df.head(10)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>count</th>
      <th>cdr3_d_aa</th>
      <th>j_d_gene</th>
      <th>v_d_gene</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>954</td>
      <td>CACEMLGHPPGDKLIF</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>1</th>
      <td>234</td>
      <td>CACDTAAGGYASSWDTRQMFF</td>
      <td>TCRDJ03-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>2</th>
      <td>369</td>
      <td>CACDYVLGAEDKLIF</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>3</th>
      <td>170</td>
      <td>CACDTLFLGEDTPTDKLIF</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>4</th>
      <td>113</td>
      <td>CACDIVLSGGLDTRQMFF</td>
      <td>TCRDJ03-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>5</th>
      <td>184</td>
      <td>CACDNLSETTDKLIF</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>6</th>
      <td>356</td>
      <td>NaN</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV01-01*01</td>
    </tr>
    <tr>
      <th>7</th>
      <td>153</td>
      <td>CACDTIRGFSSWDTRQMFF</td>
      <td>TCRDJ03-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>8</th>
      <td>143</td>
      <td>CACDTGRLLGDTADTRQMFF</td>
      <td>TCRDJ03-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
    <tr>
      <th>9</th>
      <td>249</td>
      <td>CACDPLGDNDKLIF</td>
      <td>TCRDJ01-01*01</td>
      <td>TCRDV02-01</td>
    </tr>
  </tbody>
</table>
</div>




```python
example_adaptive_df = example_adaptive_df[example_adaptive_df['cdr3_d_aa'].notnull()].reset_index()
```

### tcrdist on single-chain bulk adaptive seq


```python
tcrAdaptive = TCRrep(example_adaptive_df, chains = ["delta"])
tcrAdaptive.deduplicate()
tcrAdaptive.clones_df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cdr3_d_aa</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CAAGYTGNQFYF</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CACADLGAYPDKLIF</td>
      <td>67</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CACDDVTEVEGDKLIF</td>
      <td>43</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CACDFISPSNWGIQSGRNTDKLIF</td>
      <td>42</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CACDILLGDTADKLIF</td>
      <td>25</td>
    </tr>
  </tbody>
</table>
</div>




```python
tcrAdaptive.compute_pairwise(chain="delta", metric = "nw")
tcrAdaptive.cdr3_d_aa_pw[1:10,1:10]
```




    array([[   0.,   85.,  169.,   73.,  119.,  127.,   68.,   91.,   73.],
           [  85.,    0.,  142.,   82.,  114.,  120.,   67.,   88.,   56.],
           [ 169.,  142.,    0.,  144.,  158.,  168.,  139.,  130.,  132.],
           [  73.,   82.,  144.,    0.,   84.,  108.,   39.,   54.,   52.],
           [ 119.,  114.,  158.,   84.,    0.,   76.,   95.,  104.,  110.],
           [ 127.,  120.,  168.,  108.,   76.,    0.,   97.,  108.,  118.],
           [  68.,   67.,  139.,   39.,   95.,   97.,    0.,   69.,   45.],
           [  91.,   88.,  130.,   54.,  104.,  108.,   69.,    0.,   60.],
           [  73.,   56.,  132.,   52.,  110.,  118.,   45.,   60.,    0.]])



---

# More Information


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
python -c "import tcrdist as td; td.setup_blast.install_blast_to_externals(download_from = 'dropbox_linux);"
```

# Citing

Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas

Nature (2017) doi:10.1038/nature22383


Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z

http://dx.doi.org/10.1186/s12859-016-0930-z

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
