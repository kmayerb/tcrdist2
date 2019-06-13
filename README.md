# TCRdist pipeline, version 2.0.0 

tcrdist2 is version 2.0.0 of TCRdist developed by Phil Bradley, Jeremy Crawford, and Andrew Fiore-Gartland. 

# Dependencies

It depends on python (2.7.11), numpy(1.10.1), scipy(0.16.0), scikit-learn(0.17.1), matplotlib(1.4.3), and parasail-python(1.1.16)



The current conda python development environment can be created quickly on a macOS or linux machine:

```bash
conda create --name tcrpy27osX python=2.7.11 scipy=0.16.0 matplotlib=1.4.3 numpy=1.10.1 futures=3.2.0 pandas=0.20.3 parasail-python=1.1.16 scikit-learn=0.17.1 jupyterlab jupyter
```
The current dependencies can be found here (tcrpy27osX.yml) and (tcrpy27linux.yml). 

# Installation

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

## Version 0.0.2 (09/21/2017):

- New sequence database system that makes it easier to work with alternate gene sets

- Preliminary support for gamma-delta TCRs: edit the `db_file` field of the `pipeline_params` dictionary
stored in `basic.py`. (This is a temporary hack; likely will move to switching by command line flag sometime soon).

- New `all_genes` dictionary in `all_genes.py` indexed by organism that holds information on all the genes; it's
read from the `db_file` pointed to by `basic.py`.

- With new minor updates to the probability model and default sequence database we're no longer trying to preserve
exact numerical identity with the results from the paper. To get the classic results, you can check out the
`version_001` branch on github.
