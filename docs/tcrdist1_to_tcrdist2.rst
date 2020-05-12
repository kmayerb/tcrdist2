.. _tcrdist1totcrdist2:

tcrdist1 to tcrdist2
====================

The example illustrates how to convert tcrdist1 files to tcrdist2 inputs.


This example uses the mouse clones files from the Dash et al. 2017 study.
(`clones <https://www.dropbox.com/s/l0z12f8lc752wfx/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv?dl=1>`_ 2MB).The clones file can be downloaded 
Searching for motifs in tcrdist2 is a slow step in the pipeline.
Thus for this example, users may wish to use precomputed candidate motifs.
A file with these motifs generated for TCRs specific to the PA epitope
(found in a protein from the mouse influenza virus) can be downloaded.
(`PA motifs <https://www.dropbox.com/s/z7wwmwb1n6dpq74/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log?dl=1>`_ 32KB).


Example
-------


The steps in the code block below:

1. Load the clones file as a Pandas DataFrame.
2. Select only those clones specific to the PA and F2 epitopes. Note, we choose two here, so we can later illustrate subsetting by epitope specificity.
3. Load an appropriate mapping object (i.e. a dictionary of old and new column names).
4. Apply the mapping object to the the original DataFrame, to select and rename columns.


.. code-block:: python

    import pandas as pd
    import numpy as np
    import tcrdist as td
    import IPython
    from tcrdist import mappers
    
    tcrdist_clone_fn = 'mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
    tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")               #1


    ind = (tcrdist_clone_df.epitope == "PA") | (tcrdist_clone_df.epitope == "F2") #2
    tcrdist_clone_df = tcrdist_clone_df[ind].copy()

    mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping                     #3
    tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,         #4
                                            mapping = mapping)


Writing Test Files
------------------

.. code-block:: python

    tcrdist2_df[["subject","epitope","count",
                "v_a_gene","j_a_gene","cdr3_a_aa","cdr3_a_nucseq",
                "v_b_gene","j_b_gene","cdr3_b_aa","cdr3_b_nucseq",
                "clone_id"]]. \
                to_csv("tcrdist/test_files_compact/dash_PA_F2.csv", index  = False)