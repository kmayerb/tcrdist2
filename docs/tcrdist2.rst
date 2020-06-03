tcrdist2 
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:


This page introduces some quick ways to get started with tcrdist2. You can calculate 
tcrdistances with only a few lines of code. For some test epitope specifice, paired-chain T cell recptor data, 
download `dash.csv <https://github.com/kmayerb/tcrdist2/raw/API2/tcrdist/datasets/dash.csv>`_. 

tcrdist2 Distance
-----------------

.. code-block:: python

    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv('dash.csv')
    tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
    tr.tcrdist2()

You can then access pairwise distance attributes:

.. code-block:: python

    tr.pw_tcrdist
    tr.pw_alpha
    tr.pw_beta


tcrdist2 Basic Options
++++++++++++++++++++++


:py:meth:`tcrdist.repertoire.TCRrep.tcrdist2` provides multiple options.

.. code-block:: python
  
    import multiprocessing
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    cpu = multiprocessing.cpu_count() # determine the number of available cpu

    df = pd.read_csv('dash.csv')
    df = df[df.epitope.isin(['NP'])]

    tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
    tr.tcrdist2(processes = cpu,
                metric = 'nw',
                reduce = True,
                dump = False,
                save = True,
                dest = "default_archive", 
                dest_tar_name = "default_archive.tar.gz")        


* **processes** - number of parallel processes to use
* **metric** - (default : `nw`) specifies the Needleman Wunsch metric
* **reduce** - (default : True) If True, reduces distances from float to int16 data types to reduce file size
* **dump** - (default : False) If True, dumps pairwise distance matrics for CDR1, CDR2, and CDR2.5 after used 
* **save** - (default : False) If True, saves results to dest folder, and .tar.gz
* **dest** - Name of folder for saving results
* **dest_tar_name** - Name of archive for saving results as a .tar.gz
* **verbose** - If True, provide information to sys.stdout on each step.

.. code-block::

    Deduplicating your TCRrep.cell_df to make TCRrep.clone_df.
    Computing pairwise matrices for multiple Complementarity Determining Regions (CDRs):.
            Computing pairwise matrices for cdrs within the alpha-chain using the nw metric.
            Computing pairwise matrices for cdrs within the beta-chain using the nw metric.
    Calculating composite tcrdistance measures:
            Single chain pairwise tcrdistances are in attribute : TCRrep.pw_alpha
            Single chain pairwise tcrdistances are in attribute : TCRrep.pw_beta
            Combined pairwise tcrdistances are in attribute     : TCRrep.pw_tcrdist
            CDR specific tcrdistances are in attributes, e.g.,  : TCRrep.cdr3_b_aa_pw
    Reducing File Size: `reduce` argumment set to True:
            Reducing : paired_tcrdist to int16.
            Reducing : cdr3_a_aa_pw to int16.
            Reducing : cdr2_a_aa_pw to int16.
            Reducing : cdr1_a_aa_pw to int16.
            Reducing : pmhc_a_aa_pw to int16.
            Reducing : cdr3_b_aa_pw to int16.
            Reducing : cdr2_b_aa_pw to int16.
            Reducing : cdr1_b_aa_pw to int16.
            Reducing : pmhc_b_aa_pw to int16.
            Reducing : pw_alpha to int16.
            Reducing : pw_beta to int16.
            Reducing : pw_tcrdist to int16.
    Archiving your TCRrep using Zipdist2 (save = True)
            Saving paired_tcrdist to .csv : default_archive/paired_tcrdist.csv
            Saving paired_tcrdist to .npy : default_archive/paired_tcrdist.npy
            Saving cdr3_a_aa_pw to .csv : default_archive/cdr3_a_aa_pw.csv
            Saving cdr3_a_aa_pw to .npy : default_archive/cdr3_a_aa_pw.npy
            Saving cdr2_a_aa_pw to .csv : default_archive/cdr2_a_aa_pw.csv
            Saving cdr2_a_aa_pw to .npy : default_archive/cdr2_a_aa_pw.npy
            Saving cdr1_a_aa_pw to .csv : default_archive/cdr1_a_aa_pw.csv
            Saving cdr1_a_aa_pw to .npy : default_archive/cdr1_a_aa_pw.npy
            Saving pmhc_a_aa_pw to .csv : default_archive/pmhc_a_aa_pw.csv
            Saving pmhc_a_aa_pw to .npy : default_archive/pmhc_a_aa_pw.npy
            Saving cdr3_b_aa_pw to .csv : default_archive/cdr3_b_aa_pw.csv
            Saving cdr3_b_aa_pw to .npy : default_archive/cdr3_b_aa_pw.npy
            Saving cdr2_b_aa_pw to .csv : default_archive/cdr2_b_aa_pw.csv
            Saving cdr2_b_aa_pw to .npy : default_archive/cdr2_b_aa_pw.npy
            Saving cdr1_b_aa_pw to .csv : default_archive/cdr1_b_aa_pw.csv
            Saving cdr1_b_aa_pw to .npy : default_archive/cdr1_b_aa_pw.npy
            Saving pmhc_b_aa_pw to .csv : default_archive/pmhc_b_aa_pw.csv
            Saving pmhc_b_aa_pw to .npy : default_archive/pmhc_b_aa_pw.npy
            Saving pw_alpha to .csv : default_archive/pw_alpha.csv
            Saving pw_alpha to .npy : default_archive/pw_alpha.npy
            Saving pw_beta to .csv : default_archive/pw_beta.csv
            Saving pw_beta to .npy : default_archive/pw_beta.npy
            Saving pw_tcrdist to .csv : default_archive/pw_tcrdist.csv
            Saving pw_tcrdist to .npy : default_archive/pw_tcrdist.npy
            Saving cell_df to .csv : default_archive/cell_df.csv
            Saving cell_df to .feather : default_archive/cell_df.feather
            Saving clone_df to .csv : default_archive/clone_df.csv
            Saving clone_df to .feather : default_archive/clone_df.feather
            Saving JSON with complex attribute definitions : default_archive/complex_attributes.json
            Saving JSON with simple attribute definitions : default_archive/simple_attributes.json
            Combining saved files in : [default_archive.tar.gz].
            Archiving your TCRrep using Zipdist2 in [default_archive.tar.gz]
            Archiving your TCRrep using Zipdist2 in [default_archive.tar.gz]
    TCRrep.tcrdist2() COMPLETED SUCCESSFULLY, see the docs for Analysis steps!



TCRdist Legacy Distance
-----------------------

Calculating a legacy TCRdist using methods from Dash et al. 2017, can be achieved
with a few prepatory steps and :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_alpha_beta()`.


.. code-block:: python

    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv('dash.csv')
    tr = TCRrep(cell_df=df, chains=['alpha','beta'] organism='mouse')
    tr.infer_cdrs_from_v_gene(chain='alpha',  imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain='beta',  imgt_aligned=True)
    tr.index_cols = ['clone_id', 'subject', 'epitope',
                    'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                    'cdr3_a_aa', 'cdr3_b_aa',
                    'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                    'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                    'cdr3_b_nucseq', 'cdr3_a_nucseq',
                    'va_gene', 'vb_gene',
                    'ja_gene', 'jb_gene']
    tr.deduplicate()
    tr._tcrdist_legacy_method_alpha_beta()
    tr.pw_tcrdist



Other methods are available for single chain and paired chain analysis:

* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_alpha_beta`
* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_alpha`
* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_beta`
* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_gamma_delta`
* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_gamma`
* :py:meth:`tcrdist.repertoire.TCRrep._tcrdist_legacy_method_delta`
