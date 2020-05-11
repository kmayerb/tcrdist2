Probability of Generation
-------------------------

Estimate probability of generation of all the cdr3 sequences in a TCRrep. (tcrdist2 incorporates OLGA for generation probability estimates).

Currently probability of generations can be performed for human TCR beta- and alpha-chain cdr3s and mouse beta-chain cdr3s, using the 
TCRrep class method :py:func:`tcrdist.repertoire.TCRrep.infer_olga_aa_cdr3_pgens`.

Example 
=======

.. code-block:: python

    tr.infer_olga_aa_cdr3_pgens(chain = "beta")


.. tip::

    See complete commands for loading data in the :ref:`CompleteExample`.

.. code-block:: python

    tr = TCRrep(cell_df = tcrdist2_df, organism = "mouse")
    tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)
    tr.index_cols = ['clone_id', 'subject', 'epitope',
                    'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                    'cdr3_a_aa', 'cdr3_b_aa',
                    'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                    'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                    'cdr3_b_nucseq', 'cdr3_a_nucseq',
                    'va_countreps', 'ja_countreps',
                    'vb_countreps', 'jb_countreps',
                    'va_gene', 'vb_gene',
                    'ja_gene', 'jb_gene']
    tr.deduplicate()
    tr.clone_df['cdr3_b_pgen'] = tr.infer_olga_aa_cdr3_pgens(chain = "beta")
    tr.clone_df[['v_b_gene', 'j_b_gene', 'cdr3_b_aa', 'cdr3_b_pgen']].head()


.. code-block::

        v_b_gene    j_b_gene       cdr3_b_aa   cdr3_b_pgen
    0  TRBV19*03  TRBJ2-7*01     CASSSGGEQYF  1.688137e-05
    1   TRBV4*01  TRBJ2-1*01  CASSYPDIYAEQFF  4.584401e-08
    2  TRBV29*01  TRBJ1-5*01     CASSEGEAPLF  9.738778e-09
    3  TRBV29*01  TRBJ2-7*01     CASAGGDEQYF  0.000000e+00
    4  TRBV29*01  TRBJ2-1*01     CASSPDGEQFF  6.143096e-10



Further options for inferring a pgen using OLGA can be referenced :py:func:`tcrdist.repertoire.TCRrep.infer_olga_aa_cdr3_pgens`.

Docs for olga wrapper class can be found :py:class:`tcrdist.pgen.OlgaModel`.

.. tip::
                              
    The default makes use of Igor models which are stored in `tcrdist/olga/default_models/`; however, the `chain_folder` option 
    could be used to estimate pGen with a Igor new model. A new model should be placed into source code location `tcrdist/olga/default_models/new_model/`,
    and could be referenced using non-default `chain_folder` argument. 

    .. code-block:: python

        tr.clone_df['cdr3_b_pgen'] = tr.infer_olga_aa_cdr3_pgens(chain = "beta", chain_folder = "new_model" )



Citation
========

OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs

Zachary Sethna, Yuval Elhanati, Curtis G Callan, Aleksandra M Walczak, Thierry Mora

`Bioinformatics (2019) <https://doi.org/10.1093/bioinformatics/btz035>`_


High-throughput immune repertoire analysis with IGoR 

Marcou, Quentin, Thierry Mora, and Aleksandra M. Walczak.

`Nature communications (2018) <https://www.nature.com/articles/s41467-018-02832-w>`_


Public Method
=============

TCRrep class method :py:func:`tcrdist.repertoire.TCRrep.infer_olga_aa_cdr3_pgens`.


Private OlgaModel Class
=======================

The OlgaModel class is used to wrap olga code: :py:class:`tcrdist.pgen.OlgaModel`.
